# -*- coding: utf-8 -*-

"""
PyMOL Plugin: html_interface.py

Description:
A PyMOL plugin that opens a local static HTML file (index.html)
located alongside this script. It uses QtWebEngine to display the
HTML and QtWebChannel to establish two-way communication between
the Python environment (PyMOL) and JavaScript running in the HTML page.

Requirements:
- PyMOL built with PyQt5/PySide2 support.
- PyQt5/PySide2 installed in PyMOL's Python environment.
- The QtWebEngine component must be available (usually installed via
  `pip install PyQtWebEngine` or `conda install pyqtwebengine` if not
  bundled with PyMOL's Qt).
- An 'index.html' file in the same directory as this script.
- The 'qwebchannel.js' file (from your PyQt5/PySide2 installation)
  in the same directory as this script.

Usage:
1. Place this script, __init__.py, index.html, and qwebchannel.js
   in a directory (e.g., "html_interface_plugin") inside your
   PyMOL plugins directory.
2. Restart PyMOL.
3. Go to Plugin -> HTML Interface Plugin.
4. A window will open displaying index.html, allowing interaction.
"""

import os
import sys
import json
from pymol import cmd, plugins

# --- Attempt to set OpenGL context sharing attribute EARLY ---
QT_CORE_LOADED = False
try:
    from PyQt5 import QtCore as QtCore_Base
    QtCore_Base.QCoreApplication.setAttribute(QtCore_Base.Qt.AA_ShareOpenGLContexts)
    print("HTML Interface Plugin: Set Qt.AA_ShareOpenGLContexts using PyQt5.")
    QT_CORE_LOADED = True
except ImportError:
    try:
        from PySide2 import QtCore as QtCore_Base
        QtCore_Base.QCoreApplication.setAttribute(QtCore_Base.Qt.AA_ShareOpenGLContexts)
        print("HTML Interface Plugin: Set Qt.AA_ShareOpenGLContexts using PySide2.")
        QT_CORE_LOADED = True
    except ImportError:
        print("HTML Interface Plugin Error: Could not import QtCore from PyQt5 or PySide2 "
              "to set AA_ShareOpenGLContexts.")

# --- Import PyMOL's Qt wrapper and attempt WebEngine/WebChannel import ---
try:
    from pymol.Qt import QtWidgets, QtCore
    QT_WIDGETS_LOADED = True
except ImportError:
    print("HTML Interface Plugin Error: Could not import QtWidgets/QtCore from pymol.Qt.")
    QT_WIDGETS_LOADED = False
    QtWidgets = None
    QtCore = None

# --- Attempt to import WebEngine and WebChannel components ---
QT_WEB_AVAILABLE = False # Combined check for Engine + Channel
QtWebEngineWidgets = None
QtWebChannel = None

if QT_CORE_LOADED and QT_WIDGETS_LOADED:
    try:
        # Most common case: PyQt5
        from PyQt5 import QtWebEngineWidgets, QtWebChannel
        print("HTML Interface Plugin: Using PyQt5 for WebEngine & WebChannel.")
        QT_WEB_AVAILABLE = True
    except ImportError:
        try:
            # Less common: PySide2
            from PySide2 import QtWebEngineWidgets, QtWebChannel
            print("HTML Interface Plugin: Using PySide2 for WebEngine & WebChannel.")
            QT_WEB_AVAILABLE = True
        except ImportError:
            # Neither found
            print("HTML Interface Plugin Error: Could not import QtWebEngineWidgets "
                  "and/or QtWebChannel from PyQt5 or PySide2.")
            if QT_CORE_LOADED and QtWidgets:
                 QtWidgets.QMessageBox.critical(None, "Plugin Error",
                    "HTML Interface Plugin requires QtWebEngine & QtWebChannel support "
                    "(from PyQt5 or PySide2), which was not found.\n\n"
                    "Please ensure the necessary components are installed in PyMOL's Python environment.\n"
                    "The plugin cannot load.")
            QtWebEngineWidgets = None
            QtWebChannel = None # Ensure it's None if import failed
else:
    print("HTML Interface Plugin: Skipping WebEngine/WebChannel check due to missing core Qt components.")


# --- Global storage for window instances ---
INSTANCES = {}

# --- The Bridge Object ---
# This object's slots are exposed to JavaScript via QWebChannel.
# Only define if Web components are available.
if QT_WEB_AVAILABLE:
    class PyMolBridge(QtCore.QObject):
        """Object exposed to JavaScript for communication."""

        # Example signal (optional): Python -> JS notification
        # dataReady = QtCore.pyqtSignal(str) # Emits JSON string

        def __init__(self, web_view_widget, parent=None):
            super().__init__(parent)
            # Store reference to the parent window/widget to access webView if needed
            self.web_view_widget = web_view_widget

        # --- Slots callable FROM JavaScript ---

        @QtCore.pyqtSlot(str, result=str)
        def executePymolCommand(self, command_string):
            """Executes a PyMOL command sent from JavaScript."""
            print(f"HTML Interface: Received command: {command_string}")
            if not command_string:
                return json.dumps({"status": "error", "message": "Empty command received."})
            try:
                # Use cmd.do for simplicity. Consider cmd.async_do for long commands.
                # Security Note: Be cautious executing arbitrary strings if HTML source
                # is not fully trusted.
                cmd.do(command_string)
                print(f"HTML Interface: Executed command successfully.")
                # Return status to JS callback as JSON string
                return json.dumps({"status": "success", "command": command_string})
            except Exception as e:
                error_message = f"Error executing command '{command_string}': {e}"
                print(f"HTML Interface: {error_message}")
                # Return error details to JS callback
                return json.dumps({"status": "error", "message": error_message})

        @QtCore.pyqtSlot(str, result=str)
        def getDataFromPyMol(self, request_type):
            """Retrieves data from PyMOL based on request_type."""
            print(f"HTML Interface: Received data request: {request_type}")
            data = {}
            try:
                if request_type == "molecule_names":
                    data = {"names": cmd.get_names("objects")}
                elif request_type == "current_view":
                    data = {"view": cmd.get_view()}
                # Add more request types as needed
                else:
                    data = {"error": f"Unknown request type: {request_type}"}

                # Return data as a JSON string to the JavaScript callback
                return json.dumps(data)
            except Exception as e:
                error_message = f"Error getting data for '{request_type}': {e}"
                print(f"HTML Interface: {error_message}")
                return json.dumps({"error": error_message})

        # --- Method callable FROM Python to send data TO JavaScript ---
        def send_data_to_js(self, data):
            """Sends arbitrary data (must be JSON serializable) to the JS function 'updateFromPyMol'."""
            if not self.web_view_widget or not self.web_view_widget.webView:
                print("HTML Interface Error: WebView not available to send data.")
                return

            try:
                # Ensure data is JSON serializable before sending
                json_data = json.dumps(data)
                # Escape the JSON string for safe injection into JS and call the target function
                # Make sure the JS function `updateFromPyMol` exists in your index.html
                js_command = f"if (typeof updateFromPyMol === 'function') {{ updateFromPyMol({json_data}); }} else {{ console.error('JS function updateFromPyMol not found'); }}"
                print(f"HTML Interface: Sending data to JS function updateFromPyMol: {json_data[:100]}...") # Log truncated data
                self.web_view_widget.webView.page().runJavaScript(js_command)
            except TypeError as te:
                 print(f"HTML Interface Error: Data is not JSON serializable: {te}")
            except Exception as e:
                print(f"HTML Interface Error: Failed to send data to JS: {e}")


# --- The Main Plugin Window ---
# Only define if Web components are available.
if QT_WEB_AVAILABLE and QtWebEngineWidgets and QtWebChannel:
    class HtmlInterfaceWindow(QtWidgets.QWidget):
        """Main window hosting the QWebEngineView for local index.html."""
        def __init__(self, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Glyco.me SugarBuilder")
            self.resize(800, 600) # Default size

            # Re-verify availability at instantiation
            if not QT_WEB_AVAILABLE:
                 self.setup_error_ui("QtWebEngine/WebChannel components failed to load or became unavailable.")
                 return

            # --- Find Local HTML and JS Files ---
            self.script_dir = os.path.dirname(os.path.realpath(__file__))
            self.html_path = os.path.join(self.script_dir, "index.html")
            self.qwebchannel_js_path = os.path.join(self.script_dir, "qwebchannel.js")

            # Check if required files exist
            if not os.path.exists(self.html_path):
                error_msg = f"Cannot find index.html at:\n{self.html_path}"
                print(f"HTML Interface Plugin Error: {error_msg}")
                self.setup_error_ui(error_msg)
                #QtWidgets.QMessageBox.critical(self, "Plugin Error", error_msg) # Alternative error display
                return # Stop initialization

            if not os.path.exists(self.qwebchannel_js_path):
                 # Communication will fail without this, so warn prominently
                 error_msg = f"Cannot find qwebchannel.js at:\n{self.qwebchannel_js_path}\n\n" \
                             "Communication between PyMOL and the HTML page will likely fail.\n" \
                             "Please copy qwebchannel.js from your PyQt5/PySide2 installation " \
                             "into the plugin directory."
                 print(f"HTML Interface Plugin Warning: {error_msg}")
                 QtWidgets.QMessageBox.warning(self, "Missing File", error_msg)
                 # Continue loading, but communication might not work

            # --- Setup WebView ---
            try:
                self.webView = QtWebEngineWidgets.QWebEngineView()
                self.webView.page().profile().setHttpUserAgent("PyMOL_HTML_Interface_Plugin/1.1")
            except Exception as e:
                print(f"HTML Interface Error: Failed to create QWebEngineView: {e}")
                self.setup_error_ui(f"Failed to create Web View: {e}")
                return

            # --- Setup WebChannel ---
            self.bridge = None # Initialize to None
            self.channel = None
            try:
                self.channel = QtWebChannel.QWebChannel(self.webView.page())
                # Create the bridge object, pass self (widget) so bridge can call JS back
                self.bridge = PyMolBridge(web_view_widget=self)
                # Register the bridge object with the name JavaScript will use
                self.channel.registerObject("pyMolBridge", self.bridge)
                # Make the channel available to the JavaScript context
                self.webView.page().setWebChannel(self.channel)
                print("HTML Interface Plugin: QWebChannel setup successful.")
            except Exception as e:
                 # Should ideally not happen if imports succeeded, but good practice
                 print(f"HTML Interface Error: Failed to set up QWebChannel: {e}")
                 self.setup_error_ui(f"Failed to set up Web Channel: {e}")
                 # Ensure bridge/channel are None if setup fails
                 self.bridge = None
                 self.channel = None
                 # Optionally return here if channel is absolutely critical


            # --- Layout ---
            self.layout = QtWidgets.QVBoxLayout(self)
            self.layout.addWidget(self.webView)
            self.layout.setContentsMargins(0, 0, 0, 0)

            # --- Example Button (Python -> JS) ---
            # Add a button in the Python GUI to trigger sending data TO the webpage
            self.send_to_js_button = QtWidgets.QPushButton("Send Object Names to HTML")
            self.send_to_js_button.clicked.connect(self.send_object_names_to_js)
             # Only enable the button if the bridge was successfully created
            if self.bridge:
                 self.layout.addWidget(self.send_to_js_button)
            else:
                 print("HTML Interface Plugin: Python->JS button disabled (bridge unavailable).")


            # --- Load Local HTML File ---
            local_url = QtCore.QUrl.fromLocalFile(self.html_path)
            print(f"HTML Interface Plugin: Loading local file: {local_url.toString()}")

            # Connect signals for loading feedback (optional but useful)
            self.webView.loadFinished.connect(self._on_load_finished)
            self.webView.loadProgress.connect(self._on_load_progress)
            self.webView.loadStarted.connect(self._on_load_started)

            # Set the URL in the webView
            self.webView.setUrl(local_url)

            self.webView.settings().setAttribute(QtWebEngineWidgets.QWebEngineSettings.WebAttribute.JavascriptEnabled, True)
            self.webView.settings().setAttribute(QtWebEngineWidgets.QWebEngineSettings.WebAttribute.AllowRunningInsecureContent, True)
            self.webView.settings().setAttribute(QtWebEngineWidgets.QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, True)
            self.webView.settings().setAttribute(QtWebEngineWidgets.QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True)


        def setup_error_ui(self, message):
             """Helper to display an error message if setup fails."""
             self.webView = None # Ensure webview is None
             self.channel = None
             self.bridge = None
             # Check if layout exists before trying to add to it
             if not hasattr(self, 'layout') or self.layout is None:
                  self.layout = QtWidgets.QVBoxLayout(self) # Create layout if missing

             # Clear existing widgets if any
             while self.layout.count():
                 item = self.layout.takeAt(0)
                 widget = item.widget()
                 if widget is not None:
                     widget.deleteLater()

             error_label = QtWidgets.QLabel(f"Plugin Error:\n{message}\n\nPlease check PyMOL console/log for details.")
             error_label.setWordWrap(True)
             error_label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
             self.layout.addWidget(error_label)
             self.layout.addStretch(1)
             print(f"HTML Interface Plugin: Setup failed with error: {message}")


        # --- Handlers for Load Status (Optional) ---
        def _on_load_started(self):
            if not self.webView: return
            print(f"HTML Interface Plugin: Load started for {self.webView.url().toString()}...")
            self.setWindowTitle("Loading HTML Interface...")

        def _on_load_progress(self, progress):
            if not self.webView: return
            self.setWindowTitle(f"Loading HTML Interface ({progress}%)...")

        def _on_load_finished(self, ok):
            if not self.webView: return
            url_str = self.webView.url().toString()
            print(f"HTML Interface Plugin: Load finished for {url_str}. Success: {ok}")
            if ok:
                self.setWindowTitle(f"Glyco.me SugarBuilder - {self.webView.title()}")
                 # Check if the channel seems active (basic check)
                if self.channel and self.bridge:
                    print("HTML Interface Plugin: Page loaded. QWebChannel bridge should be available to JavaScript.")
                elif self.channel:
                    print("HTML Interface Plugin Warning: Page loaded, but bridge object seems missing.")
                else:
                    print("HTML Interface Plugin Warning: Page loaded, but QWebChannel setup may have failed earlier.")
            else:
                 self.setWindowTitle(f"Load Failed - {url_str}")


        # --- Example Method for Python -> JS Button ---
        def send_object_names_to_js(self):
            """Example method triggered by the Python button."""
            if not self.bridge:
                print("HTML Interface Plugin: Cannot send data to JS, bridge not available.")
                QtWidgets.QMessageBox.warning(self, "Bridge Error", "Communication bridge to HTML is not active.")
                return

            try:
                names = cmd.get_names("objects")
                data_to_send = {
                    "type": "object_list",
                    "source": "python_button",
                    "objects": names,
                    "count": len(names)
                }
                # Use the bridge's method to send data
                self.bridge.send_data_to_js(data_to_send)
            except Exception as e:
                error_msg = f"Could not get object names to send: {e}"
                print(f"HTML Interface Error: {error_msg}")
                # Optionally send an error message to JS as well
                # self.bridge.send_data_to_js({"error": error_msg})
                QtWidgets.QMessageBox.warning(self, "PyMOL Error", error_msg)

        def closeEvent(self, event):
            """Clean up when the window is closed."""
            print("HTML Interface Plugin: Closing window.")
            instance_key = 'html_interface_main'
            INSTANCES.pop(instance_key, None) # Remove from global dict

            # Safely disconnect signals
            if self.webView:
                 try:
                     self.webView.loadStarted.disconnect(self._on_load_started)
                     self.webView.loadProgress.disconnect(self._on_load_progress)
                     self.webView.loadFinished.disconnect(self._on_load_finished)
                 except (TypeError, RuntimeError): pass # Ignore if already disconnected

            # Clean up WebChannel and WebView
            if hasattr(self, 'webView') and self.webView:
                 if hasattr(self, 'channel') and self.channel:
                     # Disconnect page from channel first
                     self.webView.page().setWebChannel(None)
                 # Stop loading and schedule web view for deletion
                 self.webView.stop()
                 self.webView.deleteLater()
                 self.webView = None

            if hasattr(self, 'channel') and self.channel:
                 if hasattr(self, 'bridge') and self.bridge:
                     # Deregister the specific object
                     self.channel.deregisterObject(self.bridge)
                     self.bridge = None # Clear reference
                 # Schedule channel for deletion
                 self.channel.deleteLater()
                 self.channel = None

            super().closeEvent(event)


# --- PyMOL Plugin Initialization ---
def __init_plugin__(app=None):
    """Plugin initialization function called by PyMOL."""
    # Check if prerequisites loaded successfully
    if not QT_WIDGETS_LOADED:
         print("HTML Interface Plugin: Not loaded (Qt Widgets missing).")
         return
    if not QT_WEB_AVAILABLE: # Check the combined flag for Engine+Channel
        print("HTML Interface Plugin: Not loaded (QtWebEngine/QtWebChannel missing or failed).")
        return

    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Draw a sugar', run_plugin_gui)

# Function to launch the plugin GUI
def run_plugin_gui():
    """Creates and shows the main plugin window. Manages single instance."""
    # Re-check availability
    if not QT_WEB_AVAILABLE or not QtWidgets:
         if QtWidgets:
             QtWidgets.QMessageBox.critical(None, "Plugin Error",
                "HTML Interface Plugin cannot run (QtWebEngine/QtWebChannel missing or failed).")
         else:
             print("HTML Interface Plugin cannot run: Critical Qt components missing.")
         return

    instance_key = 'html_interface_main'
    if instance_key in INSTANCES and INSTANCES[instance_key].isVisible():
        INSTANCES[instance_key].raise_()
        INSTANCES[instance_key].activateWindow()
    else:
        # Create with no parent to avoid Tkinter issues
        parent = None
        # Ensure the class was defined
        if 'HtmlInterfaceWindow' in globals():
             try:
                 window_instance = HtmlInterfaceWindow(parent)
                 # Check if initialization failed internally (e.g., setup_error_ui called)
                 if hasattr(window_instance, 'layout') and window_instance.layout is not None and window_instance.webView is not None:
                     INSTANCES[instance_key] = window_instance
                     INSTANCES[instance_key].show()
                 else:
                     if window_instance: window_instance.deleteLater() # Clean up failed instance
                     print("HTML Interface Plugin: Window setup failed during initialization.")

             except Exception as e:
                 print(f"HTML Interface Plugin Error: Failed to create HtmlInterfaceWindow: {e}")
                 QtWidgets.QMessageBox.critical(None, "Plugin Error",
                    f"HTML Interface Plugin failed to create window:\n{e}")
        else:
             print("HTML Interface Plugin Error: HtmlInterfaceWindow class not defined (Web components missing).")
             QtWidgets.QMessageBox.critical(None, "Plugin Error",
                "HTML Interface Plugin failed (HtmlInterfaceWindow class missing).")

# Optional cleanup (might not be called reliably)
# def __del_plugin__():
#    instance_key = 'html_interface_main'
#    if instance_key in INSTANCES:
#        try: INSTANCES[instance_key].close()
#        except Exception as e: print(f"HTML Interface Plugin: Error during cleanup: {e}")
#    print("HTML Interface Plugin unloaded.")