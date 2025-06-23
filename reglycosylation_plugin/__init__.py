# Initialize the plugin when PyMOL starts
# This checks if the module is being imported by PyMOL's plugin loader

# Use a try-except block for compatibility
try:
    # PyMOL 2.x specific plugin loading
    import pymol.plugins
    def __init_plugin__(app=None):
        from . import reglycosylation_interface # Import the main plugin code (relative import)
        reglycosylation_interface.__init_plugin__(app)

    # Optional: Define cleanup if needed (might not be reliably called)
    # def __del_plugin__():
    #     from . import html_interface
    #     html_interface.__del_plugin__()

except ImportError:
    # Fallback or logic for older PyMOL versions if necessary
    # For simple cases, just importing might be enough if reglycosylation_interface.py
    # directly calls addmenuitemqt on import (less clean).
    # The __init_plugin__ structure is preferred.
    print("ReGlyco Plugin: Could not use __init_plugin__, attempting direct load (may fail).")
    from . import reglycosylation_interface
    # If html_interface.py called addmenuitemqt directly at module level:
    # pass
    # If not, you might need to call an init function here. The pattern above is better.