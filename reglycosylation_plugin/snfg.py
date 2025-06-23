from pymol import cmd, cgo
#from chempy import cpv # No longer used directly here
import numpy as np
import math
import warnings # To suppress potential numpy warnings

# SNFG Colors (constants moved to top for clarity)
SNFG_WHITE =      (255/255, 255/255, 255/255)
SNFG_BLUE =       (  0/255, 144/255, 188/255)
SNFG_GREEN =      (  0/255, 166/255,  81/255)
SNFG_YELLOW =     (255/255, 212/255,   0/255)
SNFG_LIGHT_BLUE = (143/255, 204/255, 233/255)
SNFG_PINK =       (246/255, 158/255, 161/255)
SNFG_PURPLE =     (165/255,  67/255, 153/255)
SNFG_BROWN =      (161/255, 122/255,  77/255)
SNFG_ORANGE =     (244/255, 121/255,  32/255)
SNFG_RED =        (237/255,  28/255,  36/255)

# Axis colors for debugging
AXIS_X_COLOR = (1.0, 0.0, 0.0) # Red
AXIS_Y_COLOR = (0.0, 1.0, 0.0) # Green
AXIS_Z_COLOR = (0.0, 0.0, 1.0) # Blue

def normalize_vector(v):
    """Normalizes a numpy vector."""
    norm = np.linalg.norm(v)
    if norm < 1e-9: # Use a small threshold instead of == 0 for float comparison
        # Return a default vector or raise an error if normalization fails
        warnings.warn(f"Attempted to normalize a near-zero vector: {v}. Returning default Z=[0,0,1].")
        # Return an arbitrary basis vector, e.g., Z-axis
        return np.array([0., 0., 1.])
    return v / norm

def get_residue_orientation(selection):
    """
    Get the geometric centers and the local coordinate system (X, Y, Z axes)
    for specified residues in the selection.
    Z-axis: points from the residue midpoint towards C4.
    Y-axis: perpendicular to Z, lies within the 'averaged' plane of C1-C5.
    X-axis: completes the right-handed coordinate system.
    """
    orientations = {}
    model = cmd.get_model(selection)
    residue_atoms = {}

    # Collect coordinates for relevant atoms per residue
    required_atoms = {'C1', 'C2', 'C3', 'C4', 'C5', 'O5'} # Atoms for center calculation
    plane_atoms = {'C1', 'C2', 'C3', 'C4', 'C5'}         # Atoms defining the plane

    for atom in model.atom:
        # Only consider atoms relevant for center, plane, or Z-axis definition
        if atom.name in required_atoms:
            res_key = (atom.chain, atom.resi, atom.resn)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = {}
            residue_atoms[res_key][atom.name] = np.array(atom.coord)

    # Calculate orientation for each residue
    for res_key, atoms in residue_atoms.items():
        resn = res_key[2]

        # Default axes in case calculation fails
        default_x = np.array([1., 0., 0.])
        default_y = np.array([0., 1., 0.])
        default_z = np.array([0., 0., 1.])

        # --- 1. Calculate Center ---
        center_coords = [atoms[name] for name in required_atoms if name in atoms]
        if len(center_coords) < 3: # Need at least 3 points for a meaningful center
             print(f"Skipping residue {res_key}: Not enough atoms ({len(center_coords)} found) among {required_atoms} for center calculation.")
             # Store default orientation with a best-guess center if possible
             available_coords = list(atoms.values())
             center = np.mean(available_coords, axis=0) if available_coords else np.array([0.,0.,0.])
             orientations[res_key] = (resn, center, default_x, default_y, default_z)
             continue

        center = np.mean(center_coords, axis=0)

        # --- 2. Calculate Z-axis (Center -> C4) ---
        if 'C4' not in atoms:
            print(f"Skipping residue {res_key}: Missing C4 atom for Z-axis. Using default orientation.")
            orientations[res_key] = (resn, center, default_x, default_y, default_z)
            continue # Skip further calculations if C4 is missing

        coord_C4 = atoms['C4']
        # Check if center and C4 are coincident
        vec_c4 = coord_C4 - center
        if np.linalg.norm(vec_c4) < 1e-6:
             print(f"Warning for residue {res_key}: Center and C4 are nearly coincident. Using default Z-axis.")
             z_axis = default_z
        else:
             z_axis = normalize_vector(vec_c4)


        # --- 3. Calculate Plane Normal (PCA on C1-C5) ---
        plane_coords_list = [atoms[name] for name in plane_atoms if name in atoms]
        plane_normal = default_z # Default plane normal if calculation fails

        if len(plane_coords_list) < 3:
            print(f"Warning for residue {res_key}: Need at least 3 atoms (C1-C5) to define a plane. Using default plane.")
            # Proceed using default plane normal, Z axis already calculated
        else:
            plane_coords = np.array(plane_coords_list)
            plane_center = np.mean(plane_coords, axis=0)
            centered_coords = plane_coords - plane_center

            # Use SVD for robustness, especially when N <= D
            try:
                # We need at least 2 dimensions for SVD on coordinates
                if centered_coords.shape[1] < 2:
                     raise np.linalg.LinAlgError("Need at least 2 dimensions for SVD")
                # Ensure matrix rank allows finding a unique normal
                if np.linalg.matrix_rank(centered_coords) < min(centered_coords.shape)-1 :
                    # Points might be collinear or coincident
                    print(f"Warning for residue {res_key}: C1-C5 atoms may be collinear or coincident. Using default plane normal.")
                    # plane_normal remains default_z
                else:
                    _, _, vh = np.linalg.svd(centered_coords) # vh contains right singular vectors as rows
                    plane_normal = normalize_vector(vh[-1]) # Normal is the last row vector
            except np.linalg.LinAlgError as e:
                print(f"Warning: SVD failed for {res_key} plane calculation ({e}). Using default plane normal.")
                plane_normal = default_z # Keep default normal if SVD fails


        # --- 4. Define X and Y axes based on Z and Plane Normal ---
        # Ensure plane_normal is not parallel to z_axis
        if np.abs(np.dot(plane_normal, z_axis)) > 0.995: # Increased tolerance slightly
            # Z-axis is (almost) parallel to the plane normal (i.e., perpendicular to the C1-C5 plane).
            # We need *any* vector roughly "in the plane" (i.e., perpendicular to Z) to define Y.
            # Let's try using C1->C2 vector projected onto the plane orthogonal to Z.
            print(f"Note: Z-axis for {res_key} is nearly parallel to C1-C5 plane normal. Defining Y axis arbitrarily orthogonal to Z.")

            temp_vec = np.array([1.0, 0.0, 0.0]) # Start with global X
            # If Z is aligned with global X, use global Y instead
            if np.abs(np.dot(temp_vec, z_axis)) > 0.99:
                temp_vec = np.array([0.0, 1.0, 0.0])

            # Define Y as orthogonal to Z using the cross product
            y_axis = normalize_vector(np.cross(z_axis, temp_vec))
            # Define X using the right-hand rule (X = Y x Z)
            x_axis = normalize_vector(np.cross(y_axis, z_axis))

        else:
            # Normal case: Z and plane_normal are not parallel.
            # Calculate X-axis (orthogonal to Z and Plane Normal) -> lies roughly in the C1-C5 plane
            x_axis = normalize_vector(np.cross(plane_normal, z_axis)) # Swapped order z, normal -> normal, z for RH system? Check! Z=X x Y, Y=Z x X, X=Y x Z. If Z is 'up' and normal is 'out', X=normal x Z points 'right'. OK.
            # Calculate Y-axis (completes right-handed system: Y = Z x X)
            y_axis = normalize_vector(np.cross(z_axis, x_axis))


        # Final check for NaN values (can happen with coincident points etc.)
        if np.isnan(x_axis).any() or np.isnan(y_axis).any() or np.isnan(z_axis).any():
            print(f"ERROR for residue {res_key}: Resulting axes contain NaN. Using default orientation.")
            orientations[res_key] = (resn, center, default_x, default_y, default_z)
        else:
            # Store the results: residue name, center, and the basis vectors
            orientations[res_key] = (resn, center, x_axis, y_axis, z_axis)

    # Return list of tuples: (resn, center, x_axis, y_axis, z_axis)
    return list(orientations.values())


def get_transformation_matrix(x_axis, y_axis, z_axis):
    """
    Creates a 3x3 transformation matrix from the basis vectors.
    This matrix transforms vectors from the local coordinate system (where the
    object is defined aligned with standard axes) to the global system.
    Columns of the matrix are the target basis vectors in the global frame.
    """
    # Ensure axes are normalized (should be already, but good practice)
    x_axis = normalize_vector(x_axis)
    y_axis = normalize_vector(y_axis)
    z_axis = normalize_vector(z_axis)

    # Check orthogonality and handedness (optional, for debugging)
    # det = np.linalg.det(np.array([x_axis, y_axis, z_axis]).T)
    # if abs(det - 1.0) > 0.01:
    #     print(f"Warning: Axes may not form a perfect right-handed system. Det={det:.3f}")
    #     print(f"  X.Y={np.dot(x_axis, y_axis):.3f}, Y.Z={np.dot(y_axis, z_axis):.3f}, Z.X={np.dot(z_axis, x_axis):.3f}")

    return np.array([x_axis, y_axis, z_axis]).T


def transform_vertices(vertices, center, rotation_matrix):
    """
    Transforms vertices defined relative to the origin in a standard orientation
    to the target orientation and position.
    """
    transformed_vertices = []
    center = np.array(center)
    for v in vertices:
        v = np.array(v)                   # Vertex relative to origin (local coords)
        v_rot = np.dot(rotation_matrix, v) # Rotate to global orientation
        transformed_vertices.append(v_rot + center) # Translate to final position
    return transformed_vertices

def transform_normals(normals, rotation_matrix):
    """ Transforms normal vectors using the rotation matrix. """
    transformed_normals = []
    for n in normals:
        n = np.array(n)
        n_rot = np.dot(rotation_matrix, n) # Apply rotation
        # Normal vectors should remain normalized, handle potential zero vectors
        transformed_normals.append(normalize_vector(n_rot))
    return transformed_normals

# --- CGO Shape Functions (Unchanged from previous version) ---

def cgo_star(x, y, z, r, color, x_axis, y_axis, z_axis):
    center = np.array([x, y, z])
    rotation_matrix = get_transformation_matrix(x_axis, y_axis, z_axis)

    r_inner = 0.45 * r
    star_cgo = []

    # Angles for the 5 points (outer and inner) in the XY plane (local coords)
    angles = [math.radians(i * 72) for i in range(5)]

    # Define the outer and inner points relative to the origin (0,0,0)
    outer_points_local = []
    inner_points_local = []
    for i, angle in enumerate(angles):
        outer_points_local.append([r * math.cos(angle), 0, r * math.sin(angle)])
        inner_angle = angle + math.radians(36)
        inner_points_local.append([r_inner * math.cos(inner_angle), 0, r_inner * math.sin(inner_angle)])

    # Define top/bottom points along the local Z axis
    volume_points_local = [[0, 0.6 * r, 0], [0, -0.6 * r,0]]

    # Transform points to global coordinates
    outer_points = transform_vertices(outer_points_local, center, rotation_matrix)
    inner_points = transform_vertices(inner_points_local, center, rotation_matrix)
    volume_points = transform_vertices(volume_points_local, center, rotation_matrix)

    # Calculate face normals (approximation - better normals would require knowing triangle vertices)
    # For now, use the transformed Z axis for the flat faces
    normal_front = transform_normals([[0,1,0]], rotation_matrix)[0]
    normal_back = transform_normals([[0,-1,0]], rotation_matrix)[0]

    # CGO object generation
    star_cgo.extend([cgo.COLOR, *color])

    # # Front face
    star_cgo.extend([cgo.BEGIN, cgo.TRIANGLE_FAN, cgo.NORMAL, *normal_front, cgo.VERTEX, *volume_points[0]])

    for i in range(5):
        star_cgo.extend([cgo.VERTEX, *outer_points[i], cgo.VERTEX, *inner_points[i]])
    star_cgo.extend([cgo.VERTEX, *outer_points[0]]) # Close the fan
    star_cgo.extend([cgo.END])

    # # Back face
    star_cgo.extend([cgo.BEGIN, cgo.TRIANGLE_FAN, cgo.NORMAL, *normal_back, cgo.VERTEX, *volume_points[1]])

    # # Add vertices in reverse order for correct winding/normal
    for i in range(5):
        idx = i
        # idx = (4-i) # index from 4 down to 0
        star_cgo.extend([cgo.VERTEX, *outer_points[idx], cgo.VERTEX, *inner_points[idx]])
    star_cgo.extend([cgo.VERTEX, *outer_points[0]]) # Close the fan
    star_cgo.extend([cgo.END])

    return star_cgo


def cgo_sphere(x, y, z, size, color):
    # Spheres are rotationally invariant, no orientation needed
    return [
        cgo.COLOR, *color,
        cgo.SPHERE, float(x), float(y), float(z), float(size)*1.3
    ]

def cgo_cone(x, y, z, r, color, x_axis, y_axis, z_axis):
    center = np.array([x, y, z])
    rotation_matrix = get_transformation_matrix(x_axis, y_axis, z_axis)

    # Define diamond vertices relative to the origin (0,0,0) along local axes
    vertices_local = [
        [0, 0, r],    # Top (+Z local)
        [0, 0, -r]   # Bottom (-Z local)
    ]

    vertices = transform_vertices(vertices_local, center, rotation_matrix)

    return [
        cgo.COLOR, *color,
        cgo.CONE, *vertices[0], *vertices[1], 0, r, *color, *color, 1, 1
    ]

def cgo_diamond(x, y, z, r, color, x_axis, y_axis, z_axis):
    center = np.array([x, y, z])
    rotation_matrix = get_transformation_matrix(x_axis, y_axis, z_axis)

    # Define diamond vertices relative to the origin (0,0,0) along local axes
    vertices_local = [
        [0, 0, r],    # Top (+Z local)
        [0, 0, -r],   # Bottom (-Z local)
        [0, r, 0],    # "+Y local" vertex
        [0, -r, 0],   # "-Y local" vertex
        [r, 0, 0],    # "+X local" vertex
        [-r, 0, 0]    # "-X local" vertex
    ]
    vertices = transform_vertices(vertices_local, center, rotation_matrix)

    # Normals for the 8 faces (approximate based on average of vertices)
    # Normal calculations should ideally be based on actual triangle vertices and cross products
    # Using transformed basis vectors scaled slightly is a simpler approximation:
    normals_local = [ # Approximations for the 8 faces pointing outwards
         normalize_vector(np.array([ 1,  1,  1])), normalize_vector(np.array([-1,  1,  1])),
         normalize_vector(np.array([-1, -1,  1])), normalize_vector(np.array([ 1, -1,  1])),
         normalize_vector(np.array([ 1,  1, -1])), normalize_vector(np.array([-1,  1, -1])),
         normalize_vector(np.array([-1, -1, -1])), normalize_vector(np.array([ 1, -1, -1])),
    ]
    normals = transform_normals(normals_local, rotation_matrix)

    diamond_cgo = [cgo.COLOR, *color]
    diamond_cgo.extend([cgo.BEGIN, cgo.TRIANGLES])
    # Face 1 (Top, +Y, +X local) - Vertices 0, 2, 4
    diamond_cgo.extend([cgo.NORMAL, *normals[0], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[2], cgo.VERTEX, *vertices[4]])
    # Face 2 (Top, +Y, -X local) - Vertices 0, 5, 2 (Reordered for winding)
    diamond_cgo.extend([cgo.NORMAL, *normals[1], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[5], cgo.VERTEX, *vertices[2]])
    # Face 3 (Top, -Y, -X local) - Vertices 0, 3, 5
    diamond_cgo.extend([cgo.NORMAL, *normals[2], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[3], cgo.VERTEX, *vertices[5]])
    # Face 4 (Top, -Y, +X local) - Vertices 0, 4, 3 (Reordered for winding)
    diamond_cgo.extend([cgo.NORMAL, *normals[3], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[4], cgo.VERTEX, *vertices[3]])
    # Face 5 (Bottom, +Y, +X local) - Vertices 1, 4, 2 (Reordered for winding)
    diamond_cgo.extend([cgo.NORMAL, *normals[4], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[4], cgo.VERTEX, *vertices[2]])
    # Face 6 (Bottom, +Y, -X local) - Vertices 1, 2, 5
    diamond_cgo.extend([cgo.NORMAL, *normals[5], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[2], cgo.VERTEX, *vertices[5]])
    # Face 7 (Bottom, -Y, -X local) - Vertices 1, 5, 3 (Reordered for winding)
    diamond_cgo.extend([cgo.NORMAL, *normals[6], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[5], cgo.VERTEX, *vertices[3]])
    # Face 8 (Bottom, -Y, +X local) - Vertices 1, 3, 4
    diamond_cgo.extend([cgo.NORMAL, *normals[7], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[3], cgo.VERTEX, *vertices[4]])
    diamond_cgo.extend([cgo.END])
    return diamond_cgo

def cgo_half_diamond(x, y, z, r, color, color_2, x_axis, y_axis, z_axis):
    center = np.array([x, y, z])
    rotation_matrix = get_transformation_matrix(x_axis, y_axis, z_axis)

    # Define diamond vertices relative to the origin (0,0,0) along local axes
    vertices_local = [
        [0, 0, r],    # Top (+Z local)
        [0, 0, -r],   # Bottom (-Z local)
        [0, r, 0],    # "+Y local" vertex
        [0, -r, 0],   # "-Y local" vertex
        [r, 0, 0],    # "+X local" vertex
        [-r, 0, 0]    # "-X local" vertex
    ]
    vertices = transform_vertices(vertices_local, center, rotation_matrix)

    # Normals for the 8 faces (approximate based on average of vertices)
    # Normal calculations should ideally be based on actual triangle vertices and cross products
    # Using transformed basis vectors scaled slightly is a simpler approximation:
    normals_local = [ # Approximations for the 8 faces pointing outwards
         normalize_vector(np.array([ 1,  1,  1])), normalize_vector(np.array([-1,  1,  1])),
         normalize_vector(np.array([-1, -1,  1])), normalize_vector(np.array([ 1, -1,  1])),
         normalize_vector(np.array([ 1,  1, -1])), normalize_vector(np.array([-1,  1, -1])),
         normalize_vector(np.array([-1, -1, -1])), normalize_vector(np.array([ 1, -1, -1])),
    ]
    normals = transform_normals(normals_local, rotation_matrix)
    diamond_cgo = []
    diamond_cgo.extend([cgo.BEGIN, cgo.TRIANGLES])
    # Face 1 (Top, +Y, +X local) - Vertices 0, 2, 4
    diamond_cgo.extend([cgo.COLOR, *color_2])
    diamond_cgo.extend([cgo.NORMAL, *normals[0], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[2], cgo.VERTEX, *vertices[4]])
    # Face 2 (Top, +Y, -X local) - Vertices 0, 5, 2 (Reordered for winding)
    diamond_cgo.extend([cgo.COLOR, *color])
    diamond_cgo.extend([cgo.NORMAL, *normals[1], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[5], cgo.VERTEX, *vertices[2]])
    # Face 3 (Top, -Y, -X local) - Vertices 0, 3, 5
    diamond_cgo.extend([cgo.COLOR, *color])
    diamond_cgo.extend([cgo.NORMAL, *normals[2], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[3], cgo.VERTEX, *vertices[5]])
    # Face 4 (Top, -Y, +X local) - Vertices 0, 4, 3 (Reordered for winding)
    diamond_cgo.extend([cgo.COLOR, *color_2])
    diamond_cgo.extend([cgo.NORMAL, *normals[3], cgo.VERTEX, *vertices[0], cgo.VERTEX, *vertices[4], cgo.VERTEX, *vertices[3]])
    # Face 5 (Bottom, +Y, +X local) - Vertices 1, 4, 2 (Reordered for winding)
    diamond_cgo.extend([cgo.COLOR, *color_2])
    diamond_cgo.extend([cgo.NORMAL, *normals[4], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[4], cgo.VERTEX, *vertices[2]])
    # Face 6 (Bottom, +Y, -X local) - Vertices 1, 2, 5
    diamond_cgo.extend([cgo.COLOR, *color])
    diamond_cgo.extend([cgo.NORMAL, *normals[5], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[2], cgo.VERTEX, *vertices[5]])
    # Face 7 (Bottom, -Y, -X local) - Vertices 1, 5, 3 (Reordered for winding)
    diamond_cgo.extend([cgo.COLOR, *color])
    diamond_cgo.extend([cgo.NORMAL, *normals[6], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[5], cgo.VERTEX, *vertices[3]])
    # Face 8 (Bottom, -Y, +X local) - Vertices 1, 3, 4
    diamond_cgo.extend([cgo.COLOR, *color_2])
    diamond_cgo.extend([cgo.NORMAL, *normals[7], cgo.VERTEX, *vertices[1], cgo.VERTEX, *vertices[3], cgo.VERTEX, *vertices[4]])
    diamond_cgo.extend([cgo.END])
    return diamond_cgo

def cgo_cube(x, y, z, r, color, x_axis, y_axis, z_axis):
    center = np.array([x, y, z])
    rotation_matrix = get_transformation_matrix(x_axis, y_axis, z_axis)
    r_half = r # Assume r is half the side length

    # Vertices relative to origin, ordered for TRIANGLE_STRIP per face
    vertices_local = [
        # Face +X (Normal: +X) - Strip order: (0, 2, 1, 3) -> Indices in list: 0, 1, 2, 3
        [+r_half, +r_half, -r_half], [+r_half, -r_half, -r_half], [+r_half, +r_half, +r_half], [+r_half, -r_half, +r_half],
        # Face -X (Normal: -X) - Strip order: (4, 6, 5, 7) -> Indices in list: 4, 5, 6, 7
        [-r_half, +r_half, +r_half], [-r_half, -r_half, +r_half], [-r_half, +r_half, -r_half], [-r_half, -r_half, -r_half],
        # Face +Y (Normal: +Y) - Strip order: (8, 10, 9, 11) -> Indices in list: 8, 9, 10, 11
        [-r_half, +r_half, +r_half], [+r_half, +r_half, +r_half], [-r_half, +r_half, -r_half], [+r_half, +r_half, -r_half],
        # Face -Y (Normal: -Y) - Strip order: (12, 14, 13, 15) -> Indices in list: 12, 13, 14, 15
        [-r_half, -r_half, -r_half], [+r_half, -r_half, -r_half], [-r_half, -r_half, +r_half], [+r_half, -r_half, +r_half],
        # Face +Z (Normal: +Z) - Strip order: (16, 18, 17, 19) -> Indices in list: 16, 17, 18, 19
        [-r_half, +r_half, +r_half], [-r_half, -r_half, +r_half], [+r_half, +r_half, +r_half], [+r_half, -r_half, +r_half],
        # Face -Z (Normal: -Z) - Strip order: (20, 22, 21, 23) -> Indices in list: 20, 21, 22, 23
        [-r_half, -r_half, -r_half], [-r_half, +r_half, -r_half], [+r_half, -r_half, -r_half], [+r_half, +r_half, -r_half]
    ]
    strip_orders = [ # Indices relative to the start of each 4-vertex block
        [0, 2, 1, 3], # +X
        [0, 2, 1, 3], # -X (Indices 4, 6, 5, 7)
        [0, 2, 1, 3], # +Y (Indices 8, 10, 9, 11)
        [0, 2, 1, 3], # -Y (Indices 12, 14, 13, 15)
        [0, 2, 1, 3], # +Z (Indices 16, 18, 17, 19)
        [0, 2, 1, 3]  # -Z (Indices 20, 22, 21, 23)
    ]

    # Define normals aligned with local axes
    normals_local = [
        [ 1,  0,  0], [-1,  0,  0], [ 0,  1,  0],
        [ 0, -1,  0], [ 0,  0,  1], [ 0,  0, -1]
    ]

    # Transform vertices and normals
    vertices = transform_vertices(vertices_local, center, rotation_matrix)
    normals = transform_normals(normals_local, rotation_matrix)

    cube_cgo = [cgo.COLOR, *color]
    for i in range(6): # Iterate through 6 faces
        cube_cgo.extend([cgo.BEGIN, cgo.TRIANGLE_STRIP])
        cube_cgo.extend([cgo.NORMAL, *normals[i]])
        # Get the 4 vertices for this face
        face_vertex_start_index = i * 4
        # Apply the strip order for this face
        for strip_index in strip_orders[i]:
            vertex_index = face_vertex_start_index + strip_index
            cube_cgo.extend([cgo.VERTEX, *vertices[vertex_index]])
        cube_cgo.append(cgo.END)

    return cube_cgo


# --- Main Function ---

def render_snfg(selection='all', transparency=0.3, scale=0.5, debug_axes=False):
    """
    Render glycans using the 3D-SNFG color and shape scheme with updated orientation.

    Parameters:
    selection (str): Atom selection for glycans (default: 'all')
    transparency (float): Transparency of sugar representations (Note: CGO transparency is complex)
    scale (float): Size of sugar representations
    debug_axes (bool): If True, draw X(red), Y(green), Z(blue) axes for each residue.
    """
    # Define SNFG color and shape scheme for common sugars
    snfg_shapes = {
        'GLC': (SNFG_BLUE, 'sphere'), # Glucose
        'MAL': (SNFG_BLUE, 'sphere'), # Alias
        'BGC': (SNFG_BLUE, 'sphere'), # Alias
        'MAN': (SNFG_GREEN, 'sphere'), # Mannose
        'BMA': (SNFG_GREEN, 'sphere'), # Alias
        'GAL': (SNFG_YELLOW, 'sphere'),
        'GLA': (SNFG_YELLOW, 'sphere'), # Alias
        'FUC': (SNFG_RED, 'cone'), # Official SNFG is Red Triangle
        'FUL': (SNFG_RED, 'cone'), # Alias
        'XYL': (SNFG_ORANGE, 'star'), # Official SNFG is Orange Star
        'XYP': (SNFG_ORANGE, 'star'), # Alias
        'XYS': (SNFG_ORANGE, 'star'), # Alias
        'ARA': (SNFG_GREEN, 'star'), # Arabinose 
        'AHR': (SNFG_GREEN, 'star'), # Alias
        'RIB': (SNFG_PINK, 'star'), # Ribose
        'NAG': (SNFG_BLUE, 'cube'), # GlcNAc
        'GLCNAC': (SNFG_BLUE, 'cube'), # Alias
        '4YS': (SNFG_BLUE, 'cube'), # Alias
        'SGN': (SNFG_BLUE, 'cube'), # Alias
        'BGLN': (SNFG_BLUE, 'cube'), # Alias
        'NDG': (SNFG_BLUE, 'cube'), # Alias
        'NGA': (SNFG_YELLOW, 'cube'), # GalNAc
        'GALNAC': (SNFG_YELLOW, 'cube'), # Alias
        'A2G' : (SNFG_YELLOW, 'cube'), # Alias
        'MANNA': (SNFG_GREEN, 'cube'), # ManNAc
        # 'SIA': (SNFG_RED, 'diamond'), # Sia
        'NEU5AC': (SNFG_PURPLE, 'diamond'), # Neu5Ac / Sia
        'SIA': (SNFG_PURPLE, 'diamond'), # Alias
        'NEU5GC': (SNFG_LIGHT_BLUE, 'diamond'), # Neu5Gc
        'NGC': (SNFG_LIGHT_BLUE, 'diamond'), # Alias
        'KDN': (SNFG_GREEN, 'diamond'), # KDN (Usually green diamond)
        'ADA': (SNFG_YELLOW, 'half_diamond'), # GalA
        'GLCA': (SNFG_BLUE, 'half_diamond'), # GlcA
        'GCU' : (SNFG_BLUE, 'half_diamond'), # Alias
        'BDP': (SNFG_BLUE, 'half_diamond'), # Alias
        'IDOA': (SNFG_BROWN, 'half_diamond_reverse'), # IdoA (Official SNFG is brown half-circle up) - Using diamond
        'IDS': (SNFG_BROWN, 'half_diamond_reverse'), # Alias
        'IDR': (SNFG_BROWN, 'half_diamond_reverse'), # Alias
        'API': (SNFG_PINK, 'diamond'), # Api (Official SNFG is pink cross rectangle) - Using diamond
        # Add other required sugars here following the pattern
    }

    # Select residues matching the keys in snfg_shapes
    target_resn = '+'.join(snfg_shapes.keys())
    # Sanitize selection string for CGO object name
    safe_selection_name = "".join(c if c.isalnum() else "_" for c in selection)
    if not safe_selection_name: safe_selection_name = "selection" # Handle empty/weird selections

    full_selection = f"({selection}) and resn {target_resn}"
    print(f"Applying render_snfg to selection: '{full_selection}'")

    # Get orientation data for the selected residues
    orientation_data = get_residue_orientation(full_selection)

    if not orientation_data:
        print("[render_snfg] No matching residues found or orientation calculation failed for all.")
        return

    cgo_objects = []
    debug_cgo_objects = [] # Separate list for debug axes
    shape_scale_factor = 4.0 * scale # Base size factor, adjust as needed

    # Define axis properties for debugging
    axis_length = shape_scale_factor * 1.5 # Length relative to shape size
    axis_radius = scale * 0.1          # Radius for cylinder representation

    print(f"[render_snfg] Processing {len(orientation_data)} residues...")
    processed_count = 0
    for (resn, center, x_axis, y_axis, z_axis) in orientation_data:
        if resn not in snfg_shapes:
            # This check might be redundant due to initial selection, but safe
            print(f"Skipping residue {resn} at {center} as it's not in snfg_shapes dict (should not happen).")
            continue

        color, shape = snfg_shapes[resn]
        x, y, z = center

        # Create CGO based on shape and orientation
        cgo_obj = []
        try:
            if shape == 'sphere':
                # Sphere size might need different scaling
                cgo_obj = cgo_sphere(x, y, z, 1.8 * scale, color) # Slightly larger base scale for spheres
            elif shape == 'star':
                cgo_obj = cgo_star(x, y, z, shape_scale_factor, color, x_axis, y_axis, z_axis)
            elif shape == 'cone':
                cgo_obj = cgo_cone(x, y, z, shape_scale_factor * 0.75, color, x_axis, y_axis, z_axis)
            elif shape == 'cube':
                cgo_obj = cgo_cube(x, y, z, shape_scale_factor * 0.55, color, x_axis, y_axis, z_axis) # Cube radius is half-side length, adjusted scale
            elif shape == 'diamond':
                cgo_obj = cgo_diamond(x, y, z, shape_scale_factor * 0.75, color, x_axis, y_axis, z_axis) # Diamond radius needs tuning
            elif shape == 'half_diamond':
                cgo_obj = cgo_half_diamond(x, y, z, shape_scale_factor * 0.75, color, SNFG_WHITE, x_axis, y_axis, z_axis) # Diamond radius needs tuning
            elif shape == 'half_diamond_reverse':
                cgo_obj = cgo_half_diamond(x, y, z, shape_scale_factor * 0.75, SNFG_WHITE, color, x_axis, y_axis, z_axis) # Diamond radius needs tuning
            else:
                 print(f"Warning: Shape '{shape}' for residue {resn} is not implemented. Using sphere.")
                 cgo_obj = cgo_sphere(x, y, z, 1.8 * scale, color)

            if cgo_obj: # Ensure something was generated
                 cgo_objects.extend(cgo_obj)
                 processed_count += 1
            else:
                 print(f"Warning: No CGO generated for {resn} at {center} with shape {shape}.")


            # --- Add Debug Axes (if requested) ---
            if debug_axes:
                # Ensure axes are valid before drawing
                if not (np.isnan(x_axis).any() or np.isnan(y_axis).any() or np.isnan(z_axis).any()):
                    # X-axis (Red)
                    end_x = center + x_axis * axis_length
                    debug_cgo_objects.extend([
                        cgo.CYLINDER, *center, *end_x, axis_radius, *AXIS_X_COLOR, *AXIS_X_COLOR # Same color start/end
                    ])
                    # Y-axis (Green)
                    end_y = center + y_axis * axis_length
                    debug_cgo_objects.extend([
                        cgo.CYLINDER, *center, *end_y, axis_radius, *AXIS_Y_COLOR, *AXIS_Y_COLOR
                    ])
                    # Z-axis (Blue)
                    end_z = center + z_axis * axis_length
                    debug_cgo_objects.extend([
                        cgo.CYLINDER, *center, *end_z, axis_radius, *AXIS_Z_COLOR, *AXIS_Z_COLOR
                    ])
                else:
                    print(f"Skipping debug axes for {resn} at {center} due to NaN axis values.")


        except Exception as e:
            print(f"ERROR generating CGO for {resn} at {center}: {e}")
            import traceback
            traceback.print_exc() # Print detailed traceback for debugging


    if not cgo_objects:
         print("[render_snfg] No valid CGO objects were generated for the main shapes.")
         # Still load debug axes if requested and generated
         if debug_axes and debug_cgo_objects:
            cgo_name = f"snfg_{safe_selection_name}" # Need a base name even if main obj is empty
            debug_cgo_name = f"{cgo_name}_axes"
            cmd.load_cgo(debug_cgo_objects, debug_cgo_name)
            print(f"[render_snfg] Debug axes loaded as CGO object: {debug_cgo_name}")
         return

    # Load the combined CGO object for shapes
    cgo_name = f"snfg_{safe_selection_name}"
    cmd.delete(cgo_name) # Delete existing object with the same name first
    cmd.load_cgo(cgo_objects, cgo_name)
    print(f"[render_snfg] Processed {processed_count} residues. Main CGO object loaded as: {cgo_name}")

    # --- Load Debug Axes CGO (if requested) ---
    if debug_axes and debug_cgo_objects:
        debug_cgo_name = f"{cgo_name}_axes"
        cmd.delete(debug_cgo_name) # Delete existing debug axes object first
        cmd.load_cgo(debug_cgo_objects, debug_cgo_name)
        print(f"[render_snfg] Debug axes loaded as CGO object: {debug_cgo_name}")
    elif debug_axes:
        print("[render_snfg] Debug axes requested, but none were generated (check for errors above).")


    # Note on transparency: CGO objects loaded this way generally don't respect cmd.set('cgo_transparency').
    # True transparency needs setting ALPHA within the CGO definition itself, which adds complexity.
    print(f"[render_snfg] Note: Transparency setting ({transparency}) is not directly applied to CGO shapes.")


# Register in PyMOL
cmd.extend("render_snfg", render_snfg)

# Example Usage (run in PyMOL console after loading this script):
# fetch 1rvz, async=0
# remove solvent
# render_snfg selection=resn NAG+MAN, scale=0.6, debug_axes=True
# show cartoon
# hide lines, (resn NAG+MAN) # Hide underlying atoms if desired
# # To hide axes later:
# # hide everything, snfg_resn_NAG_MAN_axes
# # Or delete:
# # delete snfg_resn_NAG_MAN_axes