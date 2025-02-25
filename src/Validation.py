# ========================================
# Import modules
# ========================================
import pyvista as pv
import numpy as np
import vtkmodules.util.numpy_support as vnp
import math
import PySimpleGUI as sg
import sys
import pymeshfix
import os

# ========================================
# Functions
# ========================================

def slice_z(mesh, n_slices=300):
    """
    English:
    This function receives as entry a mesh and generates slices along the z 
    axis. This is useful for visualizing and analyzing the 3D mesh along the 
    z-axis, which can be particularly helpful in understanding the structure 
    and geometry of the mesh.

    Spanish:
    Esta función recibe como entrada una malla y genera cortes en torno al 
    eje z. Esto puede ser muy útil para visualizar y analizar la malla 3D en
    torno al eje z, algo particularmente útil para comprender la estructura 
    y la geometría.
    """
    # Create vector
    vec = np.array([0, 0, 1.0])

    # Normalize the vector
    normal = vec / np.linalg.norm(vec)

    # Make points along that vector for the extent of the slices
    z_length = np.max(mesh.bounds[5]) - np.min(mesh.bounds[4])
    a = mesh.center + normal * z_length / 2.01
    b = mesh.center - normal * z_length / 2.01

    # Define the line/points for the slices
    line = pv.Line(a, b, n_slices)


    # Generate all of the slices
    slices = pv.MultiBlock()
    for point in line.points:
        slices.append(mesh.slice(normal=normal, origin=point))
    
    return slices, line

def distance2points_3d(point1, point2):
    """
    English:
    This function can be used to determine the distance between any two 
    points in a 3D space, which is a fundamental operation in many 3D 
    applications.

    Spanish:
    Esta función puede ser utilizada para determinar la distancia entre
    cualquier par de puntos en un espacio 3D, que es una operación 
    fundamental en muchas aplicaciones 3D.
    """
    # Distance between the two points using the Euclidean distance formula
    distance = math.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2 + (point2[2] - point1[2])**2)
    
    return distance

def diam_calc(slice, center):
    """
    English:
    This function is used to calculate the diameter of a slice, which is 
    the distance between the two furthest points in the slice. It also 
    calculates the minimum distance to the center of the slice.

    Spanish:
    Esta función se utiliza para calcular el diámetro de una sección, que
    es la distancia entre los dos puntos más alejados de la sección. También
    calcula la distancia mínima entre la sección y el centro.
    """
    max_dist = 0
    min_dist = 100
    point_a = None
    point_b = None

    for i in range(len(slice)):
        center_dist = distance2points_3d(slice[i], center)
        if center_dist < min_dist: 
            min_dist = center_dist
        for j in range(i+1, len(slice)):
            distancia = distance2points_3d((slice[i]), (slice[j]))
            if distancia > max_dist:
                max_dist = distancia
                point_a = (slice[i])
                point_b = (slice[j])
    
    return point_a, point_b, min_dist, max_dist

def find_blocks(list):
    """
    English:
    This function takes a list of positions as input and returns a list of 
    tuples. Each tuple represents a group of consecutive positions that are 
    within a certain distance (in this case, less than 5 units apart).

    Spanish:
    Esta función toma una lista de posiciones como entrada y devuelve una
    lista de tuplas. Cada tupla representa un grupo de posiciones consecutivas
    que están a una distancia menor a 5 unidades.
    """
    blocks = []
    current_block = [list[0]]

    for i in range(1, len(list)):
        if list[i] - list[i - 1] < 5:
            current_block.append(list[i])
        else:
            blocks.append(tuple(current_block))
            current_block = [list[i]]

    blocks.append(tuple(current_block))
    return blocks

def load_stl(filename):
    """
    English:
    This function is used to load an STL file and extract its vertices and 
    faces for further processing in the given Python script.

    Spanish:
    Esta función se utiliza para cargar un archivo STL y extraer sus
    vértices y caras para su procesamiento en el script dado.
    """
    # Open the STL file in read mode
    with open(filename, 'r') as file:
        lines = file.readlines()
    vertices = []
    faces = []

    for line in lines:
        parts = line.split()
        if len(parts) > 0:
            if parts[0] == 'vertex':
                vertex = list(map(float, parts[1:4]))
                vertices.append(vertex)
            elif parts[0] == 'endloop':
                face = [vertices[-3], vertices[-2], vertices[-1]]
                faces.append(face)

    return np.array(vertices), np.array(faces)

def save_stl(filename, faces):
    """
    English:
    This function is responsible for writing the vertex and face information
    of a 3D mesh to a new STL file.

    Spanish:
    Esta función se encarga de escribir la información de vértices y caras
    de una malla 3D a un nuevo archivo STL.
    """
    with open(filename, 'w') as file:
        file.write("solid filtered_model\n")
        for face in faces:
            file.write("  facet normal 0 0 0\n")
            file.write("    outer loop\n")
            for vertex in face:
                file.write("      vertex {} {} {}\n".format(vertex[0], vertex[1], vertex[2]))
            file.write("    endloop\n")
            file.write("  endfacet\n")
        file.write("endsolid filtered_model\n")

def add2vectors(v1, v2):
    """
    English:
    This function takes two vectors as input and returns their sum as a 
    unit vector. 

    Spanish:
    Esta función toma dos vectores como entrada y devuelve su suma como
    un vector unitario.
    """
    # Convert list to numpy array
    v1_np = np.array(v1)
    v2_np = np.array(v2)
    
    # Normalice vectors
    v1_norm = v1_np / np.linalg.norm(v1_np)
    v2_norm = v2_np / np.linalg.norm(v2_np)
    
    result = v1_norm + v2_norm
    
    # Normalice result to obtain a unit vector
    unit_vector = result / np.linalg.norm(result)
    
    # Return result as a list of Python
    return unit_vector.tolist() 

def centerline_aprox(slices, n=400):
    """
    English:
    This function takes a list of slices as input and generates a centerline
    approximation using a spline. This function is used to create a smooth 
    and continuous representation of the centerline of the 3D mesh slices, 
    which can then be used for further analysis or visualization.

    Spanish:
    Esta función toma una lista de secciones como entrada y genera una
    aproximación de la línea central usando una spline. Esta función se
    utiliza para crear una representación suave y continua de la línea
    central de las secciones de la malla 3D, que puede luego ser utilizada
    para la análisis o la visualización.
    """
    centers = []
    for slice in slices:
        center = slice.center
        centers.append(center)
        
    # Generate spline with 300 interpolation points
    centerline = pv.Spline(centers, n)

    # Centerline length estimation  
    centerline_length = 0
    spline_puntos = centerline.GetPoints() 
    for i in range(n-1):
        centerline_length += distance2points_3d(spline_puntos.GetPoint(i), spline_puntos.GetPoint(i+1))
    
    # Finally, the function returns the generated spline (centerline) 
    # and the estimated centerline length.
    return centerline, centerline_length

def ortogonal_slices(centerline, mesh, n=400):
    """
    English:
    This function takes a centerline as input and generates a list of
    orthogonal slices based on the centerline. This forthogonal slices 
    can be used for further analysis or visualization.

    Spanish:
    Esta función toma una línea central como entrada y genera una lista
    de secciones ortogonales basadas en la línea central. Estas secciones
    ortogonales pueden ser utilizadas para posterior análisis o visualización.
    """
    orto_slices = pv.MultiBlock()
    for j in range(n-1):

        # Alternative option:
        #   p0 = centerline.GetPoint(0)
        #   normal = np.array([p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]])
        #   normales.append(normal)

        point1 = centerline.GetPoint(j)
        point2 = centerline.GetPoint(j+1)
            
        normal_centerline = np.array([point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]])
        vector = add2vectors(normal_centerline, np.array([0,0,-1]))
        orto_slices.append(mesh.slice(normal=vector, origin=point1))

    return orto_slices


# ========================================
# Graphical User Interface (GUI)
# ========================================
"""
GUI: 
File Chooser with clearable history
    
    1. Copy and paste a filename into the combo element
    2. Use the last used item which will be visible when you create the window
    3. Choose an item from the list of previously used items
    4. Browse for a new name
    
To clear the list of previous entries, click the "Clear History" button.
The history is stored in a json file using the PySimpleGUI User Settings APIs   
"""

layout = [[sg.Combo(sorted(sg.user_settings_get_entry('-file-', [])), default_value=sg.user_settings_get_entry('-last file-', ''), size=(50, 1), key='-FILE-'), sg.FileBrowse(), sg.B('Clear History')],
          [sg.Button('Ok', bind_return_key=True),  sg.Button('Cancel')]]

window = sg.Window('File Chooser With History', layout)

while True:
    event, values = window.read()

    if event in (sg.WIN_CLOSED, 'Cancel'):
        sys.exit()
    if event == 'Ok':
        # If OK, then need to add the filename to the list of files and also set as the last used filename
        sg.user_settings_set_entry('-file-', list(set(sg.user_settings_get_entry('-file-', []) + [values['-FILE-'], ])))
        sg.user_settings_set_entry('-last file-', values['-FILE-'])
        entry_file = values['-FILE-']
        break
    elif event == 'Clear History':
        sg.user_settings_set_entry('-file-', [])
        sg.user_settings_set_entry('-last file-', '')
        window['-FILE-'].update(values=[], value='')
        sys.exit()

window.close()

# ========================================
# File Processing 
# ========================================

"""
SLICE ALONG VECTOR PARALLEL TO Z AXIS
"""

mesh = pv.read(entry_file)
mesh_slices, line = slice_z(mesh)

p = pv.Plotter()
p.add_mesh(line, line_width=4)
p.add_mesh(mesh, opacity=0.7)
p.show()

# Uncomment the following lines to visualize the slices:
p = pv.Plotter()
p.add_mesh(mesh.outline(), color="k")
p.add_mesh(mesh_slices, opacity=0.75)
p.show()

"""
CENTERLINE APROXIMATION 
"""
mesh_centerline, mesh_centerline_length = centerline_aprox(mesh_slices)

# Uncomment the following lines to visualize the centerlines:
p = pv.Plotter()
p.add_mesh(mesh_slices, opacity=0.75)
p.add_mesh(mesh_centerline, line_width=4, color='k')
p.show()

"""
DIAMETER CALCULUS & PLOTTING 
"""
diameters = []
orto_centers = []

for i, slice in enumerate(mesh_slices):
    np_points = vnp.vtk_to_numpy(slice.GetPoints().GetData())
    point_a, point_b, min_dist, max_dist = diam_calc(np_points, mesh_centerline.GetPoint(i))
    diameters.append(max_dist)
                
"""
PLOT RED/BLUE
"""  
limit_diameter = 30
p = pv.Plotter(shape=(1, 2))

# Limit diameter:
p.subplot(0,0)
for i, slice in enumerate(mesh_slices):  
    if diameters[i] >= limit_diameter:
        p.add_mesh(slice, opacity=0.75, color = 'red')
    else:
        p.add_mesh(slice, opacity=0.75, color = 'blue')
p.add_mesh(mesh_centerline, color='k', line_width=2)

# Adaptative diameter:
p.subplot(0,1)
dini = diameters[:5]
initial_diameter = sum(dini)/len(dini)
adaptative_diameter = initial_diameter*1.5

for i, slice in enumerate(mesh_slices):
    if diameters[i] >= limit_diameter or diameters[i] >= adaptative_diameter: 
        p.add_mesh(slice, opacity=0.75, color='red')
    else:
        p.add_mesh(slice, opacity=0.75, color='blue')
p.add_mesh(mesh_centerline, color='k', line_width=2)   
p.show()

"""
MAX DIAMETER PLOT
"""
# Uncomment the following lines to visualize the max diameter:
p = pv.Plotter()
p.add_mesh(mesh, style='wireframe', color='w', opacity=0.5)
p.add_mesh(mesh_slices[diameters.index(max(diameters))], color='purple', line_width=4)
p.show()

"""
ANEURYSM DETECTION
"""
xx = 10
aneurysm = True
positions = []
zone_of_interest = []

# Save the positions of the aneurysm:
for i in range(len(diameters)):
    if diameters[i]>limit_diameter or diameters[i] >= initial_diameter*1.50:
        positions.append(i)
    
if len(positions) < xx:
    aneurysm = False

# Consecutive positions:
else:
    blocks = find_blocks(positions)
    for i in range(len(blocks)):
        if len(blocks[i]) > xx:
            zone_of_interest.extend(blocks[i])

if len(zone_of_interest) == 0:
    aneurysm = False

if aneurysm == True:
    bounda = (mesh_slices[min(zone_of_interest)].bounds)
    bounda2 = (mesh_slices[max(zone_of_interest)].bounds)

    aneurysm_mesh = mesh.clip(normal='z', origin=[0,0,bounda2[4]], invert=False)
    aneurysm_mesh = aneurysm_mesh.clip(normal='z', origin=[0,0,bounda[5]], invert=True)


    # Uncomment the following lines to visualize the aneurysm mesh:
    p = pv.Plotter()
    p.add_mesh(aneurysm_mesh, color='r')
    p.add_mesh(mesh, style='wireframe')
    clip_plane = pv.Plane((bounda[1], bounda[3], bounda[5]), (0,0,1), i_size=50, j_size=50) 
    clip_plane2 = pv.Plane((bounda2[1], bounda2[3], bounda2[4]), (0,0,1), i_size=50, j_size=50)
    p.add_mesh(clip_plane)
    p.add_mesh(clip_plane2)
    p.show()

    fixer = pymeshfix.MeshFix(aneurysm_mesh.triangulate())
    fixer.repair()
    aneurysm = fixer.mesh
    print('Aneurysm = True')
    print(aneurysm.volume)

    # Uncomment the following lines to visualize the closed aneurysm mesh:
    p = pv.Plotter()
    p.add_mesh(aneurysm)
    p.show()

"""
RESULTS
"""
# print(diameters)
# precision = (len(diameters)*6/sum(diameters))*100
# print(precision)
print(max(diameters))
print(mesh_centerline_length)
print(mesh_centerline.length)




