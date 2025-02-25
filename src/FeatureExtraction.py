# ========================================
# Import modules
# ========================================
import pyvista as pv
import numpy as np
import vtkmodules.util.numpy_support as vnp
import math
import pymeshfix
import PySimpleGUI as sg
import sys
import json
import os
from itertools import combinations


# ========================================
# Functions
# ========================================
def split_model(input_file):
    """
    English:
    This function is used to divide the solids in an STL file into separate 
    files, which is a common preprocessing step in 3D modeling and computer 
    graphics. 

    Spanish:
    Esta función se utiliza para dividir los sólidos en un archivo STL en 
    archivos individuales, que es un paso común en el preprocesamiento de 
    modelado 3D y gráficos computacionales.
    """
    # Open the input file in read mode
    with open(input_file, 'r') as entry:
        content = entry.read()

        # Divides the content in lines
        lines = content.splitlines()

        # Iterating through the lines of an STL file, identifying 
        # the end of each solid (marked by the keyword "endsolid"), 
        # and writing the corresponding solid to a new STL file.
        text_btwn_solids = ""
        for line in lines:
            if "endsolid" in line:
                text_btwn_solids += line 
                words = line.split()
                surf = words[1]
                with open(f"{surf}.stl", 'w') as outfile:
                    outfile.write(text_btwn_solids.strip())
                    text_btwn_solids = ""
            else:
                text_btwn_solids += line + "\n"

def slice_z(mesh, n_slices=100):
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
    
    return slices 

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

    for point in slice:
        center_dist = distance2points_3d(point, center)
        min_dist = min(min_dist, center_dist)

    for p1, p2 in combinations(slice, 2):
        distance = distance2points_3d(p1, p2)
        if distance > max_dist:
            max_dist = distance
            point_a, point_b = p1, p2

    """
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
    """
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

def centerline_aprox(slices, n=100):
    """
    English:
    This function takes a list of slices as input and generates a centerline
    approximation using a spline. This function is used to create a smooth 
    and continuous representation of the centerline of the 3D mesh slices, 
    which can then be used for further analysis or visualization.

    Spanish:
    Esta función toma una lista de secciones como entrada y genera una
    aproximación de la línea central usando una spline. Se utiliza para 
    crear una representación suave y continua de la línea central de las 
    secciones de la malla 3D, que puede luego ser utilizada para el análisis 
    o la visualización.
    """
    centers = []
    for slice in slices:
        center = slice.center
        centers.append(center)
        
    # Generate spline with n interpolation points
    centerline = pv.Spline(centers, n)

    # Centerline length estimation  
    centerline_length = 0
    spline_puntos = centerline.GetPoints() 
    for i in range(n-1):
        centerline_length += distance2points_3d(spline_puntos.GetPoint(i), spline_puntos.GetPoint(i+1))
    
    # Finally, the function returns the generated spline (centerline) 
    # and the estimated centerline length.
    return centerline, centerline_length

def ortogonal_slices(centerline, mesh, n=100):
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
Folder Chooser with clearable history
    
    1. Copy and paste a filename into the combo element
    2. Use the last used item which will be visible when you create the window
    3. Choose an item from the list of previously used items
    4. Browse for a new name
    
To clear the list of previous entries, click the "Clear History" button.
The history is stored in a json file using the PySimpleGUI User Settings APIs   
"""

layout = [[sg.Combo(sorted(sg.user_settings_get_entry('-folder-', [])), default_value=sg.user_settings_get_entry('-last folder-', ''), size=(50, 1), key='-FOLDER-'), sg.FolderBrowse(), sg.B('Clear History')],
          [sg.Button('Ok', bind_return_key=True),  sg.Button('Cancel')]]

window = sg.Window('Folder Chooser With History', layout)

while True:
    event, values = window.read()

    if event in (sg.WIN_CLOSED, 'Cancel'):
        sys.exit()
    if event == 'Ok':
        # If OK, then need to add the filename to the list of files and also set as the last used filename
        sg.user_settings_set_entry('-folder-', list(set(sg.user_settings_get_entry('-folder-', []) + [values['-FOLDER-'], ])))
        sg.user_settings_set_entry('-last folder-', values['-FOLDER-'])
        entry_file = values['-FOLDER-']
        break
    elif event == 'Clear History':
        sg.user_settings_set_entry('-folder-', [])
        sg.user_settings_set_entry('-last folder-', '')
        window['-FOLDER-'].update(values=[], value='')
        sys.exit()

window.close()

# ========================================
# File Processing 
# ========================================
files = os.listdir(entry_file)
final_data = []

"""
IMPORT & SPLIT:
Use the function split_stl.py to divide our solid, it generates the following stl files:

    1. Surf_Lumen.stl
    2. Wall.stl
    3. Thrombus.stl
    4. Surf_Exterior.stl
    5. Surf_Thrombus.stl
"""

for id in files:
    split_model(entry_file + '/' + id)
    print(id)

    # Loading the four different STL files into five separate PyVista meshes.
    lumen_mesh = pv.read('Surf_Lumen.stl')
    wall_mesh = pv.read('Wall.stl')
    surf_exterior_mesh = pv.read('Surf_Exterior.stl')
    surf_thrombus_mesh = pv.read('Surf_Thrombus.stl')
    thrombus_mesh = pv.read('Thrombus.stl')

    """
    INITIAL PLOT:
    Create an interactive 3D visualization of different mesh surfaces. 
    The visualization includes a checkbox button widget to each mesh surface,
    allowing the user to toggle the visibility of each surface.
    """

    class SetVisibilityCallback:
        """Helper callback to keep a reference to the actor being modified."""
        def __init__(self, actor):
            self.actor = actor
        def __call__(self, state):
            self.actor.SetVisibility(state)

    size = 40
    Startpos = 12
    # colors = ['w', 'r', 'EC7063', 'c']
    # labels = ["Surf_Lumen", "Wall", "Surf_Exterior", "Surf_Thrombus"]
    # surfaces = [lumen_mesh, wall_mesh, surf_exterior_mesh, surf_thrombus_mesh]

    # If including the Thrombus mesh, uncomment the following lines:
    colors = ['w', 'r', 'y', 'EC7063', 'c']
    labels = ["Surf_Lumen", "Wall", "Thrombus", "Surf_Exterior", "Surf_Thrombus"]
    surfaces = [lumen_mesh, wall_mesh, thrombus_mesh, surf_exterior_mesh, surf_thrombus_mesh]

    p = pv.Plotter()
    for i, color in enumerate(colors):
        actor = p.add_mesh(surfaces[i], color=color, opacity=0.5)
        callback = SetVisibilityCallback(actor)
        p.add_checkbox_button_widget(
            callback,
            value=True,
            position=(5.0, Startpos),
            size=size,
            border_size=3,
            color_on=color,
            color_off='D5DBDB',
            background_color='D5DBDB',
        )
        p.add_text(
            labels[i], 
            position=(47, Startpos),
            font_size=15,
            font='times'
        )
        Startpos = Startpos + size + (size // 10)
    p.show()

    """
    SLICE ALONG VECTOR PARALLEL TO Z AXIS
    """
    
    lumen_slices = slice_z(lumen_mesh)
    exterior_slices = slice_z(surf_exterior_mesh)
    # Remove the first 3 slices from the lumen_slices and exterior_slices
    # becacause they tend to have a lot of noise in the surface exterior.
    for i in range(3): exterior_slices.pop(0)

    # Uncomment the following lines to visualize the slices:
    p = pv.Plotter()
    p.add_mesh(lumen_slices, opacity=0.75)
    # p.add_mesh(exterior_slices, opacity=0.75, color='r')
    # p.show()

    p = pv.Plotter(off_screen=True)
    p.add_mesh(lumen_mesh, color = '#dc143c')
    p.add_mesh(surf_exterior_mesh, style='wireframe', color = '#FFFFF0')
    path = p.generate_orbital_path(n_points=55)
    p.show(auto_close=False)
    p.open_gif('foo.gif')
    p.orbit_on_path(path, write_frames=True)
    p.close()


    p = pv.Plotter(off_screen=True)
    p.add_mesh(lumen_slices, color = 'blue')
    path = p.generate_orbital_path(n_points=55)
    p.show(auto_close=False)
    p.open_gif('foo2.gif')
    p.orbit_on_path(path, write_frames=True)
    p.close()


    """
    CENTERLINE APROXIMATION 
    """
    lumen_centerline, lumen_centerline_length = centerline_aprox(lumen_slices)
    ext_centerline, ext_centerline_length = centerline_aprox(exterior_slices, 100)


    p = pv.Plotter(off_screen=True)
    p.add_mesh(lumen_slices, color = 'blue')
    p.add_mesh(lumen_centerline, line_width=4, color='k')
    path = p.generate_orbital_path(n_points=55)
    p.show(auto_close=False)
    p.open_gif('foo3.gif')
    p.orbit_on_path(path, write_frames=True)
    p.close()


    """
    USING CENTERLINE TO GENERATE SLICES 
    """
    lumen_ortoslices = ortogonal_slices(lumen_centerline, lumen_mesh)
    exterior_ortoslices = ortogonal_slices(ext_centerline, surf_exterior_mesh, 200)

    filtered_ext_ortoslices = pv.MultiBlock()
    delete = pv.MultiBlock()
    for i, slice in enumerate(exterior_ortoslices): 
        if abs(slice.bounds[5]-slice.bounds[4]) < 20: 
            filtered_ext_ortoslices.append(slice)
        else: delete.append(slice)

    # Uncomment the following lines to visualize the ortogonal slices:
    # p = pv.Plotter()
    # p.add_mesh(lumen_ortoslices, opacity=0.75)
    # p.add_mesh(delete, opacity=0.70, color='r')
    # p.add_mesh(filtered_ext_ortoslices, opacity=0.50, color='b')
    # p.add_mesh(lumen_centerline, color="red", line_width=2)
    # p.show()

    """
    DIAMETER CALCULUS & PLOTTING 
    """
    diameters = []
    center_distance = []
    ext_diameters = []
    orto_centers = []

    for i, slice in enumerate(lumen_ortoslices):
        slice_center = slice.center
        orto_centers.append(slice_center)
        np_points = vnp.vtk_to_numpy(slice.GetPoints().GetData())
        point_a, point_b, min_dist, max_dist = diam_calc(np_points, lumen_centerline.GetPoint(i))
        diameters.append(max_dist)
        center_distance.append(min_dist)
    
    #########################################################################
    # for slice in filtered_ext_ortoslices:
    #     slice_center = slice.center
    #     np_points = vnp.vtk_to_numpy(slice.GetPoints().GetData())
    #     point_a, point_b, min_dist, max_dist = diam_calc(np_points, slice_center)
    #     ext_diameters.append(max_dist)

    """
    LOOK FOR BIFURCATION - IF BIFURCATION CLIP MESH
    """
    for position, distance in enumerate(center_distance):
        if distance < 3 and position > len(center_distance)/2:
            iliac_point = lumen_centerline.GetPoint(position)
            aux_point = lumen_centerline.GetPoint(position-5) 

            # Uncomment the following lines to visualize the iliac plane and aux plane:
            p = pv.Plotter()
            iliac_plane = pv.Plane(iliac_point, (0,0,1), i_size=50, j_size=50)
            aux_plane = pv.Plane(aux_point, (0,0,1), i_size=50, j_size=50)
            p.add_mesh(iliac_plane, color='k', opacity=0.25)
            p.add_mesh(aux_plane, color='b', opacity=0.5)
            # p.add_mesh(lumen_centerline, color='k', line_width=2)
            p.add_mesh(lumen_ortoslices)
            p.show()

            original_lumen_mesh = lumen_mesh
            original_ext_mesh = surf_exterior_mesh

            lumen_mesh = lumen_mesh.clip('z', origin=aux_point, invert=False)
            surf_thrombus_mesh = surf_thrombus_mesh.clip('z', origin=aux_point, invert=False)
            surf_exterior_mesh = surf_exterior_mesh.clip('z', origin=aux_point, invert=False)

            # Uncomment the following lines to visualize the clipped meshes:
            # p = pv.Plotter()
            # p.add_mesh(lumen_mesh, color='w')
            # p.show()

            p = pv.Plotter(off_screen=True)
            p.add_mesh(lumen_mesh, color = '#dc143c')
            p.add_mesh(surf_exterior_mesh, style='wireframe', color = '#FFFFF0')
            path = p.generate_orbital_path(n_points=55)
            p.show(auto_close=False)
            p.open_gif('foo1.gif')
            p.orbit_on_path(path, write_frames=True)
            p.close()

            # p = pv.Plotter()
            # p.add_mesh(original_lumen_mesh, color='w')
            # p.show()

            # Slice along z axis
            lumen_slices = slice_z(lumen_mesh)
            exterior_slices = slice_z(surf_exterior_mesh)
            for i in range(3): exterior_slices.pop(0)

            # Aproximation of centerline
            lumen_centerline, lumen_centerline_length = centerline_aprox(lumen_slices)
            ext_centerline, ext_centerline_length = centerline_aprox(exterior_slices)

            # Extract ortogonal slices
            lumen_ortoslices = ortogonal_slices(lumen_centerline, lumen_mesh)
            exterior_ortoslices = ortogonal_slices(ext_centerline, surf_exterior_mesh)
            filtered_ext_ortoslices = pv.MultiBlock()
            for i, slice in enumerate(exterior_ortoslices): 
                if abs(slice.bounds[5]-slice.bounds[4]) < 20: filtered_ext_ortoslices.append(slice)

            # Uncomment the following lines to visualize the ortogonal slices:
            p = pv.Plotter()
            p.add_mesh(lumen_ortoslices, opacity=0.75)
            # p.add_mesh(filtered_ext_ortoslices, opacity=0.50, color='r')
            p.add_mesh(lumen_centerline, color="r", line_width=2)
            # p.add_mesh(ext_centerline, color="k", line_width=2)
            p.show()

            diameters = []
            ext_diameters = []
            center_distance = []
            orto_centers = []
            
            for slice in lumen_ortoslices:
                slice_center = slice.center
                orto_centers.append(slice_center)
                np_points = vnp.vtk_to_numpy(slice.GetPoints().GetData())
                point_a, point_b, min_dist, max_dist = diam_calc(np_points, slice_center)
                diameters.append(max_dist)
                center_distance.append(min_dist)

            break 

    for slice in filtered_ext_ortoslices:
                slice_center = slice.center
                np_points = vnp.vtk_to_numpy(slice.GetPoints().GetData())
                point_a, point_b, min_dist, max_dist = diam_calc(np_points, slice_center)
                ext_diameters.append(max_dist)
                
    """
    PLOT RED/BLUE
    """  
    limit_diameter = 40
    p = pv.Plotter(shape=(1, 2))

    # Limit diameter:
    p.subplot(0,0)
    for i, slice in enumerate(lumen_ortoslices):  
        if diameters[i] >= limit_diameter:
            p.add_mesh(slice, opacity=0.75, color = 'red')
        else:
            p.add_mesh(slice, opacity=0.75, color = 'blue')
    p.add_mesh(lumen_centerline, color='k', line_width=2)

    # Adaptative diameter:
    p.subplot(0,1)
    dini = diameters[:5]
    initial_diameter = sum(dini)/len(dini)
    adaptative_diameter = initial_diameter*1.5

    for i, slice in enumerate(lumen_ortoslices):
        if diameters[i] >= limit_diameter or diameters[i] >= adaptative_diameter: 
            p.add_mesh(slice, opacity=0.75, color='red')
        else:
            p.add_mesh(slice, opacity=0.75, color='blue')
    p.add_mesh(lumen_centerline, color='k', line_width=2)   
    p.show()

    """
    MAX DIAMETER PLOT
    """
    # Uncomment the following lines to visualize the max diameter:
    p = pv.Plotter()
    p.add_mesh(lumen_mesh, color='w', opacity=0.5)
    p.add_mesh(surf_exterior_mesh, style='wireframe', color='k')
    p.add_mesh(lumen_ortoslices[diameters.index(max(diameters))], color='#8A2BE2', line_width=4)
    p.add_mesh(filtered_ext_ortoslices[ext_diameters.index(max(ext_diameters))], color='#ffa500', line_width=4)
    p.show()

    p = pv.Plotter(off_screen=True)
    p.add_mesh(lumen_mesh, color='w', opacity=0.5)
    p.add_mesh(surf_exterior_mesh, style='wireframe', color='k')
    p.add_mesh(lumen_ortoslices[diameters.index(max(diameters))], color='#8A2BE2', line_width=4)
    p.add_mesh(filtered_ext_ortoslices[ext_diameters.index(max(ext_diameters))], color='#ffa500', line_width=4)
    path = p.generate_orbital_path(n_points=55)
    p.show(auto_close=False)
    p.open_gif('foo4.gif')
    p.orbit_on_path(path, write_frames=True)
    p.close()

    """
    ANEURYSM VOLUME
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
        bounda = (lumen_ortoslices[min(zone_of_interest)].bounds)
        bounda2 = (lumen_ortoslices[max(zone_of_interest)].bounds)

        aneurysm_mesh = lumen_mesh.clip(normal='z', origin=[0,0,bounda2[4]], invert=False)
        aneurysm_mesh = aneurysm_mesh.clip(normal='z', origin=[0,0,bounda[5]], invert=True)

        # Other alternative:
        #aneurysm_mesh = lumen_mesh.clip(normal=normales[min(zona_interes)], origin=lumen_centerline.GetPoint(min(zona_interes)), invert=False) 
        #aneurysm_mesh = aneurysm_mesh.clip(normal=normales[max(zona_interes)], origin=lumen_centerline.GetPoint(max(zona_interes)), invert=True) 

        # Uncomment the following lines to visualize the aneurysm mesh:
        # p = pv.Plotter()
        # p.add_mesh(aneurysm_mesh, color='w')
        # p.add_mesh(lumen_mesh, style='wireframe')
        # clip_plane = pv.Plane((bounda[1], bounda[3], bounda[5]), (0,0,1), i_size=50, j_size=50) 
        # clip_plane2 = pv.Plane((bounda2[1], bounda2[3], bounda2[4]), (0,0,1), i_size=50, j_size=50)
        # p.add_mesh(clip_plane)
        # p.add_mesh(clip_plane2)
        # p.show()

        fixer = pymeshfix.MeshFix(aneurysm_mesh.triangulate())
        fixer.repair()
        repaired = fixer.mesh

        # Uncomment the following lines to visualize the closed aneurysm mesh:
        p = pv.Plotter()
        p.add_mesh(repaired, style='wireframe', color = 'k')
        # p.add_mesh(lumen_mesh, style='wireframe')
        p.show()

    """
    THROMBUS VOLUME
    """
    aux_mesh = pymeshfix.MeshFix(surf_thrombus_mesh.triangulate())
    aux_mesh.repair()
    thrombus = aux_mesh.mesh

    aux_mesh = pymeshfix.MeshFix(lumen_mesh.triangulate())
    aux_mesh.repair()
    lumen = aux_mesh.mesh

    # Alternatively:
    #aneu_thrombus_mesh = surf_thrombus_mesh.clip(normal=normales[min(zona_interes)], origin=lumen_centerline.GetPoint(min(zona_interes)), invert=False)
    #aneu_thrombus_mesh = aneu_thrombus_mesh.clip(normal=normales[max(zona_interes)], origin=lumen_centerline.GetPoint(max(zona_interes)), invert=True)
    # or:
    #aneu_thrombus_mesh = surf_thrombus_mesh.clip(normal='z', origin=[0,0,bounda2[4]], invert=False)
    #aneu_thrombus_mesh = aneu_thrombus_mesh.clip(normal='z', origin=[0,0,bounda[5]], invert=True)


    # Uncomment the following lines to visualize the thrombus mesh:
    # p = pv.Plotter()
    # p.add_mesh(lumen, color='w')
    # p.add_mesh(thrombus, color='y', style='wireframe')
    # p.show()

    """
    SIMULATION 
    """
    aux_mesh = pymeshfix.MeshFix(lumen_mesh.triangulate())
    aux_mesh.repair()
    sim_mesh = aux_mesh.mesh

    # Save the simulation meshes in other directory: 
    sim_mesh.save('sim/{}_sim.stl'.format(id), binary=False)

    # Load the simulation mesh: 
    vertices, faces = load_stl('sim/{}_sim.stl'.format(id))

    # Filter the top and bottom faces of the simulation mesh: 
    top_z = sim_mesh.bounds[5]
    bot_z = sim_mesh.bounds[4]
    
    top_faces = []
    for face in faces:
        if all(vertex[2] == top_z for vertex in face):
            top_faces.append(face)

    bot_faces = []
    for face in faces:
        if all(vertex[2] == bot_z for vertex in face):
            bot_faces.append(face)

    # Save the inlet and outlet faces in other files: 
    save_stl('inlet.stl', np.array(top_faces))
    save_stl('outlet.stl', np.array(bot_faces))

    inlet = pv.read('inlet.stl')
    outlet = pv.read('outlet.stl')

    # Uncomment the following lines to visualize the inlet and outlet faces:
    p = pv.Plotter()
    p.add_mesh(sim_mesh, style='wireframe', color = '#dc143c')
    # p.add_mesh(inlet, color = 'green')
    # p.add_mesh(outlet, color = 'r')
    p.show()

    total_volume = sim_mesh.volume

    """
    GEOMETRICAL AAA PARAMETERS 
    """

    ## HYPOTHETIC AORTA LENGTH:
    hypothesis_length = distance2points_3d(lumen_centerline.GetPoint(0), lumen_centerline.GetPoint(lumen_centerline.n_points-1))

    ## ANEURYSM LENGTH - LAAA :
    LAAA = 0
    p = pv.Plotter()
    if aneurysm == True:
        for i in range(max(zone_of_interest) - min(zone_of_interest)):
                LAAA += distance2points_3d(lumen_centerline.GetPoint(min(zone_of_interest) + i), lumen_centerline.GetPoint(min(zone_of_interest) + i + 1))
                point_array = np.array(lumen_centerline.GetPoint(min(zone_of_interest) + i + 1))
                # p.add_points(point_array, render_points_as_spheres=True)
        p.add_mesh(aneurysm_mesh, style='wireframe', color='r')
        p.show()
        # print('ANEURYSM LENGTH (LAAA) = ', LAAA)
    else:
        LAAA = 0
        # print('ANEURYSM LENGTH (LAAA) = ', LAAA)

    ## MAXIMUM DIAMETER OF THE ANEURYSM - DAMAX:
    DAMAX = max(diameters)
    # print('MAXIMUM DIAMETER OF THE ANEURYSM (DAMAX) = ', DAMAX)

    ## MAXIMUM EXTERIOR DIAMETER - D_EXT_MAX:
    D_EXT_MAX = max(ext_diameters)
    # print('MAXIMUM EXTERIOR DIAMETER (D_EXT_MAX) = ', D_EXT_MAX)

    ## SACCULAR INDEX:
    if aneurysm == True:
        sac_index = DAMAX/LAAA
    else:
        sac_index = 0
    # print('SACCULAR INDEX = ', sac_index)

    ## CENTERLINE LENGTH:
    # print('CENTERLINE LENGTH = ', lumen_centerline_length)

    ## TORTUOSITY:
    tortuosity = lumen_centerline_length/lumen_centerline.length 
    # print('TORTUOSITY = ', tortuosity)

    ## VOLUME AAA:
    if aneurysm == True:
        volume_AAA = repaired.volume
    else:
        volume_AAA = 0
    # print("ANEURYSM VOLUME = {:.3f}".format(volume_AAA))

    ## THROMBUS VOLUME:
    if aneurysm == True:
        th_volume = thrombus.volume - lumen.volume
    else:
        th_volume = 0
    # print("ANEURYSM THROMBUS VOLUME = {:.3f}".format(th_volume))


    test_data = {'id' : id, 
                'HYPOTHETIC AORTA LENGTH' : hypothesis_length, 
                'ANEURYSM LENGTH' : LAAA,
                'DAMAX' : DAMAX,
                'D_EXT_MAX' : D_EXT_MAX,
                'SACCULAR INDEX' : sac_index, 
                'CENTERLINE LENGTH' : lumen_centerline_length,
                'TORTUOSITY' : tortuosity,
                'ANEURYSM VOLUME' : volume_AAA,
                'THROMBUS VOLUME' : th_volume,
                'TOTAL VOLUME' : total_volume
                }
    
    final_data.append(test_data)

    # Delete the temporary files: 
    os.remove('Surf_Exterior.stl')
    os.remove('Surf_Lumen.stl')
    os.remove('Surf_Thrombus.stl')
    os.remove('Wall.stl')
    # os.remove('Thrombus.stl')


# Save the data in a json file: 
with open('datos.json', "w") as file:
    json.dump(final_data, file, indent=4)