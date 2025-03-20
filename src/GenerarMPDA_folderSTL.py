import os
import shutil
import pyvista as pv
import numpy as np
import tetgen

def stl_to_mdpa_pyvista(stl_file, mdpa_file):
    """ Convierte un archivo STL en MDPA usando PyVista. """
    # Cargar el STL
    mesh = pv.read(stl_file)

    # Obtener nodos y elementos
    nodes = {}
    elements = []
    node_id = 1

    for cell in mesh.faces.reshape(-1, 4):  # Cada celda tiene (3 vértices + 1 tamaño)
        element_nodes = []
        for j in range(1, 4):  # Omitir el primer valor (tamaño del triángulo)
            vertex = tuple(np.round(mesh.points[cell[j]], decimals=6))  # Redondeo para evitar duplicados
            if vertex not in nodes:
                nodes[vertex] = node_id
                node_id += 1
            element_nodes.append(nodes[vertex])
        elements.append(element_nodes)

    # Escribir archivo MDPA
    with open(mdpa_file, "w") as f:
        f.write("Begin ModelPartData\n// Propiedades globales\nEnd ModelPartData\n\n")
        f.write("Begin Properties 1\n// Definir materiales aquí\nEnd Properties\n\n")

        # Escribir nodos
        f.write("Begin Nodes\n")
        for (x, y, z), id in sorted(nodes.items(), key=lambda item: item[1]):
            f.write(f" {id} {x} {y} {z}\n")
        f.write("End Nodes\n\n")

        # Escribir elementos
        f.write("Begin Elements SurfaceCondition3D3N\n")
        for i, element in enumerate(elements, start=1):
            f.write(f" {i} 1 {element[0]} {element[1]} {element[2]}\n")
        f.write("End Elements\n")

    print(f"MDPA generado: {mdpa_file}")

def generar_malla_volumetrica(stl_file, vtk_file):
    """ Genera una malla de volumen tetraédrica a partir de un STL y la guarda en VTK. """
    
    print(f"Procesando {stl_file}...")
    
    # Cargar la superficie STL
    surface_mesh = pv.read(stl_file)

    # Convertir la malla a formato adecuado para TetGen
    tet = tetgen.TetGen(surface_mesh)

    # Parámetros de calidad del mallado
    tet_mesh = tet.tetrahedralize(
        quality=True,       # Activa el control de calidad
        minratio=1.5,       # Relación mínima de calidad
        maxvolume=0.1,     # Volumen máximo de cada tetraedro
        mindihedral=10      # Ángulo mínimo permitido en los tetraedros
    )

    # Guardar la malla en formato VTK
    tet_mesh.save(vtk_file)
    print(f"Malla volumétrica guardada en: {vtk_file}")

def process_all_stl(directory):
    """ Procesa todos los STL en la carpeta especificada. """
    # Obtener todos los archivos STL en el directorio
    for filename in os.listdir(directory):
        if filename.lower().endswith(".stl"):
            stl_path = os.path.join(directory, filename)
            base_name = os.path.splitext(filename)[0]  # Obtener nombre sin extensión
            
            # Crear carpeta específica si no existe
            folder_path = os.path.join(directory, base_name)
            os.makedirs(folder_path, exist_ok=True)

            # Copiar el STL a su carpeta (en lugar de moverlo)
            new_stl_path = os.path.join(folder_path, filename)
            shutil.copy(stl_path, new_stl_path)

            # Ruta del archivo MDPA
            mdpa_path = os.path.join(folder_path, f"{base_name}.mdpa")
            
            # Generar MDPA
            stl_to_mdpa_pyvista(new_stl_path, mdpa_path)
            print(f"Procesado: {filename}")

            # Ruta para la malla volumétrica en VTK
            vtk_path = os.path.join(folder_path, f"{base_name}.vtk")

            # Generar malla volumétrica  ==> Demomento esto no se ejecuta
            #generar_malla_volumetrica(new_stl_path, vtk_path)
            print(f"Procesado: {filename}")

# Directorio actual (puedes cambiarlo a otra ruta si es necesario)
directorio_trabajo = "/Users/mariarodriguez/solids/STL_rupture"
process_all_stl(directorio_trabajo)
