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
        for j in range(1, 4):  # Omitir el primer valor (tamaño del 'face')
            node = mesh.points[cell[j]]
            nodes[node_id] = node
            element_nodes.append(node_id)
            node_id += 1
        elements.append(element_nodes)

    # Guardar el archivo MDPA
    with open(mdpa_file, 'w') as f:
        f.write("Begin Nodes\n")
        for node_id, node in nodes.items():
            f.write(f"{node_id} {node[0]} {node[1]} {node[2]}\n")
        f.write("End Nodes\n")
        f.write("Begin Elements\n")
        for element in elements:
            f.write(f"1 {' '.join(map(str, element))}\n")
        f.write("End Elements\n")

def generate_mdpa_from_stl_folder(stl_folder, output_folder):
    """Genera archivos MDPA para todos los archivos STL en una carpeta y los guarda en la carpeta de salida."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for stl_file in os.listdir(stl_folder):
        if stl_file.endswith(".stl"):
            stl_path = os.path.join(stl_folder, stl_file)
            mdpa_filename = os.path.splitext(stl_file)[0] + ".mdpa"
            mdpa_path = os.path.join(output_folder, mdpa_filename)
            stl_to_mdpa_pyvista(stl_path, mdpa_path)

# Ejemplo de uso
stl_folder = '/Users/mariarodriguez/solids/STL_controls'  # Reemplaza con tu carpeta de STL
output_folder = '/Users/mariarodriguez/solids/mdpa_controls'  # Reemplaza con la carpeta donde deseas guardar los MDPA
generate_mdpa_from_stl_folder(stl_folder, output_folder)@
