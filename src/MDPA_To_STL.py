import numpy as np
from stl import mesh

def stl_to_mdpa(stl_file, mdpa_file):
    # Cargar la malla STL
    stl_mesh = mesh.Mesh.from_file(stl_file)
    
    # Extraer nodos y triangular la superficie
    nodes = {}
    elements = []
    node_id = 1

    for i, triangle in enumerate(stl_mesh.vectors):
        element_nodes = []
        for vertex in triangle:
            vertex_tuple = tuple(np.round(vertex, decimals=6))  # Redondeo para evitar duplicados
            if vertex_tuple not in nodes:
                nodes[vertex_tuple] = node_id
                node_id += 1
            element_nodes.append(nodes[vertex_tuple])
        elements.append(element_nodes)
    
    # Escribir archivo MDPA
    with open(mdpa_file, "w") as f:
        f.write("Begin ModelPartData\n//  Propiedades globales de la malla\nEnd ModelPartData\n\n")
        f.write("Begin Properties 1\n// Definir materiales y constantes aquí\nEnd Properties\n\n")
        
        # Escribir nodos
        f.write("Begin Nodes\n")
        for (x, y, z), id in sorted(nodes.items(), key=lambda item: item[1]):
            f.write(f" {id} {x} {y} {z}\n")
        f.write("End Nodes\n\n")

        # Escribir elementos (triángulos de superficie)
        f.write("Begin Elements SurfaceCondition3D3N\n")
        for i, element in enumerate(elements, start=1):
            f.write(f" {i} 1 {element[0]} {element[1]} {element[2]}\n")
        f.write("End Elements\n")

    print(f"Archivo MDPA generado: {mdpa_file}")

# Uso del script
stl_to_mdpa("/Users/mariarodriguez/k3014.stl", "modelo.mdpa")
