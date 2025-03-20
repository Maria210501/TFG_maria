import numpy as np
from stl import mesh

def read_stl(file_path):
    """Lee un archivo STL y devuelve nodos y elementos."""
    stl_mesh = mesh.Mesh.from_file(file_path)
    
    # Extraer los nodos únicos
    nodes, indices = np.unique(stl_mesh.vectors.reshape(-1, 3), axis=0, return_inverse=True)
    
    # Crear los elementos triangulares
    elements = indices.reshape(-1, 3)
    
    return nodes, elements

def write_mdpa(nodes, elements, output_file):
    """Escribe el archivo MDPA con los nodos y elementos."""
    with open(output_file, 'w') as f:
        f.write('Begin ModelPartData\nEnd ModelPartData\n\n')
        f.write('Begin Properties 1\nEnd Properties\n\n')
        
        # Escribir nodos
        f.write('Begin Nodes\n')
        for i, node in enumerate(nodes, start=1):
            f.write(f' {i} {node[0]} {node[1]} {node[2]}\n')
        f.write('End Nodes\n\n')
        
        # Escribir elementos
        f.write('Begin Elements Element3D3N\n')
        for i, element in enumerate(elements, start=1):
            f.write(f' {i} 1 {element[0]+1} {element[1]+1} {element[2]+1}\n')
        f.write('End Elements\n')

# Uso del código
stl_file = "modelo.stl"  # Cambia esto con el nombre de tu STL
mdpa_file = "modelo.mdpa"
nodes, elements = read_stl(stl_file)
write_mdpa(nodes, elements, mdpa_file)
print(f"Archivo MDPA generado: {mdpa_file}")