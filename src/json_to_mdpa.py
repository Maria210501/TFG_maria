import json

def json_to_mdpa(json_file, mdpa_file):
    """
    Converts a JSON file to an MDPA file for KRATOS.

    Parameters:
    json_file (str): Path to the input JSON file.
    mdpa_file (str): Path to the output MDPA file.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    with open(mdpa_file, 'w') as f:
        f.write("Begin ModelPartData\n")
        f.write("//  Variable Name value\n")
        f.write("End ModelPartData\n\n")

        f.write("Begin Properties 1\n")
        f.write("// Definir, por ejemplo, las propiedades del material (módulo de Young, coeficiente de Poisson, densidad, etc.) o parámetros específicos que afecten el comportamiento del modelo durante la simulación.\n")
        f.write("End Properties\n\n")

        f.write("Begin Nodes\n")
        for item in data:
            f.write(f"{item['id']} {item['HYPOTHETIC AORTA LENGTH']} {item['ANEURYSM LENGTH']} {item['DAMAX']} {item['D_EXT_MAX']} {item['SACCULAR INDEX']} {item['CENTERLINE LENGTH']} {item['TORTUOSITY']} {item['ANEURYSM VOLUME']} {item['THROMBUS VOLUME']} {item['TOTAL VOLUME']}\n")
        f.write("End Nodes\n\n")

        f.write("Begin Elements Element2D3N\n")
        for item in data:
            f.write(f"1 1 {item['id']} {item['HYPOTHETIC AORTA LENGTH']} {item['ANEURYSM LENGTH']}\n")
        f.write("End Elements\n\n")

        f.write("Begin NodalData DISPLACEMENT_X\n")
        for item in data:
            f.write(f"{item['id']} 0 0.0\n")
        f.write("End NodalData\n\n")

        f.write("Begin NodalData DISPLACEMENT_Y\n")
        for item in data:
            f.write(f"{item['id']} 0 0.0\n")
        f.write("End NodalData\n\n")

        f.write("Begin NodalData DISPLACEMENT_Z\n")
        for item in data:
            f.write(f"{item['id']} 0 0.0\n")
        f.write("End NodalData\n")

# Example usage
json_to_mdpa('/Users/mariarodriguez/Desktop/TFG/TFG_maria/src/data/datos_control.json', '/Users/mariarodriguez/Desktop/TFG/TFG_maria/src/data/datos_control.mdpa')

json_to_mdpa('/Users/mariarodriguez/Desktop/TFG/TFG_maria/src/data/datos_rupture.json', '/Users/mariarodriguez/Desktop/TFG/TFG_maria/src/data/datos_rupture.mdpa')