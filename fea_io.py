import json

class FEAIO:
    """Handles file input/output operations."""
    
    def write_conversion_file(self, node_list, element_list, bc_list, load_list,
                              material_props, analysis_type, filename="INPUT_FEA_PROTUS_3.txt"):
        try:
            with open(filename, "w") as f:
                f.write("#---------------------MATERIAL DATA\n")
                f.write(f"{analysis_type}\n")
                f.write(f"{material_props['E']}\n")
                f.write(f"{material_props['Poisson']}\n")
                f.write(f"{material_props['Density']}\n")
                f.write("1.0\n")  # Thickness
                f.write("0.0\n")  # Acceleration X
                f.write("-9.81\n")  # Acceleration Y
                f.write(f"{len(node_list)}\n")
                f.write(f"{len(element_list)}\n")
                f.write(f"{len(bc_list)}\n")
                f.write(f"{len(load_list)}\n")
                f.write("#---------------------NODE DATA\n")
                for node in node_list:
                    f.write(f"{node[0]},{node[1]},{node[2]}\n")
                f.write("#---------------------ELEMENT DATA\n")
                for elem in element_list:
                    f.write(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]}\n")
                f.write("#---------------------BC, 1-Fixed or 0-Free\n")
                for bc in bc_list:
                    f.write(f"{bc[0]},{bc[1]},{bc[2]}\n")
                f.write("#---------------------LOAD AND BC\n")
                for load in load_list:
                    f.write(f"{load[0]},{load[1]},{load[2]}\n")
            return True
        except Exception as e:
            print(f"Error writing conversion file: {e}")
            return False
    
    def export_file(self, data, save_path):
        try:
            with open(save_path, "w") as f:
                json.dump(data, f, indent=4)
            return True
        except Exception as e:
            print(f"Error exporting file: {e}")
            return False
    
    def load_file(self, load_path):
        try:
            with open(load_path, "r") as f:
                data = json.load(f)
            return data
        except Exception as e:
            print(f"Error loading file: {e}")
            return None
