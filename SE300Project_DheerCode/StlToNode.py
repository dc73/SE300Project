import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import trimesh
import numpy as np

class MeshToFEAConverter:
    def __init__(self, master):
        self.master = master
        master.title("Mesh to FEA Converter")
        master.geometry("600x700")

        # Aerospace Materials Dictionary
        self.materials = {
            "Aluminum 6061": {
                "E": 68900000000.0,  # Pa
                "Poisson": 0.33,
                "Density": 2700.0,   # kg/m³
            },
            "Titanium Ti-6Al-4V": {
                "E": 113800000000.0,  # Pa
                "Poisson": 0.342,
                "Density": 4430.0,   # kg/m³
            },
            "Carbon Fiber Composite": {
                "E": 181000000000.0,  # Pa
                "Poisson": 0.30,
                "Density": 1600.0,   # kg/m³
            },
            "Stainless Steel 304": {
                "E": 193000000000.0,  # Pa
                "Poisson": 0.29,
                "Density": 7900.0,   # kg/m³
            },
            "Inconel 718": {
                "E": 179000000000.0,  # Pa
                "Poisson": 0.31,
                "Density": 8190.0,   # kg/m³
            }
        }

        # Create Widgets
        self.create_widgets()

    def create_widgets(self):
        # Mesh File Selection
        tk.Label(self.master, text="Mesh File (STL/OBJ):").pack(pady=(10,0))
        self.mesh_frame = tk.Frame(self.master)
        self.mesh_frame.pack(pady=(0,10))
        
        self.mesh_path = tk.StringVar()
        self.mesh_entry = tk.Entry(self.mesh_frame, textvariable=self.mesh_path, width=50)
        self.mesh_entry.pack(side=tk.LEFT, padx=(0,10))
        
        tk.Button(self.mesh_frame, text="Browse", command=self.browse_mesh_file).pack(side=tk.LEFT)

        # Material Selection
        tk.Label(self.master, text="Select Material:").pack()
        self.material_var = tk.StringVar()
        self.material_dropdown = ttk.Combobox(
            self.master, 
            textvariable=self.material_var, 
            values=list(self.materials.keys()),
            state="readonly",
            width=50
        )
        self.material_dropdown.pack(pady=(0,10))
        self.material_dropdown.set("Select Material")
        self.material_dropdown.bind("<<ComboboxSelected>>", self.update_material_properties)

        # Material Properties Frame
        self.prop_frame = tk.LabelFrame(self.master, text="Material Properties")
        self.prop_frame.pack(padx=20, pady=10, fill="x")

        # Young's Modulus
        tk.Label(self.prop_frame, text="Young's Modulus (Pa):").pack()
        self.e_var = tk.StringVar()
        self.e_entry = tk.Entry(self.prop_frame, textvariable=self.e_var, state="readonly", width=50)
        self.e_entry.pack()

        # Poisson's Ratio
        tk.Label(self.prop_frame, text="Poisson's Ratio:").pack()
        self.poisson_var = tk.StringVar()
        self.poisson_entry = tk.Entry(self.prop_frame, textvariable=self.poisson_var, state="readonly", width=50)
        self.poisson_entry.pack()

        # Density
        tk.Label(self.prop_frame, text="Density (kg/m³):").pack()
        self.density_var = tk.StringVar()
        self.density_entry = tk.Entry(self.prop_frame, textvariable=self.density_var, state="readonly", width=50)
        self.density_entry.pack()

        # Analysis Type
        tk.Label(self.master, text="Analysis Type:").pack()
        self.analysis_var = tk.StringVar(value="Plane Stress")
        self.analysis_dropdown = ttk.Combobox(
            self.master, 
            textvariable=self.analysis_var, 
            values=["Plane Stress", "Plane Strain"],
            state="readonly",
            width=50
        )
        self.analysis_dropdown.pack(pady=(0,10))

        # Force Input
        tk.Label(self.master, text="Force Application:").pack()
        self.force_frame = tk.Frame(self.master)
        self.force_frame.pack(pady=(0,10))

        tk.Label(self.force_frame, text="Node:").pack(side=tk.LEFT)
        self.force_node_var = tk.StringVar()
        self.force_node_entry = tk.Entry(self.force_frame, textvariable=self.force_node_var, width=10)
        self.force_node_entry.pack(side=tk.LEFT, padx=(0,10))

        tk.Label(self.force_frame, text="Force (N):").pack(side=tk.LEFT)
        self.force_value_var = tk.StringVar()
        self.force_value_entry = tk.Entry(self.force_frame, textvariable=self.force_value_var, width=10)
        self.force_value_entry.pack(side=tk.LEFT)

        # Convert Button
        tk.Button(self.master, text="Convert Mesh to FEA", command=self.convert_mesh).pack(pady=20)

    def browse_mesh_file(self):
        """Open file dialog to select mesh file"""
        filename = filedialog.askopenfilename(
            title="Select Mesh File",
            filetypes=[("STL files", "*.stl"), ("OBJ files", "*.obj"), ("All files", "*.*")]
        )
        if filename:
            self.mesh_path.set(filename)

    def update_material_properties(self, event=None):
        """Update material properties based on selected material"""
        material = self.material_var.get()
        if material in self.materials:
            props = self.materials[material]
            self.e_var.set(f"{props['E']:.2e}")
            self.poisson_var.set(str(props['Poisson']))
            self.density_var.set(str(props['Density']))

    def convert_mesh(self):
        """Convert mesh to FEA input file"""
        # Validate inputs
        if not self.mesh_path.get():
            messagebox.showerror("Error", "Please select a mesh file!")
            return

        if self.material_var.get() == "Select Material":
            messagebox.showerror("Error", "Please select a material!")
            return

        try:
            # Load mesh
            mesh = trimesh.load(self.mesh_path.get())

            # Select output file
            output_filename = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
            )
            if not output_filename:
                return

            # Extract vertices and faces
            vertices = mesh.vertices
            nodes2d = vertices[:, :2]

            # Create node list
            node_list = []
            for i, coord in enumerate(nodes2d):
                node_list.append((i+1, coord[0], coord[1]))

            # Create element list
            faces = mesh.faces
            element_list = []
            for i, face in enumerate(faces):
                face = list(face)
                if len(face) == 3:
                    element_list.append((i+1, face[0]+1, face[1]+1, face[2]+1, face[2]+1))
                elif len(face) == 4:
                    element_list.append((i+1, face[0]+1, face[1]+1, face[2]+1, face[3]+1))
                else:
                    print(f"Warning: Face {i+1} has {len(face)} vertices; skipping it.")

            # Material properties from selected material
            material_props = self.materials[self.material_var.get()]
            
            # Analysis type
            Analysis_Type = 1 if self.analysis_var.get() == "Plane Stress" else 2
            
            # Create boundary conditions and loads
            max_node_id = max([node[0] for node in node_list])
            bc_list = [(nid, 1, 1) for nid in range(1, min(14, max_node_id + 1))]
            
            # Process force input
            load_list = []
            try:
                force_node = int(self.force_node_var.get())
                force_value = float(self.force_value_var.get())
                
                if 1 <= force_node <= max_node_id:
                    # 0 for x-direction, 1 for y-direction
                    load_list.append((force_node, force_value, 1))
            except (ValueError, TypeError):
                pass

            # Write FEA input file
            with open(output_filename, 'w') as f:
                # Material and analysis data
                f.write("#---------------------MATERIAL DATA\n")
                f.write(f"{Analysis_Type}                  # 1 for plane stress or 2 for plane strain\n")
                f.write(f"{material_props['E']}           # Young's Modulus\n")
                f.write(f"{material_props['Poisson']}               # Poisson Ratio\n")
                f.write(f"{material_props['Density']}         # Material Density\n")
                f.write("1.0                  # Thickness (Leave as 1 if plane strain)\n")
                f.write("0.0                  # Acceleration X\n")
                f.write("-9.81              # Acceleration Y\n")
                f.write(f"{len(node_list)}                  # Number of nodes\n")
                f.write(f"{len(element_list)}                  # Number of elements\n")
                f.write(f"{len(bc_list)}                  # Number of BC's applied\n")
                f.write(f"{len(load_list)}                  # Number of Loads Applied\n")

                # Node data
                f.write("#---------------------NODE DATA\n")
                for node in node_list:
                    f.write(f"{node[0]},{node[1]},{node[2]}\n")

                # Element data
                f.write("#---------------------ELEMENT DATA\n")
                for elem in element_list:
                    f.write(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]}\n")

                # Boundary conditions
                f.write("#---------------------BC, 1-Fixed or 0-Free\n")
                for bc in bc_list:
                    f.write(f"{bc[0]},{bc[1]},{bc[2]}\n")

                # Loads
                f.write("#---------------------LOAD AND BC\n")
                for load in load_list:
                    f.write(f"{load[0]},{load[1]},{load[2]}\n")

            messagebox.showinfo("Success", f"Conversion complete! Output written to: {output_filename}")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

def main():
    root = tk.Tk()
    app = MeshToFEAConverter(root)
    root.mainloop()

if __name__ == "__main__":
    main()