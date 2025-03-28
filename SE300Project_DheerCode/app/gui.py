import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sys
import io
import json

from fea_calculations import RunAnalysis, plot_2d_result
from stl_to_node import convert_stl_to_node

class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.configure(bg="white")
        self.title("FEA Software")
        self.geometry("1800x800")
        self.analysis_results = None  # Will hold (node_coords, element_data, results)
        self.create_widgets()

    def create_widgets(self):
        # Left panel (Mesh conversion & material selection)
        left_frame = tk.Frame(self, bg="white", width=500)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        # Right panel (Analysis & visualization)
        right_frame = tk.Frame(self, bg="white")
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ------------------ Left Section ------------------
        tk.Label(left_frame, text="Mesh Conversion & Material Selection", bg="white", fg="black",
                 font=("Helvetica", 14, "bold")).pack(pady=10)
        tk.Label(left_frame, text="Mesh File (STL/OBJ):", bg="white", fg="black",
                 font=("Helvetica", 12)).pack(anchor=tk.W)
        self.mesh_path_var = tk.StringVar()
        mesh_frame = tk.Frame(left_frame, bg="white")
        mesh_frame.pack(fill=tk.X, pady=5)
        tk.Entry(mesh_frame, textvariable=self.mesh_path_var, width=40, bg="white", fg="black",
                 bd=0, highlightthickness=1, relief=tk.FLAT).pack(side=tk.LEFT, padx=5)
        tk.Button(mesh_frame, text="Browse", command=self.browse_mesh_file, bg="white", fg="black").pack(side=tk.LEFT)
        
        tk.Label(left_frame, text="Mesh Density (1 = Coarse, 10 = Fine):", bg="white", fg="black",
                 font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        self.mesh_density_var = tk.IntVar(value=5)
        tk.Scale(left_frame, from_=1, to=10, orient=tk.HORIZONTAL, variable=self.mesh_density_var,
                 bg="white").pack(fill=tk.X, pady=5)
        
        tk.Label(left_frame, text="Select Material:", bg="white", fg="black",
                 font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        self.material_var = tk.StringVar()
        materials = ["Aluminum 6061", "Titanium Ti-6Al-4V", "Carbon Fiber Composite",
                     "Stainless Steel 304", "Inconel 718", "Mild Steel", "Copper", "Brass"]
        self.material_dropdown = ttk.Combobox(left_frame, textvariable=self.material_var, values=materials,
                                              state="readonly", width=40)
        self.material_dropdown.pack(pady=5)
        self.material_dropdown.set("Select Material")
        self.material_dropdown.bind("<<ComboboxSelected>>", self.update_material_properties)
        
        prop_frame = tk.LabelFrame(left_frame, text="Material Properties", bg="white", fg="black",
                                   font=("Helvetica", 12))
        prop_frame.pack(fill=tk.X, pady=5, padx=5)
        tk.Label(prop_frame, text="Young's Modulus (Pa):", bg="white", fg="black", font=("Helvetica", 10)).pack(anchor=tk.W)
        self.e_var = tk.StringVar()
        tk.Entry(prop_frame, textvariable=self.e_var, state="readonly", width=40, bg="black", fg="white").pack(pady=2)
        tk.Label(prop_frame, text="Poisson's Ratio:", bg="white", fg="black", font=("Helvetica", 10)).pack(anchor=tk.W)
        self.poisson_var = tk.StringVar()
        tk.Entry(prop_frame, textvariable=self.poisson_var, state="readonly", width=40, bg="black", fg="white").pack(pady=2)
        tk.Label(prop_frame, text="Density (kg/m³):", bg="white", fg="black", font=("Helvetica", 10)).pack(anchor=tk.W)
        self.density_var = tk.StringVar()
        tk.Entry(prop_frame, textvariable=self.density_var, state="readonly", width=40, bg="black", fg="white").pack(pady=2)
        
        tk.Label(left_frame, text="Analysis Type:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        self.analysis_type_var = tk.StringVar(value="Plane Stress")
        self.analysis_dropdown = ttk.Combobox(left_frame, textvariable=self.analysis_type_var,
                                              values=["Plane Stress", "Plane Strain"], state="readonly", width=40)
        self.analysis_dropdown.pack(pady=5)
        
        tk.Label(left_frame, text="Force Application:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        force_frame = tk.Frame(left_frame, bg="white")
        force_frame.pack(pady=5)
        tk.Label(force_frame, text="Node:", bg="white", fg="black", font=("Helvetica", 10)).pack(side=tk.LEFT)
        self.force_node_var = tk.StringVar()
        tk.Entry(force_frame, textvariable=self.force_node_var, width=5, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Label(force_frame, text="Force (N):", bg="white", fg="black", font=("Helvetica", 10)).pack(side=tk.LEFT)
        self.force_value_var = tk.StringVar()
        tk.Entry(force_frame, textvariable=self.force_value_var, width=5, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        
        tk.Button(left_frame, text="Convert Mesh to FEA", command=self.convert_mesh, bg="white", fg="black").pack(pady=15)
        tk.Button(left_frame, text="Export Conversion File", command=self.export_conversion, bg="white", fg="black").pack(pady=5)
        
        # ------------------ Right Section ------------------
        top_right = tk.Frame(right_frame, bg="white")
        top_right.pack(fill=tk.X)
        tk.Label(top_right, text="FEA Analysis & Visualization", bg="white", fg="black",
                 font=("Helvetica", 14, "bold")).pack(side=tk.LEFT, padx=10, pady=10)
        tk.Button(top_right, text="Run Analysis", command=self.run_analysis, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Export Results", command=self.export_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Load Old Results", command=self.load_old_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Label(top_right, text="Visualization Mode:", bg="white", fg="black", font=("Helvetica", 12)).pack(side=tk.LEFT, padx=10)
        self.vis_mode_var = tk.StringVar(value="2D Heatmap")
        vis_modes = ["2D Heatmap"]
        self.vis_dropdown = ttk.Combobox(top_right, textvariable=self.vis_mode_var, values=vis_modes, state="readonly", width=15)
        self.vis_dropdown.pack(side=tk.LEFT, padx=5)
        tk.Label(top_right, text="Result Type:", bg="white", fg="black", font=("Helvetica", 12)).pack(side=tk.LEFT, padx=10)
        self.result_type_var = tk.StringVar(value="Stress X")
        result_types = ["Stress X", "Stress Y", "Strain X", "Strain Y", "Von Mises Stress", "Integration"]
        self.result_dropdown = ttk.Combobox(top_right, textvariable=self.result_type_var, values=result_types, state="readonly", width=15)
        self.result_dropdown.pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="View Results", command=self.view_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Export VTK for Paraview", command=self.export_vtk, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        
        self.log_text = tk.Text(right_frame, height=15, bg="white", fg="black")
        self.log_text.pack(fill=tk.X, padx=10, pady=5)
        self.vis_frame = tk.Frame(right_frame, bg="white", relief=tk.SUNKEN, bd=2)
        self.vis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def browse_mesh_file(self):
        filename = filedialog.askopenfilename(
            title="Select Mesh File",
            filetypes=[("STL files", "*.stl"), ("OBJ files", "*.obj"), ("All files", "*.*")]
        )
        if filename:
            self.mesh_path_var.set(filename)

    def update_material_properties(self, event=None):
        materials = {
            "Aluminum 6061": {"E": 68900000000.0, "Poisson": 0.33, "Density": 2700.0},
            "Titanium Ti-6Al-4V": {"E": 113800000000.0, "Poisson": 0.342, "Density": 4430.0},
            "Carbon Fiber Composite": {"E": 181000000000.0, "Poisson": 0.30, "Density": 1600.0},
            "Stainless Steel 304": {"E": 193000000000.0, "Poisson": 0.29, "Density": 7900.0},
            "Inconel 718": {"E": 179000000000.0, "Poisson": 0.31, "Density": 8190.0},
            "Mild Steel": {"E": 210e9, "Poisson": 0.3, "Density": 7850},
            "Copper": {"E": 110e9, "Poisson": 0.34, "Density": 8960},
            "Brass": {"E": 100e9, "Poisson": 0.34, "Density": 8500}
        }
        material = self.material_var.get()
        if material in materials:
            props = materials[material]
            self.e_var.set(f"{props['E']:.2e}")
            self.poisson_var.set(str(props['Poisson']))
            self.density_var.set(str(props['Density']))

    def convert_mesh(self):
        if not self.mesh_path_var.get():
            messagebox.showerror("Error", "Please select a mesh file!")
            return
        if self.material_var.get() == "Select Material":
            messagebox.showerror("Error", "Please select a material!")
            return
        try:
            material_keys = {
                "Aluminum 6061": {"E": 68900000000.0, "Poisson": 0.33, "Density": 2700.0},
                "Titanium Ti-6Al-4V": {"E": 113800000000.0, "Poisson": 0.342, "Density": 4430.0},
                "Carbon Fiber Composite": {"E": 181000000000.0, "Poisson": 0.30, "Density": 1600.0},
                "Stainless Steel 304": {"E": 193000000000.0, "Poisson": 0.29, "Density": 7900.0},
                "Inconel 718": {"E": 179000000000.0, "Poisson": 0.31, "Density": 8190.0},
                "Mild Steel": {"E": 210e9, "Poisson": 0.3, "Density": 7850},
                "Copper": {"E": 110e9, "Poisson": 0.34, "Density": 8960},
                "Brass": {"E": 100e9, "Poisson": 0.34, "Density": 8500}
            }
            material_props = material_keys[self.material_var.get()]
            Analysis_Type = 1 if self.analysis_type_var.get() == "Plane Stress" else 2
            # Use the convert_stl_to_node function to create the input file.
            input_filename = convert_stl_to_node(
                self.mesh_path_var.get(),
                self.mesh_density_var.get(),
                material_props,
                Analysis_Type,
                self.force_node_var.get(),
                self.force_value_var.get()
            )
            messagebox.showinfo("Success", f"Mesh conversion complete!\n'{input_filename}' has been written.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

    def export_conversion(self):
        try:
            with open("INPUT_FEA_3.txt", "r") as f:
                data = f.read()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read conversion file: {str(e)}")
            return
        save_path = filedialog.asksaveasfilename(title="Export Conversion File", defaultextension=".txt",
                                                 filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if save_path:
            try:
                with open(save_path, "w") as f:
                    f.write(data)
                messagebox.showinfo("Success", f"Conversion file exported to:\n{save_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Error saving file: {str(e)}")

    def run_analysis(self):
        self.log_text.delete("1.0", tk.END)
        for widget in self.vis_frame.winfo_children():
            widget.destroy()
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        node_coords, element_data, results = RunAnalysis()
        output = sys.stdout.getvalue()
        sys.stdout = old_stdout
        self.log_text.insert(tk.END, output)
        if node_coords is not None:
            self.analysis_results = (node_coords, element_data, results)
        else:
            self.analysis_results = None
        messagebox.showinfo("FEA", "Analysis complete!")
    
    def export_vtk(self):
        try:
            with open("SE300Paraview.vtk", "r") as f:
                vtk_data = f.read()
            save_path = filedialog.asksaveasfilename(title="Export VTK File", defaultextension=".vtk",
                                                     filetypes=[("VTK files", "*.vtk"), ("All files", "*.*")])
            if save_path:
                with open(save_path, "w") as f:
                    f.write(vtk_data)
                messagebox.showinfo("Success", f"VTK file exported to:\n{save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Error exporting VTK file: {str(e)}")
    
    def export_results(self):
        if self.analysis_results is None:
            messagebox.showerror("Error", "No analysis results available. Run analysis first.")
            return
        node_coords, element_data, results = self.analysis_results
        export_data = {
            "node_coords": [list(coord) for coord in node_coords],
            "element_data": element_data,
            "results": results
        }
        save_path = filedialog.asksaveasfilename(title="Export Results File", defaultextension=".json",
                                                 filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if save_path:
            try:
                with open(save_path, "w") as f:
                    json.dump(export_data, f, indent=4)
                messagebox.showinfo("Success", f"Results exported to:\n{save_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Error saving results: {str(e)}")
    
    def load_old_results(self):
        load_path = filedialog.askopenfilename(title="Load Old Results File",
                                               filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if load_path:
            try:
                with open(load_path, "r") as f:
                    data = json.load(f)
                node_coords = [tuple(coord) for coord in data.get("node_coords", [])]
                element_data = data.get("element_data", [])
                results = data.get("results", {})
                self.analysis_results = (node_coords, element_data, results)
                messagebox.showinfo("Success", "Old results loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading results: {str(e)}")
    
    def view_results(self):
        for widget in self.vis_frame.winfo_children():
            widget.destroy()
        if self.analysis_results is None:
            messagebox.showerror("Error", "No analysis results available. Please run analysis or load old results.")
            return
        node_coords, element_data, results = self.analysis_results
        result_key = self.result_type_var.get()
        if result_key not in results:
            messagebox.showerror("Error", f"Result type '{result_key}' not found.")
            return
        # For demonstration, we use a 2D heatmap to display the chosen result.
        plot_2d_result(self.vis_frame, node_coords, element_data, results[result_key], result_key)

# Additional class "App" as a simple entry point for running the program.
class App:
    def __init__(self):
        self.gui = MainApp()
    def run(self):
        self.gui.mainloop()

if __name__ == "__main__":
    app = App()
    app.run()
