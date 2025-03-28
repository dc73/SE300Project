import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import trimesh
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
import io
import sys

from fea_computation import FEAComputation
from fea_data import FEAData
from fea_io import FEAIO

class FEAGUI(tk.Tk):
    """Main GUI application class for FEA Software."""
    def __init__(self):
        super().__init__()
        self.configure(bg="white")
        self.title("FEA Software")
        self.geometry("1800x800")
        
        # Initialize supporting classes
        self.io = FEAIO()
        self.computation = FEAComputation()
        self.data_storage = FEAData()
        
        self.analysis_results = None  # (node_coords, element_data, results)
        self.create_widgets()
    
    def create_widgets(self):
        # Left panel for mesh conversion and material selection
        left_frame = tk.Frame(self, bg="white", width=500)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        # Right panel for analysis and visualization
        right_frame = tk.Frame(self, bg="white")
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # ------------------ Left Section ------------------
        tk.Label(left_frame, text="Mesh Conversion & Material Selection", bg="white",
                 fg="black", font=("Helvetica", 14, "bold")).pack(pady=10)
        tk.Label(left_frame, text="Mesh File (STL):", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W)
        self.mesh_path_var = tk.StringVar()
        mesh_frame = tk.Frame(left_frame, bg="white")
        mesh_frame.pack(fill=tk.X, pady=5)
        tk.Entry(mesh_frame, textvariable=self.mesh_path_var, width=40, bg="white",
                 fg="black", bd=0, highlightthickness=1, relief=tk.FLAT).pack(side=tk.LEFT, padx=5)
        tk.Button(mesh_frame, text="Browse", command=self.browse_mesh_file, bg="white", fg="black").pack(side=tk.LEFT)
        tk.Label(left_frame, text="Select Material:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        self.material_var = tk.StringVar()
        materials = ["Aluminum 6061", "Titanium Ti-6Al-4V", "Carbon Fiber Composite",
                     "Stainless Steel 304", "Inconel 718", "Mild Steel", "Copper", "Brass"]
        self.material_dropdown = ttk.Combobox(left_frame, textvariable=self.material_var,
                                              values=materials, state="readonly", width=40)
        self.material_dropdown.pack(pady=5)
        self.material_dropdown.set("Select Material")
        self.material_dropdown.bind("<<ComboboxSelected>>", self.update_material_properties)
        prop_frame = tk.LabelFrame(left_frame, text="Material Properties", bg="white", fg="black", font=("Helvetica", 12))
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
        tk.Label(top_right, text="FEA Analysis & Visualization", bg="white", fg="black", font=("Helvetica", 14, "bold")).pack(side=tk.LEFT, padx=10, pady=10)
        tk.Button(top_right, text="Run Analysis", command=self.run_analysis, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Export Results", command=self.export_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Load Old Results", command=self.load_old_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Label(top_right, text="Visualization Mode:", bg="white", fg="black", font=("Helvetica", 12)).pack(side=tk.LEFT, padx=10)
        self.vis_mode_var = tk.StringVar(value="2D Heatmap")
        vis_modes = ["2D Heatmap", "3D Heatmap"]
        self.vis_dropdown = ttk.Combobox(top_right, textvariable=self.vis_mode_var,
                                         values=vis_modes, state="readonly", width=15)
        self.vis_dropdown.pack(side=tk.LEFT, padx=5)
        tk.Label(top_right, text="Result Type:", bg="white", fg="black", font=("Helvetica", 12)).pack(side=tk.LEFT, padx=10)
        self.result_type_var = tk.StringVar(value="Stress X")
        result_types = ["Stress X", "Stress Y", "Strain X", "Strain Y", "Von Mises Stress", "Integration"]
        self.result_dropdown = ttk.Combobox(top_right, textvariable=self.result_type_var,
                                            values=result_types, state="readonly", width=15)
        self.result_dropdown.pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="View Results", command=self.view_results, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        self.log_text = tk.Text(right_frame, height=15, bg="white", fg="black")
        self.log_text.pack(fill=tk.X, padx=10, pady=5)
        self.vis_frame = tk.Frame(right_frame, bg="white", relief=tk.SUNKEN, bd=2)
        self.vis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
    
    def browse_mesh_file(self):
        filename = filedialog.askopenfilename(
            title="Select Mesh File",
            filetypes=[("STL files", "*.stl"), ("All files", "*.*")]
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
            mesh = trimesh.load(self.mesh_path_var.get())
            vertices = mesh.vertices
            nodes2d = vertices[:, :2]
            node_list = []
            for i, coord in enumerate(nodes2d):
                node_list.append((i + 1, coord[0], coord[1]))
            faces = mesh.faces
            element_list = []
            for i, face in enumerate(faces):
                face = list(face)
                if len(face) == 3:
                    element_list.append((i + 1, face[0] + 1, face[1] + 1, face[2] + 1, face[2] + 1))
                elif len(face) == 4:
                    element_list.append((i + 1, face[0] + 1, face[1] + 1, face[2] + 1, face[3] + 1))
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
            material_props = materials[self.material_var.get()]
            Analysis_Type = 1 if self.analysis_type_var.get() == "Plane Stress" else 2
            max_node_id = max(node[0] for node in node_list)
            bc_list = [(nid, 1, 1) for nid in range(1, min(14, max_node_id + 1))]
            load_list = []
            try:
                force_node = int(self.force_node_var.get())
                force_value = float(self.force_value_var.get())
                if 1 <= force_node <= max_node_id:
                    load_list.append((force_node, force_value, 1))
            except:
                pass
            if self.io.write_conversion_file(node_list, element_list, bc_list, load_list,
                                             material_props, Analysis_Type):
                messagebox.showinfo("Success", "Mesh conversion complete! 'INPUT_FEA_PROTUS_3.txt' has been written.")
            else:
                messagebox.showerror("Error", "Failed to write conversion file.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
    
    def export_conversion(self):
        try:
            with open("INPUT_FEA_PROTUS_3.txt", "r") as f:
                data = f.read()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read conversion file: {str(e)}")
            return
        save_path = filedialog.asksaveasfilename(title="Export Conversion File",
                                                 defaultextension=".txt",
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
        node_coords, element_data, results = FEAComputation.run_analysis()
        output = sys.stdout.getvalue()
        sys.stdout = old_stdout
        self.log_text.insert(tk.END, output)
        if node_coords is not None:
            self.analysis_results = (node_coords, element_data, results)
            self.data_storage.node_coords = node_coords
            self.data_storage.element_data = element_data
            self.data_storage.results = results
        else:
            self.analysis_results = None
    
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
        save_path = filedialog.asksaveasfilename(title="Export Results File",
                                                 defaultextension=".json",
                                                 filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if save_path:
            if self.io.export_file(export_data, save_path):
                messagebox.showinfo("Success", f"Results exported to:\n{save_path}")
            else:
                messagebox.showerror("Error", "Error saving results.")
    
    def load_old_results(self):
        load_path = filedialog.askopenfilename(title="Load Old Results File",
                                               filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if load_path:
            data = self.io.load_file(load_path)
            if data:
                node_coords = [tuple(coord) for coord in data.get("node_coords", [])]
                element_data = data.get("element_data", [])
                results = data.get("results", {})
                self.analysis_results = (node_coords, element_data, results)
                self.data_storage.node_coords = node_coords
                self.data_storage.element_data = element_data
                self.data_storage.results = results
                messagebox.showinfo("Success", "Old results loaded successfully!")
            else:
                messagebox.showerror("Error", "Error loading results.")
    
    def view_results(self):
        for widget in self.vis_frame.winfo_children():
            widget.destroy()
        if self.analysis_results is None:
            messagebox.showerror("Error", "No analysis results available. Please run analysis or load old results.")
            return
        node_coords, element_data, results = self.analysis_results
        result_type = self.result_type_var.get()
        if result_type not in results:
            messagebox.showerror("Error", f"Result type '{result_type}' not found.")
            return
        if self.vis_mode_var.get() == "2D Heatmap":
            self.plot_2d_result(self.vis_frame, node_coords, element_data, results[result_type], result_type)
        else:
            self.plot_3d_result(self.vis_frame, node_coords, element_data, results[result_type], result_type)
    
    def plot_2d_result(self, parent, node_coords, element_data, result_values, result_label):
        fig, ax = plt.subplots(figsize=(6, 6))
        fig.patch.set_facecolor('#4A4A4A')
        ax.set_facecolor('#4A4A4A')
        for elem in element_data:
            indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
            poly_x = [node_coords[i][0] for i in indices]
            poly_y = [node_coords[i][1] for i in indices]
            poly_x.append(poly_x[0])
            poly_y.append(poly_y[0])
            ax.plot(poly_x, poly_y, color='black', linewidth=1)
        x_vals = [pt[0] for pt in node_coords]
        y_vals = [pt[1] for pt in node_coords]
        sc = ax.scatter(x_vals, y_vals, c=result_values, cmap='viridis', edgecolors='black')
        ax.set_title(f"{result_label} (2D)", color='white')
        ax.set_xlabel("X Axis", color='white')
        ax.set_ylabel("Y Axis", color='white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        cbar = fig.colorbar(sc, ax=ax)
        cbar.ax.yaxis.set_tick_params(color='white')
        cbar.outline.set_edgecolor('white')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
        ax.set_aspect('equal', 'box')
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.draw()
        widget = canvas.get_tk_widget()
        widget.pack(fill=tk.BOTH, expand=True)
        return canvas
    
    def plot_3d_result(self, parent, node_coords, element_data, result_values, result_label):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        fig.patch.set_facecolor('#4A4A4A')
        ax.set_facecolor('#4A4A4A')
        for elem in element_data:
            indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
            poly_x = [node_coords[i][0] for i in indices]
            poly_y = [node_coords[i][1] for i in indices]
            poly_z = [0] * len(poly_x)
            poly_x.append(poly_x[0])
            poly_y.append(poly_y[0])
            poly_z.append(poly_z[0])
            ax.plot(poly_x, poly_y, poly_z, color='black', linewidth=1)
        x_vals = [pt[0] for pt in node_coords]
        y_vals = [pt[1] for pt in node_coords]
        z_vals = result_values
        sc = ax.scatter(x_vals, y_vals, z_vals, c=result_values, cmap='viridis', edgecolors='black')
        ax.set_title(f"{result_label} (3D)", color='white')
        ax.set_xlabel("X Axis", color='white')
        ax.set_ylabel("Y Axis", color='white')
        ax.set_zlabel(result_label, color='white')
        ax.xaxis.set_tick_params(colors='white')
        ax.yaxis.set_tick_params(colors='white')
        ax.zaxis.set_tick_params(colors='white')
        ax.w_xaxis.line.set_color('white')
        ax.w_yaxis.line.set_color('white')
        ax.w_zaxis.line.set_color('white')
        cbar = fig.colorbar(sc, ax=ax, pad=0.1)
        cbar.ax.yaxis.set_tick_params(color='white')
        cbar.outline.set_edgecolor('white')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
        ax.view_init(elev=25, azim=-60)
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.draw()
        widget = canvas.get_tk_widget()
        widget.pack(fill=tk.BOTH, expand=True)
        return canvas

if __name__ == "__main__":
    app = FEAGUI()
    app.mainloop()
