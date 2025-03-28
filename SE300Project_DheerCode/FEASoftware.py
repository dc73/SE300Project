import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import re
import math
import io
import sys
import trimesh
from copy import deepcopy
import json
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D


# ------------------------- MATRIX OPERATIONS -------------------------
def MatrixMultiply(MatrixA, MatrixB):
    RowsA = len(MatrixA)
    ColsA = len(MatrixA[0])
    ColsB = len(MatrixB[0])
    Result = [[0 for _ in range(ColsB)] for _ in range(RowsA)]
    for i in range(RowsA):
        for j in range(ColsB):
            for k in range(ColsA):
                Result[i][j] += MatrixA[i][k] * MatrixB[k][j]
    return Result


def MatrixAdd(MatrixA, MatrixB):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    Result = [[0 for _ in range(Cols)] for _ in range(Rows)]
    for i in range(Rows):
        for j in range(Cols):
            Result[i][j] = MatrixA[i][j] + MatrixB[i][j]
    return Result


def MatrixScalarMultiply(MatrixA, Scalar):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    return [[MatrixA[i][j] * Scalar for j in range(Cols)] for i in range(Rows)]


def MatrixTranspose(MatrixA):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    return [[MatrixA[i][j] for i in range(Rows)] for j in range(Cols)]


# ------------------------- FEA ANALYSIS FUNCTION -------------------------
def RunAnalysis():
    """
    Reads the fixed FEA input file "INPUT_FEA_PROTUS_3.txt" (generated during mesh conversion),
    performs (dummy) FEA calculations, and returns:
      - node_coords: list of (x, y) for each node,
      - element_data: list of element connectivity (each element: [id, n1, n2, n3, n4]),
      - results: a dictionary mapping result type to a list of nodal values.
    """
    print("*******************************************************************************")
    print("Processing FEA Analysis...")
    try:
        with open("INPUT_FEA_PROTUS_3.txt", "r") as f:
            InputLines = f.readlines()
    except FileNotFoundError:
        print("Error: 'INPUT_FEA_PROTUS_3.txt' not found. Please convert a mesh first.")
        return None, None, None

    # Read header parameters (using raw strings to avoid escape warnings)
    AnalysisType = int(re.split(r'\s+', InputLines[1])[0])
    ElasticModulus = float(re.split(r'\s+', InputLines[2])[0])
    PoissonsRatio = float(re.split(r'\s+', InputLines[3])[0])
    MaterialDensity = float(re.split(r'\s+', InputLines[4])[0])
    Thickness = float(re.split(r'\s+', InputLines[5])[0])
    AccelX = float(re.split(r'\s+', InputLines[6])[0])
    AccelY = float(re.split(r'\s+', InputLines[7])[0])
    NumNodes = int(re.split(r'\s+', InputLines[8])[0])
    NumElements = int(re.split(r'\s+', InputLines[9])[0])
    NumBC = int(re.split(r'\s+', InputLines[10])[0])
    NumLoads = int(re.split(r'\s+', InputLines[11])[0])

    DOFTotal = 2 * NumNodes
    DOFPerNode = 2

    print("END OF USER INPUT PROCESSING")
    print("================= SOLVER ======================")
    print("NUMBER OF NODES IS -", NumNodes)
    print("NUMBER OF ELEMENTS IS -", NumElements)
    print("TOTAL NUMBER OF VARIABLES IN THE MODEL -", DOFTotal)

    # ------------------------- READ NODE DATA -------------------------
    NodeData = [[0, 0, 0] for _ in range(NumNodes)]
    Index = 13
    while Index < (NumNodes + 13):
        NodeValues = [float(val.strip()) for val in re.split(r',', InputLines[Index])]
        NodeData[Index - 13][0] = NodeValues[0]  # node id
        NodeData[Index - 13][1] = NodeValues[1]  # x coordinate
        NodeData[Index - 13][2] = NodeValues[2]  # y coordinate
        Index += 1

    # ------------------------- READ ELEMENT DATA -------------------------
    ElementData = []
    Index = 14 + NumNodes
    while Index < (14 + NumNodes + NumElements):
        ElementValues = [int(val.strip()) for val in re.split(r',', InputLines[Index])]
        ElementData.append(ElementValues[:5])
        Index += 1

    # (BC and Loads reading omitted here for brevity)

    # ------------------------- Dummy Post-Processing -------------------------
    # Create a dummy results matrix. In your real FEA code, you would fill these values from post-processing.
    NumNodes = len(NodeData)
    ResultsMatrix = [[None] * 72 for _ in range(NumNodes)]
    for i in range(NumNodes):
        # Dummy assignments: these could be replaced with real computed values.
        ResultsMatrix[i][1] = NodeData[i][1] * 0.2  # Stress X
        ResultsMatrix[i][6] = NodeData[i][2] * 0.2  # Stress Y
        ResultsMatrix[i][66] = NodeData[i][1] * 0.3  # Strain X
        ResultsMatrix[i][70] = NodeData[i][2] * 0.3  # Strain Y
        ResultsMatrix[i][30] = NodeData[i][1] * 0.4  # Von Mises Stress
        ResultsMatrix[i][15] = NodeData[i][1] * 0.5  # Integration dummy

    results = {
        "Stress X": [ResultsMatrix[i][1] for i in range(NumNodes)],
        "Stress Y": [ResultsMatrix[i][6] for i in range(NumNodes)],
        "Strain X": [ResultsMatrix[i][66] for i in range(NumNodes)],
        "Strain Y": [ResultsMatrix[i][70] for i in range(NumNodes)],
        "Von Mises Stress": [ResultsMatrix[i][30] for i in range(NumNodes)],
        "Integration": [ResultsMatrix[i][15] for i in range(NumNodes)]
    }

    print("THE ANALYSIS HAS COMPLETED SUCCESSFULLY")
    node_coords = [(NodeData[i][1], NodeData[i][2]) for i in range(NumNodes)]
    return node_coords, ElementData, results


# ------------------------- UPDATED PLOTTING FUNCTIONS -------------------------
def plot_2d_result(parent, node_coords, element_data, result_values, result_label):
    """
    2D wireframe with a colored heatmap for the nodal results.
    """
    # Create figure/axis
    fig, ax = plt.subplots(figsize=(6, 6))

    # --- Set background to gray ---
    fig.patch.set_facecolor('#4A4A4A')  # A medium-dark gray for the figure background
    ax.set_facecolor('#4A4A4A')  # Same gray for the axes background

    # --- Plot wireframe for each element ---
    for elem in element_data:
        # elem = [elem_id, n1, n2, n3, n4], each node index is 1-based
        indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
        poly_x = [node_coords[i][0] for i in indices]
        poly_y = [node_coords[i][1] for i in indices]
        # Close the polygon
        poly_x.append(poly_x[0])
        poly_y.append(poly_y[0])
        ax.plot(poly_x, poly_y, color='black', linewidth=1)

    # --- Scatter the nodal results ---
    x_vals = [pt[0] for pt in node_coords]
    y_vals = [pt[1] for pt in node_coords]
    sc = ax.scatter(x_vals, y_vals, c=result_values, cmap='viridis', edgecolors='black')

    # --- Format axis, labels, and colorbar ---
    ax.set_title(f"{result_label} (2D)", color='white')
    ax.set_xlabel("X Axis", color='white')
    ax.set_ylabel("Y Axis", color='white')
    # Change tick colors to white
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

    # Add colorbar on the right
    cbar = fig.colorbar(sc, ax=ax)
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.outline.set_edgecolor('white')
    # Change colorbar tick labels to white
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    # Optional: fix aspect ratio so the mesh doesn't look skewed
    ax.set_aspect('equal', 'box')

    # Embed in Tkinter
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.draw()
    widget = canvas.get_tk_widget()
    widget.pack(fill=tk.BOTH, expand=True)
    return canvas


def plot_3d_result(parent, node_coords, element_data, result_values, result_label):
    """
    3D wireframe with colored scatter (Z = result_values).
    """
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    # --- Set background to gray ---
    fig.patch.set_facecolor('#4A4A4A')
    ax.set_facecolor('#4A4A4A')  # 3D axes also accept facecolor

    # --- Plot wireframe in XY plane (z=0) ---
    for elem in element_data:
        indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
        poly_x = [node_coords[i][0] for i in indices]
        poly_y = [node_coords[i][1] for i in indices]
        poly_z = [0] * len(poly_x)
        # Close the polygon
        poly_x.append(poly_x[0])
        poly_y.append(poly_y[0])
        poly_z.append(poly_z[0])
        ax.plot(poly_x, poly_y, poly_z, color='black', linewidth=1)

    # --- Plot nodes with Z = result_values, color mapped ---
    x_vals = [pt[0] for pt in node_coords]
    y_vals = [pt[1] for pt in node_coords]
    z_vals = result_values
    sc = ax.scatter(x_vals, y_vals, z_vals, c=result_values, cmap='viridis', edgecolors='black')

    # --- Format axis, labels, and colorbar ---
    ax.set_title(f"{result_label} (3D)", color='white')
    ax.set_xlabel("X Axis", color='white')
    ax.set_ylabel("Y Axis", color='white')
    ax.set_zlabel(result_label, color='white')

    # Change tick colors to white
    ax.xaxis.set_tick_params(colors='white')
    ax.yaxis.set_tick_params(colors='white')
    ax.zaxis.set_tick_params(colors='white')

    # 3D axis lines can also be recolored, if needed:
    ax.w_xaxis.line.set_color('white')
    ax.w_yaxis.line.set_color('white')
    ax.w_zaxis.line.set_color('white')

    cbar = fig.colorbar(sc, ax=ax, pad=0.1)
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.outline.set_edgecolor('white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    # Optionally set a decent view angle
    ax.view_init(elev=25, azim=-60)

    # Embed in Tkinter
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.draw()
    widget = canvas.get_tk_widget()
    widget.pack(fill=tk.BOTH, expand=True)
    return canvas


# ------------------------- MAIN GUI APPLICATION -------------------------
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
        tk.Label(left_frame, text="Mesh Conversion & Material Selection", bg="white", fg="black",font=("Helvetica", 14, "bold")).pack(pady=10)
        tk.Label(left_frame, text="Mesh File (STL/OBJ):", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W)
        self.mesh_path_var = tk.StringVar()
        mesh_frame = tk.Frame(left_frame, bg="white")
        mesh_frame.pack(fill=tk.X, pady=5)
        tk.Entry(mesh_frame, textvariable=self.mesh_path_var, width=40, bg="white", fg="black",bd=0,highlightthickness=1,relief=tk.FLAT).pack(side=tk.LEFT,padx=5)
        tk.Button(mesh_frame, text="Browse", command=self.browse_mesh_file, bg="white", fg="black").pack(side=tk.LEFT)
        tk.Label(left_frame, text="Select Material:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W,pady=(10, 0))
        self.material_var = tk.StringVar()
        materials = ["Aluminum 6061", "Titanium Ti-6Al-4V", "Carbon Fiber Composite", "Stainless Steel 304", "Inconel 718", "Mild Steel", "Copper", "Brass"]
        self.material_dropdown = ttk.Combobox(left_frame, textvariable=self.material_var,values=materials, state="readonly", width=40)
        self.material_dropdown.pack(pady=5)
        self.material_dropdown.set("Select Material")
        self.material_dropdown.bind("<<ComboboxSelected>>", self.update_material_properties)
        prop_frame = tk.LabelFrame(left_frame, text="Material Properties", bg="white", fg="black",font=("Helvetica", 12))
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
        tk.Label(left_frame, text="Analysis Type:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W,pady=(10, 0))
        self.analysis_type_var = tk.StringVar(value="Plane Stress")
        self.analysis_dropdown = ttk.Combobox(left_frame, textvariable=self.analysis_type_var,values=["Plane Stress", "Plane Strain"], state="readonly", width=40)
        self.analysis_dropdown.pack(pady=5)
        tk.Label(left_frame, text="Force Application:", bg="white", fg="black", font=("Helvetica", 12)).pack(anchor=tk.W, pady=(10, 0))
        force_frame = tk.Frame(left_frame, bg="white")
        force_frame.pack(pady=5)
        tk.Label(force_frame, text="Node:", bg="white", fg="black", font=("Helvetica", 10)).pack(side=tk.LEFT)
        self.force_node_var = tk.StringVar()
        tk.Entry(force_frame, textvariable=self.force_node_var, width=5, bg="white", fg="black").pack(side=tk.LEFT,padx=5)
        tk.Label(force_frame, text="Force (N):", bg="white", fg="black", font=("Helvetica", 10)).pack(side=tk.LEFT)
        self.force_value_var = tk.StringVar()
        tk.Entry(force_frame, textvariable=self.force_value_var, width=5, bg="white", fg="black").pack(side=tk.LEFT,padx=5)
        tk.Button(left_frame, text="Convert Mesh to FEA", command=self.convert_mesh, bg="white", fg="black").pack(pady=15)
        tk.Button(left_frame, text="Export Conversion File", command=self.export_conversion, bg="white",fg="black").pack(pady=5)

        # ------------------ Right Section ------------------
        top_right = tk.Frame(right_frame, bg="white")
        top_right.pack(fill=tk.X)
        tk.Label(top_right, text="FEA Analysis & Visualization", bg="white", fg="black", font=("Helvetica", 14, "bold")).pack(side=tk.LEFT, padx=10, pady=10)
        tk.Button(top_right, text="Run Analysis", command=self.run_analysis, bg="white", fg="black").pack(side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Export Results", command=self.export_results, bg="white", fg="black").pack(
            side=tk.LEFT, padx=5)
        tk.Button(top_right, text="Load Old Results", command=self.load_old_results, bg="white", fg="black").pack(
            side=tk.LEFT, padx=5)
        tk.Label(top_right, text="Visualization Mode:", bg="white", fg="black", font=("Helvetica", 12)).pack(
            side=tk.LEFT, padx=10)
        self.vis_mode_var = tk.StringVar(value="2D Heatmap")
        vis_modes = ["2D Heatmap", "3D Heatmap"]
        self.vis_dropdown = ttk.Combobox(top_right, textvariable=self.vis_mode_var, values=vis_modes, state="readonly", width=15)
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
            # Write the FEA input file using the fixed filename "INPUT_FEA_PROTUS_3.txt"
            with open("INPUT_FEA_PROTUS_3.txt", "w") as f:
                f.write("#---------------------MATERIAL DATA\n")
                f.write(f"{Analysis_Type}\n")
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
            messagebox.showinfo("Success", "Mesh conversion complete! 'INPUT_FEA_PROTUS_3.txt' has been written.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

    def export_conversion(self):
        # Allows the user to choose a location to save a copy of the conversion file.
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
        node_coords, element_data, results = RunAnalysis()
        output = sys.stdout.getvalue()
        sys.stdout = old_stdout
        self.log_text.insert(tk.END, output)
        if node_coords is not None:
            self.analysis_results = (node_coords, element_data, results)
        else:
            self.analysis_results = None

    def export_results(self):
        if self.analysis_results is None:
            messagebox.showerror("Error", "No analysis results available. Run analysis first.")
            return
        # Prepare the results data as a dictionary.
        node_coords, element_data, results = self.analysis_results
        # Convert any tuples to lists for JSON serialization.
        export_data = {
            "node_coords": [list(coord) for coord in node_coords],
            "element_data": element_data,
            "results": results
        }
        save_path = filedialog.asksaveasfilename(title="Export Results File",
                                                 defaultextension=".json",
                                                 filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if save_path:
            try:
                import json
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
                import json
                with open(load_path, "r") as f:
                    data = json.load(f)
                # Restore data (convert lists back to tuples if desired)
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
        result_type = self.result_type_var.get()
        if result_type not in results:
            messagebox.showerror("Error", f"Result type '{result_type}' not found.")
            return
        if self.vis_mode_var.get() == "2D Heatmap":
            plot_2d_result(self.vis_frame, node_coords, element_data, results[result_type], result_type)
        else:
            plot_3d_result(self.vis_frame, node_coords, element_data, results[result_type], result_type)


if __name__ == "__main__":
    app = MainApp()
    app.mainloop()
