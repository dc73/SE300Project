# =============================================================================
# Title: Mesh to FEA Conversion Script
# Author: Dheer Chhabria
# Course: SE 300
# Date: 01-03-2025
#
# Description:
#
# Basically this script converts a 3D mesh file (OBJ or STL) into a FEA
# input file. It extracts 2D node data and element connectivity from the mesh, applies
# default material properties, boundary conditions, and loads, and then wrtes the data
# in a specific format required for FEA.
# =============================================================================

import tkinter as tk
from tkinter import simpledialog, filedialog
import trimesh
import numpy as np

def convert_mesh_to_fea(input_filename, output_filename): 


    mesh = trimesh.load(input_filename)
    

    vertices = mesh.vertices
    nodes2d = vertices[:, :2]


    node_list = []
    for i, coord in enumerate(nodes2d):
        node_list.append((i+1, coord[0], coord[1]))


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

    # Material and analysis parameters
    Analysis_Type = 1         # 1 for plane stress, 2 for plane strain
    E = 41000000.0            # Young's Modulus
    Poisson = 0.33            # Poisson's Ratio
    density = 0.00000783      # Material Density
    Thickness = 1.0           # Thickness of the element (use 1 for plane strain)
    AccelX = 0.0              # Acceleration in the X direction
    AccelY = -9.81            # Acceleration in the Y direction

    NNodes = len(node_list)
    NElements = len(element_list)
    
    
    max_node_id = max([node[0] for node in node_list])
    bc_list = [
        (nid, 1, 1) for nid in range(1, min(14, max_node_id + 1))
    ]
    NBC = len(bc_list)
    
    
    load_list = [
        (nid, 2000, 0) for nid in range(1, min(6, max_node_id + 1))
    ]
    
    
    root = tk.Tk()
    root.withdraw()
    

    force_value = simpledialog.askfloat("Input", "Enter the force value:")
    force_node = simpledialog.askinteger("Input", "Enter the node number to apply the force:")

    if force_value is not None and force_node is not None and 1 <= force_node <= max_node_id:
        load_list.append((force_node, force_value, 0))
    
    NLOAD = len(load_list)
    

    with open(output_filename, 'w') as f:

        f.write("#---------------------MATERIAL DATA\n")
        f.write(f"{Analysis_Type}                  # 1 for plane stress or 2 for plane strain\n")
        f.write(f"{E}           # Young's Modulus\n")
        f.write(f"{Poisson}               # Poisson Ratio\n")
        f.write(f"{density}         # Material Density\n")
        f.write(f"{Thickness}                  # Thickness (Leave as 1 if it is a plane strain Analysis)\n")
        f.write(f"{AccelX}                  # Acceleration X\n")
        f.write(f"{AccelY}              # Acceleration Y\n")
        f.write(f"{NNodes}                  # Number of nodes\n")
        f.write(f"{NElements}                  # Number of elements\n")
        f.write(f"{NBC}                  # Number of BC's applied\n")
        f.write(f"{NLOAD}                  # Number of Loads Applied\n")
        

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
    
    
    print(f"Conversion complete! Output written to: {output_filename}")

def main():
    """
    Main function to prompt the user for file inputs and initiate the mesh conversion.
    """
    
    root = tk.Tk()
    root.withdraw()
    
    input_filename = filedialog.askopenfilename(
        title="Select Mesh File",
        filetypes=[("OBJ files", "*.obj"), ("STL files", "*.stl"), ("All files", "*.*")]
    )
    if not input_filename:
        print("No input file selected. Exiting.")
        return

    # Open a dialog to select the output file location and name
    output_filename = filedialog.asksaveasfilename(
        title="Save FEA Input File",
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not output_filename:
        print("No output file selected. Exiting.")
        return

    # Convert the mesh to an FEA input file
    convert_mesh_to_fea(input_filename, output_filename)

if __name__ == "__main__":
    main()
