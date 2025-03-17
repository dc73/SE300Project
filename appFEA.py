import tkinter as tk
from tkinter import simpledialog, filedialog

def main():
    """
    Main function to prompt the user for file inputs and initiate the mesh conversion.
    """
    
    start = tk.Tk()
    start.withdraw()
    
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
    
if __name__ == "__main__":
    main()