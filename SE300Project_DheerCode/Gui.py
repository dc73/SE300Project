import tkinter as tk
from tkinter import ttk

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SE300 FEA Software")
        self.geometry("1200x700")
        self.configure(bg="#ECF0F1")  # Light background for the main window

        # Set up ttk styling for a modern look
        style = ttk.Style(self)
        style.theme_use("clam")
        
        # Configure styles for the right panel and its widgets
        style.configure("Right.TFrame", background="#ECF0F1")
        style.configure("Modern.TButton",
                        font=("Helvetica", 12),
                        padding=10,
                        foreground="white",
                        background="#3498DB",
                        borderwidth=0)
        style.map("Modern.TButton", background=[('active', '#2980B9')])
        style.configure("TLabel", background="#ECF0F1", font=("Helvetica", 12))

        # ---- Left Panel ----
        # Use tk.Frame for more flexible background styling
        self.left_frame = tk.Frame(self, bg="#2C3E50", width=250)
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y)
        
        # ---- Right Panel ----
        self.right_frame = ttk.Frame(self, style="Right.TFrame", padding=(20,20))
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.create_left_panel()
        self.create_right_panel()

    def create_left_panel(self):
        # Custom style for left panel buttons (using tk.Button for full control)
        btn_font = ("Helvetica", 12)
        btn_style = {
            "font": btn_font,
            "bg": "#34495E",
            "fg": "white",
            "activebackground": "#3E5871",
            "activeforeground": "white",
            "relief": tk.FLAT,
            "bd": 0,
            "highlightthickness": 0,
            "padx": 10,
            "pady": 10
        }
        
        # List of button texts for the left panel
        buttons = ["Import CAD File", "Import Mesh File", "Assign Material", "Compute Geometry"]
        for text in buttons:
            btn = tk.Button(self.left_frame, text=text, **btn_style)
            btn.pack(fill=tk.X, pady=10, padx=20)

    def create_right_panel(self):
        # Top row frame for buttons in the right panel
        top_frame = ttk.Frame(self.right_frame)
        top_frame.pack(fill=tk.X, pady=(0, 15))
        
        # Top row buttons with the modern style
        btn_set_boundary = ttk.Button(top_frame, text="Set Boundary Conditions", style="Modern.TButton")
        btn_set_boundary.pack(side=tk.LEFT, padx=5)
        
        btn_save_calc = ttk.Button(top_frame, text="Save Calculation", style="Modern.TButton")
        btn_save_calc.pack(side=tk.LEFT, padx=5)
        
        btn_post_process = ttk.Button(top_frame, text="Choose Post Process", style="Modern.TButton")
        btn_post_process.pack(side=tk.LEFT, padx=5)
        
        # "Results" label with larger, bold font
        results_label = ttk.Label(self.right_frame, text="Results", font=("Helvetica", 16, "bold"))
        results_label.pack(anchor=tk.NW, pady=(0, 5))
        
        # Placeholder frame for 3D plot or results visualization
        plot_frame = ttk.Frame(self.right_frame, relief=tk.SUNKEN)
        plot_frame.pack(fill=tk.BOTH, expand=True, pady=(0,15))
        
        # Error message label at bottom
        error_label = ttk.Label(self.right_frame,
                                text="Solving Error - Please Check Input Conditions",
                                foreground="red", font=("Helvetica", 12, "bold"))
        error_label.pack(side=tk.BOTTOM, pady=10)

if __name__ == "__main__":
    app = App()
    app.mainloop()
