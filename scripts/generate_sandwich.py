import os
import sys

# 添加项目根目录到Python路径
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
from carbon_easy_builder.graphene import Graphene
from carbon_easy_builder.nanotube import CarbonNanotube
from carbon_easy_builder.lammps_writer import LAMMPSWriter
from carbon_easy_builder.box import Box

class SandwichStructureGenerator:
    def __init__(self):
        self.output_dir = "scripts/output"
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Create GUI
        self.root = tk.Tk()
        self.root.title("CNT-Graphene Sandwich Structure Generator")
        self.root.geometry("600x800")
        
        # Create input fields
        self.create_input_fields()
        
        # Create buttons
        self.create_buttons()
        
        # Create output text area
        self.create_output_area()
        
    def create_input_fields(self):
        # CNT Parameters
        cnt_frame = ttk.LabelFrame(self.root, text="CNT Parameters", padding=10)
        cnt_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(cnt_frame, text="CNT Radius (Å):").grid(row=0, column=0, sticky="w")
        self.cnt_radius = ttk.Entry(cnt_frame)
        self.cnt_radius.grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(cnt_frame, text="CNT Length (Å):").grid(row=1, column=0, sticky="w")
        self.cnt_length = ttk.Entry(cnt_frame)
        self.cnt_length.grid(row=1, column=1, padx=5, pady=2)
        
        # Graphene Parameters
        graphene_frame = ttk.LabelFrame(self.root, text="Graphene Parameters", padding=10)
        graphene_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(graphene_frame, text="Graphene Length (Å):").grid(row=0, column=0, sticky="w")
        self.graphene_length = ttk.Entry(graphene_frame)
        self.graphene_length.grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(graphene_frame, text="Graphene Width (Å):").grid(row=1, column=0, sticky="w")
        self.graphene_width = ttk.Entry(graphene_frame)
        self.graphene_width.grid(row=1, column=1, padx=5, pady=2)
        
        # Tilt Parameters
        tilt_frame = ttk.LabelFrame(self.root, text="Tilt Parameters", padding=10)
        tilt_frame.pack(fill="x", padx=10, pady=5)
        
        self.use_tilt = tk.BooleanVar(value=False)
        ttk.Checkbutton(tilt_frame, text="Use Tilted CNT", variable=self.use_tilt).grid(row=0, column=0, columnspan=2, sticky="w")
        
        ttk.Label(tilt_frame, text="Tilt Angle (degrees):").grid(row=1, column=0, sticky="w")
        self.tilt_angle = ttk.Entry(tilt_frame)
        self.tilt_angle.grid(row=1, column=1, padx=5, pady=2)
        self.tilt_angle.insert(0, "45")
        
    def create_buttons(self):
        button_frame = ttk.Frame(self.root)
        button_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(button_frame, text="Generate Structure", command=self.generate_structure).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear Output", command=self.clear_output).pack(side="left", padx=5)
        
    def create_output_area(self):
        output_frame = ttk.LabelFrame(self.root, text="Output", padding=10)
        output_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.output_text = tk.Text(output_frame, height=20, width=60)
        self.output_text.pack(fill="both", expand=True)
        
    def clear_output(self):
        self.output_text.delete(1.0, tk.END)
        
    def log_output(self, message):
        self.output_text.insert(tk.END, message + "\n")
        self.output_text.see(tk.END)
        
    def calculate_graphene_dimensions(self, length_angstrom, width_angstrom):
        # Convert physical dimensions to number of unit cells
        # Graphene unit cell size: a = 2.46 Å
        nx = int(np.ceil(length_angstrom / 2.46))
        ny = int(np.ceil(width_angstrom / 2.46))
        return nx, ny
        
    def calculate_cnt_dimensions(self, radius_angstrom, length_angstrom):
        # Convert physical dimensions to CNT indices
        # CNT radius = a * sqrt(n^2 + m^2) / (2*pi)
        # where a = 2.46 Å is the graphene lattice constant
        # For simplicity, we'll use (n,0) type CNT
        n = int(np.ceil(2 * np.pi * radius_angstrom / 2.46))
        # CNT length = p * 2.46 Å
        p = int(np.ceil(length_angstrom / 2.46))
        return n, p
        
    def generate_structure(self):
        try:
            # Get input parameters
            cnt_radius = float(self.cnt_radius.get())
            cnt_length = float(self.cnt_length.get())
            graphene_length = float(self.graphene_length.get())
            graphene_width = float(self.graphene_width.get())
            use_tilt = self.use_tilt.get()
            tilt_angle = float(self.tilt_angle.get()) if use_tilt else 0
            
            # Calculate dimensions
            nx, ny = self.calculate_graphene_dimensions(graphene_length, graphene_width)
            n, p = self.calculate_cnt_dimensions(cnt_radius, cnt_length)
            
            self.log_output(f"\nCalculated dimensions:")
            self.log_output(f"Graphene: {nx}x{ny} unit cells")
            self.log_output(f"CNT: (n,p) = ({n},{p})")
            
            # Create graphene sheets
            bottom_graphene = Graphene(nx=nx, ny=ny)
            top_graphene = Graphene(nx=nx, ny=ny)
            
            # Create CNT
            cnt = CarbonNanotube(n=n, p=p)
            
            # Get initial CNT dimensions
            cnt_min_before, cnt_max_before = cnt.get_bounding_box()
            actual_cnt_length = cnt_max_before[2] - cnt_min_before[2]
            actual_cnt_radius = cnt.r
            
            self.log_output(f"\nActual CNT dimensions:")
            self.log_output(f"Length: {actual_cnt_length:.2f} Å")
            self.log_output(f"Radius: {actual_cnt_radius:.2f} Å")
            
            # Create box
            vacuum = max(200.0, cnt_length * 2)  # Ensure enough vacuum
            box = Box.init_from_graphene(bottom_graphene, vacuum=vacuum)
            
            # Calculate graphene positions
            bottom_min, bottom_max = bottom_graphene.get_bounding_box()
            top_min, top_max = top_graphene.get_bounding_box()
            
            if use_tilt:
                # Calculate distance based on tilt angle
                distance = cnt_length * np.cos(np.radians(tilt_angle))
                x_offset = cnt_length * np.sin(np.radians(tilt_angle))
                
                # Position bottom graphene
                bottom_center = (bottom_min + bottom_max) / 2
                bottom_translation = -bottom_center.copy()
                bottom_translation[2] = -distance/2 - bottom_center[2]
                bottom_graphene.translate(bottom_translation)
                
                # Position top graphene
                top_center = (top_min + top_max) / 2
                top_translation = -top_center.copy()
                top_translation[0] = x_offset - top_center[0]
                top_translation[2] = distance/2 - top_center[2]
                top_graphene.translate(top_translation)
                
                # Position and rotate CNT
                cnt_min, cnt_max = cnt.get_bounding_box()
                cnt_center = (cnt_min + cnt_max) / 2
                cnt.translate(-cnt_center)
                
                # Rotate CNT
                rotation_axis = np.array([0.0, 1.0, 0.0])
                rotation_angle = np.radians(tilt_angle)
                rotation_center = np.array([0.0, 0.0, 0.0])
                cnt.rotate(rotation_axis, rotation_angle, rotation_center)
            else:
                # Position bottom graphene
                bottom_center = (bottom_min + bottom_max) / 2
                bottom_translation = -bottom_center.copy()
                bottom_translation[2] = -actual_cnt_length/2 - bottom_center[2]
                bottom_graphene.translate(bottom_translation)
                
                # Position top graphene
                top_center = (top_min + top_max) / 2
                top_translation = -top_center.copy()
                top_translation[2] = actual_cnt_length/2 - top_center[2]
                top_graphene.translate(top_translation)
                
                # Position CNT
                cnt_min, cnt_max = cnt.get_bounding_box()
                cnt_center = (cnt_min + cnt_max) / 2
                cnt.translate(-cnt_center)
            
            # Add structures to box
            box.add_cluster(bottom_graphene)
            box.add_cluster(top_graphene)
            box.add_cluster(cnt)
            
            # Delete atoms in CNT region
            box.delete_atoms_in_cnt(bottom_graphene, cnt)
            box.delete_atoms_in_cnt(top_graphene, cnt)
            
            # Delete CNT atoms in graphene planes
            bottom_min, bottom_max = bottom_graphene.get_bounding_box()
            top_min, top_max = top_graphene.get_bounding_box()
            
            bottom_z_min = bottom_min[2] - 0.1
            bottom_z_max = bottom_max[2] + 0.1
            top_z_min = top_min[2] - 0.1
            top_z_max = top_max[2] + 0.1
            
            box.delete_atoms_in_region(axis=2, min_val=bottom_z_min, max_val=bottom_z_max, cluster=cnt)
            box.delete_atoms_in_region(axis=2, min_val=top_z_min, max_val=top_z_max, cluster=cnt)
            
            # Generate output filename
            filename = f"sandwich_r{cnt_radius:.1f}_l{cnt_length:.1f}_g{graphene_length:.1f}x{graphene_width:.1f}"
            if use_tilt:
                filename += f"_t{tilt_angle:.1f}"
            filename += ".data"
            
            output_file = os.path.join(self.output_dir, filename)
            
            # Write output file
            writer = LAMMPSWriter(box)
            writer.write_data_file(output_file)
            
            # Calculate and display final dimensions
            bottom_min, bottom_max = bottom_graphene.get_bounding_box()
            top_min, top_max = top_graphene.get_bounding_box()
            cnt_min, cnt_max = cnt.get_bounding_box()
            
            self.log_output(f"\nFinal structure dimensions:")
            self.log_output(f"Bottom graphene Z range: {bottom_min[2]:.2f} to {bottom_max[2]:.2f} Å")
            self.log_output(f"Top graphene Z range: {top_min[2]:.2f} to {top_max[2]:.2f} Å")
            self.log_output(f"CNT Z range: {cnt_min[2]:.2f} to {cnt_max[2]:.2f} Å")
            self.log_output(f"Distance between graphene sheets: {top_min[2] - bottom_max[2]:.2f} Å")
            self.log_output(f"Total atoms: {len(box.get_all_atom_ids())}")
            self.log_output(f"\nOutput file: {output_file}")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.log_output(f"\nError: {str(e)}")
    
    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = SandwichStructureGenerator()
    app.run() 