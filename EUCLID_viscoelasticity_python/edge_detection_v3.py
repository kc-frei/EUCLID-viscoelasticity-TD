# Automatic detection of the outer edges and holes in a 2D mesh from CSV data.
# This script uses Delaunay triangulation to identify edges and holes in a mesh.
# Parameters for hole size and edge margin can be adjusted, and results can be saved in separate or combined CSV files.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import os
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog

def edge_hole_detection(csv_path, min_hole_size=0.5, max_hole_size=20.0, combine_output=False):
   # Load data
   data = pd.read_csv(csv_path, skiprows=6, sep=';')
   node_ids = data.iloc[:, 0].values
   x = data.iloc[:, 1].values
   y = data.iloc[:, 2].values
   points = np.column_stack((x, y))
   
   # Get boundary nodes
   edge_node_ids = get_boundary_nodes(points, node_ids, margin=0.25)
   
   # Get hole nodes
   hole_node_ids = find_holes(points, node_ids, min_hole_size, max_hole_size)
   
   # Remove boundary nodes from hole nodes
   hole_node_ids = np.setdiff1d(hole_node_ids, edge_node_ids)
   
   # Save results
   output_dir = os.path.dirname(csv_path)
   
   if combine_output:
       combined_ids = np.concatenate([edge_node_ids, hole_node_ids])
       pd.DataFrame({'node_id': combined_ids}).to_csv(
           os.path.join(output_dir, 'combined_nodes.csv'), index=False)
   else:
       pd.DataFrame({'node_id': edge_node_ids}).to_csv(
           os.path.join(output_dir, 'edge_nodes.csv'), index=False)
       pd.DataFrame({'node_id': hole_node_ids}).to_csv(
           os.path.join(output_dir, 'hole_nodes.csv'), index=False)
   
   # Visualize results
   plt.figure(figsize=(12, 10))
   plt.scatter(x, y, s=1, color='lightblue')
   
   # Plot edge nodes
   edge_mask = np.isin(node_ids, edge_node_ids)
   plt.scatter(x[edge_mask], y[edge_mask], s=8, color='red', label='Edge Nodes')
   
   # Plot hole nodes
   hole_mask = np.isin(node_ids, hole_node_ids)
   plt.scatter(x[hole_mask], y[hole_mask], s=8, color='green', label='Hole Nodes')
   
   plt.title('Edge and Hole Detection')
   plt.xlabel('X (mm)')
   plt.ylabel('Y (mm)')
   plt.legend()
   plt.axis('equal')
   plt.grid(True)
   plt.savefig(os.path.join(output_dir, 'node_detection.png'), dpi=300)
   plt.show()
   
   print(f"Found {len(edge_node_ids)} edge nodes and {len(hole_node_ids)} hole nodes")
   
   return edge_node_ids, hole_node_ids

def get_boundary_nodes(points, node_ids, margin=0.2):
   x_min, x_max = np.min(points[:, 0]), np.max(points[:, 0])
   y_min, y_max = np.min(points[:, 1]), np.max(points[:, 1])
   
   boundary_mask = (
       (points[:, 0] <= x_min + margin) | 
       (points[:, 0] >= x_max - margin) | 
       (points[:, 1] <= y_min + margin) | 
       (points[:, 1] >= y_max - margin)
   )
   
   return node_ids[boundary_mask]

def find_holes(points, node_ids, min_hole_size=0.5, max_hole_size=20.0):
   # Create Delaunay triangulation
   tri = Delaunay(points)
   
   # Calculate circumradius for each triangle
   hole_candidates = []
   
   for i, simplex in enumerate(tri.simplices):
       # Get triangle vertices
       vertices = points[simplex]
       
       # Calculate lengths of triangle sides
       a = np.linalg.norm(vertices[1] - vertices[0])
       b = np.linalg.norm(vertices[2] - vertices[1])
       c = np.linalg.norm(vertices[0] - vertices[2])
       
       # Calculate semi-perimeter
       s = (a + b + c) / 2
       
       # Calculate area using Heron's formula
       area = np.sqrt(max(0, s * (s - a) * (s - b) * (s - c)))
       
       # Calculate circumradius
       if area > 0:
           circum_r = (a * b * c) / (4 * area)
           
           # Check if circumradius is within range
           if min_hole_size <= circum_r <= max_hole_size:
               hole_candidates.extend(simplex)
   
   # Get unique points
   hole_indices = np.unique(hole_candidates)
   return node_ids[hole_indices]

def main():
   # Create GUI
   root = tk.Tk()
   root.withdraw()
   
   # Select folder
   folder_path = filedialog.askdirectory(title="Select Folder with CSV Files")
   if not folder_path:
       print("No folder selected. Exiting.")
       return
   
   # Find all CSV files in the folder
   csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
   if not csv_files:
       messagebox.showerror("Error", "No CSV files found in the selected folder.")
       return
   
   # Get first CSV file
   csv_path = os.path.join(folder_path, csv_files[0])
   
   # Get parameters
   min_hole_size = simpledialog.askfloat("Input", 
                                        "Enter minimum hole size:", 
                                        initialvalue=0.5, 
                                        minvalue=0.1, 
                                        maxvalue=20.0)
   if min_hole_size is None:
       min_hole_size = 0.5
       
   max_hole_size = simpledialog.askfloat("Input", 
                                        "Enter maximum hole size:", 
                                        initialvalue=10.0, 
                                        minvalue=min_hole_size, 
                                        maxvalue=100.0)
   if max_hole_size is None:
       max_hole_size = 20.0
   
   # Ask about combined output
   combine_output = messagebox.askyesno("Output Option", 
                                       "Combine edge and hole nodes in one CSV file?")
   
   # Run detection
   edge_node_ids, hole_node_ids = edge_hole_detection(csv_path, min_hole_size, max_hole_size, combine_output)
   
   messagebox.showinfo("Complete", 
                      f"Processing complete!\nFound {len(edge_node_ids)} edge nodes and {len(hole_node_ids)} hole nodes")

if __name__ == "__main__":
   main()