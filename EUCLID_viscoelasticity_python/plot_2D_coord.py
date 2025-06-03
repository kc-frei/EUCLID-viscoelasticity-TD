# Visualization of the DIC coordinates from a CSV file using pandas and matplotlib

import pandas as pd
import matplotlib.pyplot as plt

# Load with specific column separator
data = pd.read_csv('Creep_Perforated_WJ_3holes2_0001.csv', 
                 skiprows=6, 
                 sep=';')

# Print column names for debugging
print("Column names:", data.columns.tolist())

# Extract coordinates
x = data.iloc[:, 1].values  # Second column (x coordinates)
y = data.iloc[:, 2].values  # Third column (y coordinates)

# Create plot
plt.figure(figsize=(10, 8))
plt.scatter(x, y, s=1)
plt.title('2D Coordinate Plot')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.axis('equal')
plt.grid(True)
plt.savefig('coordinates_plot.png', dpi=300)
plt.show()