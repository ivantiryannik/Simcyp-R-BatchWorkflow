import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D

# Path to the folder containing matrices
folder_path = r'-------'

# Get a list of all files in the folder
file_list = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# Sort the file list based on the matrix number (assuming the file names are in the format "matrix_X.csv")
file_list.sort(key=lambda x: int(x.split('_')[1].split('.')[0]))

# Define the number of rows and columns for the subplot grid
num_rows = 5  
num_cols = 5  

# Create a figure and subplots
fig = plt.figure(figsize=(20, 18))

#Define the indices of the matrices to be visible
visible_indices = {0, 6, 12, 18, 24}  # 0-based indexing

# Loop through each matrix file and create a subplot
for i, file_name in enumerate(file_list):
    # Load the matrix from the CSV file
    matrix = np.genfromtxt(os.path.join(folder_path, file_name), delimiter=',', skip_header=1, filling_values=0)

    # Define custom colors based on values
    colors = np.empty_like(matrix, dtype=str)
    colors[matrix < 0.8] = 'yellow'
    colors[(matrix >= 0.8) & (matrix <= 1.25)] = 'green'
    colors[matrix > 1.25] = 'red'

    # Create meshgrid from matrix indices
    x, y = np.meshgrid(np.arange(matrix.shape[0]), np.arange(matrix.shape[1]))

    # Calculate subplot position dynamically (bottom-up order)
    row_index = i % num_cols
    col_index = i // num_cols

    # Add a subplot to the figure at the specified position
    ax = fig.add_subplot(num_rows, num_cols, (num_rows - 1 - row_index) * num_cols + col_index + 1, projection='3d')
    
    # Add surface
    surf = ax.plot_surface(x, y, matrix, facecolors=colors, cmap='viridis', rcount=100, ccount=100, alpha=0.8)

    # Set z-axis limits
    ax.set_zlim(0.2, 2.0)

    # Set labels
    # ax.set_xlabel('Ki')
    # ax.set_ylabel('Ka')
    ax.zaxis.set_rotate_label(False)
    # ax.set_zlabel('IMDRA', rotation=90)
    
     # Position x and y axes closer to the graph
    ax.xaxis.labelpad = -5  
    ax.yaxis.labelpad = -5  

    # Set custom x and y axis values
    custom_xticks = [0.015, 0.03, 0.06, 0.12, 0.24, 0.48, 0.96, 1.02, 3.84, 7.68]  # Replace with your desired x-axis values
    custom_yticks = [0.05, 0.1, 0.2, 0.4, 0.8, 1, 2, 3, 4, 6]   # Replace with your desired y-axis values
    ax.set_xticks(np.arange(len(custom_xticks)))
    ax.set_yticks(np.arange(len(custom_yticks)))
    ax.set_xticklabels(custom_xticks)
    ax.set_yticklabels(custom_yticks)
    ax.tick_params(axis='x', rotation=-20, labelsize=4)  # Adjust the rotation angle and label size
    ax.tick_params(axis='y', rotation=20, labelsize=4)  # Adjust the rotation angle and label size
    # ax.tick_params(axis='both', which='major', labelsize=8)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    ax.set_zticklabels([])
    
    # Set rotation angles (elevation, azimuthal)
    ax.view_init(elev=32, azim=53)
    
    plt.subplots_adjust(wspace=0, hspace=0)
    
    # Make specific subplots invisible except the ones in visible_indices
    #if i not in visible_indices:
        #ax.set_visible(False)

# Adjust layout to prevent overlapping
plt.tight_layout()

# X Axis
fig.text(0.2, 0.09, '0.1', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.355, 0.09, '0.3', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.51, 0.09, '0.5', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.665, 0.09, '0.7', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.82, 0.09, '1', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.51, 0.05, '(Substrate) Fg', ha='center', va='center', fontsize=20, color='Black')

# Y Axis
fig.text(0.10, 0.19, '0.05', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.10, 0.345, '0.25', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.10, 0.50, '1.25', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.10, 0.65, '6.25', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.10, 0.80, '31.25', ha='center', va='center', fontsize=16, color='Black')
fig.text(0.05, 0.50, '(Inhibitor) Vss', ha='center', va='center', fontsize=20, color='Black', rotation=90)

# Create custom legend handles
legend_colors = ['red', 'green', 'yellow']
legend_labels = ['> 1.25', '0.8 - 1.25', '< 0.8']
handles = [mpatches.Patch(color=legend_colors[i], label=legend_labels[i]) for i in range(len(legend_colors))]

# Add the legend
fig.legend(handles=handles, loc='upper right', title='Discrepancy Intervals', frameon=False, fontsize='16', title_fontsize='16')

#fig.text(0.1, 0.9, 'b', ha='left', fontsize=30, color='black', fontweight='bold')

# Export Figure
plt.savefig('X', dpi=600)

# Show the plot
plt.show()

