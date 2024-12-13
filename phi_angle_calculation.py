# Calculating the Phi_Angle

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles
import matplotlib.pyplot as plt
import numpy as np

def get_phi_atoms(residue):
    prev_residue_idx = residue.resindices[0] - 1
    C_prev = u.select_atoms(f"resid {prev_residue_idx} and name C").positions
    N = residue.atoms.select_atoms("name N").positions
    CA = residue.atoms.select_atoms("name CA").positions
    return C_prev, N, CA

# Load your protein structure and trajectory
u = mda.Universe('your_file.gro', 'your_file.xtc')

# Select the residue for which you want to calculate the phi angle
residue = u.select_atoms('resid 00')

# Initialize a list to store phi angle data
phi_data = []

# Iterate over the frames in the trajectory
frame_indices = []
for i, ts in enumerate(u.trajectory):
    # Calculate the phi angle
    phi_atoms = get_phi_atoms(residue)
    if len(phi_atoms[0]) > 0:
        phi_angle = np.rad2deg(calc_angles(*phi_atoms))
        # Append the angle to the data list
        phi_data.append(phi_angle)
        frame_indices.append(i)

# Convert the lists to arrays
phi_data = np.array(phi_data)
frame_indices = np.array(frame_indices)

# Save the frame indices and phi angles to a .dat file
plot_data = np.column_stack((frame_indices, phi_data))
np.savetxt('phi_angles_plot_data_resname.dat', plot_data)

# Plot the phi angles
plt.figure()
plt.plot(frame_indices, phi_data)
plt.xlabel('Frame')
plt.ylabel('Phi Angle (degrees)')
plt.title('Phi Angle of RES00')
plt.show()

