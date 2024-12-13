# Calculating the Psi_Angle

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles
import matplotlib.pyplot as plt
import numpy as np

def get_psi_atoms(residue):
    next_residue_idx = residue.resindices[0] + 1
    N_next = u.select_atoms(f"resid {next_residue_idx} and name N").positions
    CA = residue.atoms.select_atoms("name CA").positions
    C = residue.atoms.select_atoms("name C").positions
    return CA, C, N_next

# Load your protein structure and trajectory
u = mda.Universe('your_file.gro', 'your_file.xtc')

# Select the residue for which you want to calculate the psi angle
residue = u.select_atoms('resid 00')

# Initialize a list to store psi angle data
psi_data = []

# Iterate over the frames in the trajectory
frame_indices = []
for i, ts in enumerate(u.trajectory):
    # Calculate the psi angle
    psi_atoms = get_psi_atoms(residue)
    if len(psi_atoms[2]) > 0:
        psi_angle = np.rad2deg(calc_angles(*psi_atoms))
        # Append the angle to the data list
        psi_data.append(psi_angle)
        frame_indices.append(i)

# Convert the lists to arrays
psi_data = np.array(psi_data)
frame_indices = np.array(frame_indices)

# Save the frame indices and psi angles to a .dat file
plot_data = np.column_stack((frame_indices, psi_data))
np.savetxt('psi_angles_plot_data_resname.dat', plot_data)

# Plot the psi angles
plt.figure()
plt.plot(frame_indices, psi_data)
plt.xlabel('Frame')
plt.ylabel('Psi Angle (degrees)')
plt.title('Psi Angle of RES00')
plt.show()


