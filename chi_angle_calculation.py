# Calculating the Chi1_Angle
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import matplotlib.pyplot as plt
import numpy as np

def get_chi1_atoms(residue):
    N = residue.atoms.select_atoms("name N").positions
    CA = residue.atoms.select_atoms("name CA").positions
    CB = residue.atoms.select_atoms("name CB").positions
    CG = residue.atoms.select_atoms("name CG").positions
    return N, CA, CB, CG

# Load your protein structure and trajectory
u = mda.Universe('your_file.gro', 'your_file.xtc')

# Select the residue for which you want to calculate the chi1 angle
residue = u.select_atoms('resid 00')

# Initialize a list to store chi1 angle data
chi1_data = []

# Iterate over the frames in the trajectory
frame_indices = []
for i, ts in enumerate(u.trajectory):
    # Calculate the chi1 angle
    chi1_atoms = get_chi1_atoms(residue)
    
    # Make sure all required atoms are present
    if all(len(atom_group) > 0 for atom_group in chi1_atoms):
        chi1_angle = np.rad2deg(calc_dihedrals(*chi1_atoms))
        # Append the angle to the data list
        chi1_data.append(chi1_angle)
        frame_indices.append(i)

# Convert the lists to arrays
chi1_data = np.array(chi1_data)
frame_indices = np.array(frame_indices)

# Save the frame indices and chi1 angles to a .dat file
plot_data = np.column_stack((frame_indices, chi1_data))
np.savetxt('chi1_angles_plot_data_resname.dat', plot_data)

# Plot the chi1 angles
plt.figure()
plt.plot(frame_indices, chi1_data)
plt.xlabel('Frame')
plt.ylabel('Chi1 Angle (degrees)')
plt.title('Chi1 Angle of RES00')
plt.show()


