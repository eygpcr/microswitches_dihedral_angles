# Micro-switches Dihedral Angles

## Project Description
This project focuses on analyzing the angular changes occurring in the micro-switch regions of GPCRs based on the common activation mechanism of class A GPCRs. It allows the calculation of **phi**, **psi**, and **chi** dihedral angles throughout molecular dynamics simulations. The code is designed to compare angular changes before and after ligand binding in both active and inactive structures of GPCRs.

---

## Essential Tools
- **Python**
- **MDAnalysis**
- **Google Colab**

---

## Installation Steps
1. Prepare your **GROMACS** trajectory files in `.gro` and `.xtc` formats.
2. Clone the project or download the scripts.
3. Install the required Python libraries:
   ```bash
   pip install MDAnalysis matplotlib numpy

---

## Features
- Analyze angular changes in **PIF**, **E/DRY**, **CWxP**, and **NPxxY** micro-switch regions.
- Compare GPCR activation states before and after ligand binding.
- Generate graphical outputs for the **phi**, **psi**, and **chi1** angles over the course of a molecular dynamics simulation.

---

## Usage

### Calculating the Phi Angle
To calculate the **phi angle**, use the following resources:
- **Google Colab Notebook:** [Open in Colab](https://colab.research.google.com/github/eygpcr/microswitches-dihedral-angles/blob/main/phi_angle_calculation.ipynb)
- **Python Script:** [Download phi_angle_calculation.py](https://github.com/eygpcr/microswitches-dihedral-angles/raw/main/phi_angle_calculation.py)

### Calculating the Psi Angle
To calculate the **psi angle**, use the following resources:
- **Google Colab Notebook:** [Open in Colab](https://colab.research.google.com/github/eygpcr/microswitches-dihedral-angles/blob/main/psi_angle_calculation.ipynb)
- **Python Script:** [Download psi_angle_calculation.py](https://github.com/eygpcr/microswitches-dihedral-angles/raw/main/psi_angle_calculation.py)

### Calculating the Chi Angle
To calculate the **chi angle**, use the following resources:
- **Google Colab Notebook:** [Open in Colab](https://colab.research.google.com/github/eygpcr/microswitches-dihedral-angles/blob/main/chi_angle_calculation.ipynb)
- **Python Script:** [Download chi_angle_calculation.py](https://github.com/eygpcr/microswitches-dihedral-angles/raw/main/chi_angle_calculation.py)

---

## Contribution

We welcome contributions to this project! Here's how you can contribute:

1. **Fork** the repository.
2. **Clone** your fork to your local machine:
   ```bash
   git clone https://github.com/eygpcr/microswitches-dihedral-angles.git

---

## License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute this project under the terms of the license.
