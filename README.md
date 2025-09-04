# pyIBM
A user-friendly Python code for IBM/IBFM-1 calculations

The code is divided into the following packages:

1. pyIBM: A stand-alone package with all CG coefficients up to $N$=10 stored in a library. One can do calculations with user-defined Hamiltonian parameters, optimize the parameters by providing a list of energy levels, or study the sensitivity of the solution to the choice of parameters by performing a 2D contour plot.

2. generate_basis: A code package coupled to an executive file of a Fortran code which generates CG coefficients. The package first tries to read from the CG library, if the CG coefficients are not found, it calls the Fortran CG code to generate the coefficients and updates the library.

## pyIBM
This package provides three main functions: 

- `Energy_Bvalues.py`
- `fitting_parameters.py` 
- `contour.py`

More details about these functions are listed below

### Energy_Bvalues.py

Run in order to calculate energy levels, reduced transition probabilities (B values) and B ratios.

In order to run make sure to update:
- Variables in the "variables_library directory" directory
- Experimental data in the "experimental_data_library" directory
- If needed, update read_paths in `functions/read.py`

Run with command `python Energy_Bvalues.py`
Input name of element and mass number that matches the corresponding folders in "variables_library", e.g. "Te" and "120".

Spins should be given as angular momentum above the band head (such as 2 for a 15/2 state above the 11/2 band head)

Edit `Energy_Bvalues.py` to customize "settings" or output. Selected B values can be obtained with the following command:

`print(calculate_BE2(eigenvectors["A_B"], eigenvectors["C_D"], Li, Lf, N, cg, chi))`

where A and C are the spin of the initial and final states, respectively. B and D are the orders of the state (1 = yrast).
Li and Lf are the inital and final state angular momentum above the spin of the band head.
Ex: `calculate_BE2(eigenvectors["19/2_1"], eigenvectors["15/2_1"], 4, 2, N, cg, chi)` if the band head is 11/2.


### fitting_parameters.py

Run in order to fit the coefficients of the Hamiltonian, in order to optimize the RMSE between the
calculated and experimental energy levels.

In order to run make sure to update:
- "Settings" in fitting_parameters.py
- Optional to change the parameters of the algorithm in fitting_parameters.py in order to optimize the fitting
- Experimental data in the "experimental_data_library" directory
- If needed, update read_paths in functions/read.py

Run with command "python fitting_parameters.py" after editing the file to select nucleus and parameters

### contour.py

Creates a 2-D contour plot.
 
In order to run make sure to update:
- "Settings" in contour.py
- Variables in the "variables_library directory" directory (only boson number N is needed)
- If needed, update read_paths in functions/read.py
- Experimental data in the "experimental_data_library" directory (if RMSE function is selected)

Run with command `python contour.py` after editing the file to select nucleus and parameters

## create_basis
### creat_cg_coefficients.py

`create_cg_coefficients.py` is a script designed to calculate Clebsch-Gordan coefficients needed in
the calculations of the IBM. The code calls a Fortran script in SU3cgvcs to do the calculations and saves the output in the basis folder

to run type `python create_cg_coefficients.py` and input the boson number.

Copy the results in the basis folder to the pyIBM folder, or change the path to save directly in the basis folder of pyIBM, to use the coefficients in calculations. 



