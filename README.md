## Overview

This project implements a Density Functional Theory (DFT) simulation framework in Python. It numerically solves for the equilibrium density profile of hard spheres confined in a nanopore using self-consistent iterative methods (e.g., Picard or Newton schemes). This work represents the first step in a larger study of lithium nanofiltration, aiming to model ion transport and selectivity in nanoporous membranes from a fundamental, density-functional perspective. The framework provides a foundation for extending toward multi-component electrolyte systems and coupling with electrostatic or dynamic transport models in future research.

## Features

* **Modular Architecture**: Each physical operation (e.g., density update, correlation functional, convolution kernel) is implemented as a separate function.
* **Self-consistent Iteration Loop**: Iteratively updates densities until convergence within a specified tolerance.
* **Physical Subroutines**: Includes methods like:
    * `DDRC()` - compute weighted densities and derivatives 
    * `DFEXRC()` - compute excess free energy derivatives
    * `DPS()` - update potential and chemical potential terms
    * `EXTERN()` - compute external potential field
    * `WD()` - calculate weight functions for the fundamental measure theory (FMT)
    * `CP()` - calculate bulk chemical potential
* **Formatted Output**: Produces well-structured excel files (e.g., `Concentration_Profile.csv`, `Reduced_Density_Profile.csv`).
* **Configurable Parameters**: All key constants and numerical grid parameters are defined at the top of the script for easy modification.

## File Structure

```
model.py
# Main DFT program
output/
├── Concentration_Profile.csv
├── Concentration_Profile.png
├── Reduced_Density_Profile.csv
├── Reduced_Density_Profile.png
└── Output.dat
```

## Workflow

1.  **Initialize parameters**: Grid size, density bounds, chemical potentials, etc.
2.  **Call hierarchy**:
    ```
    main() -> CP() 
           -> DDRC() -> WD() 
           -> DPS()  -> EXTERN(), DFEXRC() 
    ```
3.  **Check convergence**: Evaluate maximum and cumulative errors until fall below the specified tolerance.
4.  **Write results**: Output formatted profiles for postprocessing.

## Example Usage

To run the simulation, execute the main script:

```bash
python model.py
```

Optional arguments can be modified inside the file (e.g., number of grid points `NP`, step size `DR`, tolerance `TOLE`).

## Output Format Example

### Concentration_Profile.txt 

```csv
"R","RH01","RH02","RH03"
"0.00","0.0001","0.0001","0.0001"
"0.01","0.0002","0.0002","0.0002"
...
```

### Reduced_Density_Profile.txt 

```
R         RHO_Reduced
0.00      0.010
0.01      0.012
...
```

## Dependencies

* Python 3.10+
* NumPy
* SciPy (for integration)
* Matplotlib (optional, for visualization)

Install dependencies:
```bash
pip install numpy scipy matplotlib
```

## Visualization

Example snippet for plotting results:

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Concentration_Profile.txt", skiprows=1)
R = data[:,0]
RHO = data[:,1:]

plt.plot(R, RHO[:, 0], label='Component 1')
plt.plot(R, RHO [:, 1], label='Component 2')
plt.xlabel("Radial Distance R")
plt.ylabel("Concentration")
plt.legend()
plt.show()
```

## Future Improvements

* [x] Parallelization with joblib or multiprocessing
* [x] Integration of advanced exchange-correlation functionals
* [ ] Adaptive grid refinement
* [ ] GUI or web interface for parameter tuning
