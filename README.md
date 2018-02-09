# mbpol_calc_ase
MBPol calculator for the atomic simulation environment (ASE).

## Installing dependencies
Assuming anaconda is installed, run the following commands in your terminal
```
conda install -c omnia openmm
conda config --add channels https://conda.anaconda.org/omnia
conda config --add channels https://conda.anaconda.org/paesanilab
conda install mbpol
conda install -c conda-forge ase
```

## How-to

MbpolCalculator can be used like any other ASE calculator (see their documentation). A quick example of how to initialize it is given in driver.py

### Restrictions 

- Only orthorhombic lattices supported (by OpenMM)
- get_stress() not implemented 
- Only vacuum boundary conditions or PBC in all euclid. directions supported 

