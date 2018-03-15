# mbpol_calculator
MB-pol calculator for the atomic simulation environment (ASE). The calculator implements an interface to OpenMM from which the MB-pol plugin can be used.

## Installation

### Dependencies
Assuming anaconda is installed, run the following commands in your terminal
```
conda install -c omnia openmm
conda config --add channels https://conda.anaconda.org/omnia
conda config --add channels https://conda.anaconda.org/paesanilab
conda install mbpol
conda install -c conda-forge ase
```
### Install library
Change into the directory that contains ```setup.py``` and execute
```
pip install -e .
```


## How-to

MbpolCalculator can be used like any other ASE calculator (see their documentation). A quick example is given in ```driver.py```.

### Restrictions 

- Only orthorhombic lattices supported (by OpenMM)
- ```get_stress()``` not implemented 
- Only vacuum boundary conditions or PBC in all euclid. directions supported 

