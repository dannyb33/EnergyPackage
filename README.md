energy_package
==============================
[//]: # (Badges)
[![CI](https://github.com/dannyb33/EnergyPackage/actions/workflows/CI.yaml/badge.svg)](https://github.com/dannyb33/EnergyPackage/actions/workflows/CI.yaml)
[![codecov](https://codecov.io/gh/dannyb33/EnergyPackage/graph/badge.svg?token=4M1TFSOPJE)](https://codecov.io/gh/dannyb33/EnergyPackage)


# Ising Model Simulation Package

A modular Python package for simulating bitstring-based configurations and interactions using custom BitString and IsingHamiltonian classes.

Supports both exact and Monte Carlo simulation methods to compute average physical quantities such as energy, magnetization, heat capacity, and magnetic susceptibility.

## Features:

- Exact thermodynamic computation across all states

- Monte Carlo sampling for approximate calculations

- Graph-based systems via NetworkX

- Efficient state updates and delta E calculations

- Unit tests using Pytest

## Installation: 
### Run: 
pip install -e . 

### Required packages: 
numpy, networkx, pytest

## Usage Instructions:

### 1. Import Packages
```python
import energy_package as ep
import networkx as nx
import numpy as np
```

### 2. Create Graph and Hamiltonian
```python
G = nx.Graph()
G.add_nodes_from(range(6))
G.add_edges_from([(i, (i+1)%6) for i in range(6)])
for edge in G.edges:
    G.edges[edge]['weight'] = 2.0
ham = ep.IsingHamiltonian(G)
ham.set_mu(np.array([0.1] * 6))
```

### 3. Create and Configure BitString
```python
bs = ep.BitString(6)
bs.set_integer_config(3)
bs.flip_site(2)
```

### 4. Compute Exact Thermodynamic Properties 
```python
E, M, HC, MS = ham.compute_average_values(T=1.0)
```

### 5. Run Monte Carlo Simulation
```python
mc = ep.MonteCarlo(ham)
E_list, M_list = mc.run(T=1.0, n_samples=1000, n_burn=500)
```

### Utility Functions
```python
bs.on() # Number of 1s 
bs.off() # Number of 0s 
bs.integer() # Convert bitstring to integer 
ham.energy(bs) 
ham.magnetization(bs) 
ham.delta_e(bs, i) # Energy difference for flipping site i
```

### Testing: 
#### Run: 
pytest

## License & Contributions: 
Feel free to fork!

## Copyright

Copyright (c) 2025, Daniel Benish

## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.
