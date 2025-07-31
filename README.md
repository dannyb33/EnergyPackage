energy_package
==============================
[//]: # (Badges)
[![CI](https://github.com/dannyb33/EnergyPackage/actions/workflows/CI.yaml/badge.svg)](https://github.com/dannyb33/EnergyPackage/actions/workflows/CI.yaml)
[![codecov](https://codecov.io/gh/dannyb33/EnergyPackage/graph/badge.svg?token=4M1TFSOPJE)](https://codecov.io/gh/dannyb33/EnergyPackage)


# Ising Model Simulation Package

**Technologies:** Python, NumPy, NetworkX, Simulation & Modeling

---

## Overview

Built a modular Python package to simulate bitstring-based configurations and interactions using custom `BitString` and `IsingHamiltonian` classes.

Implemented both exact and Monte Carlo algorithms to estimate average values like energy and magnetization of systems at varying temperatures.

Used NetworkX to define and analyze graph-based systems, with optimized update methods for efficient simulation performance.

---

## Features

- Custom bitstring configuration handling with `BitString` class
- Hamiltonian calculations for Ising models via `IsingHamiltonian` class
- Exact and Monte Carlo methods for statistical estimation
- Graph-based system modeling using NetworkX
- Performance optimizations for simulation updates

---

## Usage Instructions

1. **Create a BitString configuration**

```python
import energy_package as ep

bitstr = ep.BitString(5)
bitstr.set_integer_config(3)

```

2. **Create a Network Graph**

```python
import networkx as nx

N = 6
Jval = 2.0
G = nx.Graph()
G.add_nodes_from(range(N))
G.add_edges_from([(i, (i + 1) % N) for i in range(N)])
for edge in G.edges:
    G.edges[edge]['weight'] = Jval

```

3. **Initialize Ising Hamiltonian and set external field**

```python
ham = ep.IsingHamiltonian(G)
ham.set_mu([0.1] * N)

```

4. **Calculate energy and magnetization**

```python
bitstr.flip_site(2)
bitstr.flip_site(3)

energy = ham.energy(bitstr)
magnetization = ham.magnetization(bitstr)

```

5. **Compute average values at a given temperature**

```python
E, M, HC, MS = ham.compute_average_values(1) # Temperature = 1
print(f"Energy: {E}, Magnetization: {M}")

```

### Copyright

Copyright (c) 2025, Daniel Benish

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.
