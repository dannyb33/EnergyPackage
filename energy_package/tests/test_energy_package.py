"""
Unit and regression test for the energy_package package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import energy_package as ep
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt



def test_energy_package_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "energy_package" in sys.modules
    
def test_bitstr():
    """Test functions of bitstring class"""
    bit1 = ep.BitString(5)
    bit2 = ep.BitString(5)
    assert(len(bit1) == 5)
    assert(str(bit1) == "00000")
    assert(bit1.__eq__(bit2))
    
    bit1.set_integer_config(3)
    bit2.set_config([0,0,0,1,1])
    assert(bit1.__eq__(bit2))
    
    assert(bit1.integer == 3)
    
    
def test_energy():
    """Test energy() Function"""
    
    #Create a graph definining ising interactions
    N = 6
    Jval = 2.0
    G = nx.Graph()
    G.add_nodes_from([i for i in range(N)])
    G.add_edges_from([(i,(i+1)% G.number_of_nodes() ) for i in range(N)])
    for e in G.edges:
        G.edges[e]['weight'] = Jval
    
    conf = ep.BitString(N)
    ham = ep.IsingHamiltonian(G)
    ham.set_mu(np.array([.1 for i in range(N)]))

    conf.flip_site(2)
    conf.flip_site(3)
    e = ham.energy(conf)
    assert(np.isclose(e, 3.8))

def test_average():
    """Test compute_average_values() function"""
    
    #Create a graph definining ising interactions
    N = 6
    Jval = 2.0
    G = nx.Graph()
    G.add_nodes_from([i for i in range(N)])
    G.add_edges_from([(i,(i+1)% G.number_of_nodes() ) for i in range(N)])
    for e in G.edges:
        G.edges[e]['weight'] = Jval
    
    # Define a new configuration instance for a 6-site lattice
    conf = ep.BitString(N)
    ham = ep.IsingHamiltonian(G)

    # Compute the average values for Temperature = 1
    E, M, HC, MS = ham.compute_average_values(1)

    assert(np.isclose(E,  -11.95991923))
    assert(np.isclose(M,   -0.00000000))
    assert(np.isclose(HC,   0.31925472))
    assert(np.isclose(MS,   0.01202961))
    
    

