import numpy as np
import math      
import copy as cp 
import networkx as nx
import random as rand

class BitString:
    """
    Simple class to implement a config of bits
    """
    
    def __init__(self, N):
        """
        Set bitstring length, default all values to 0
        
        Parameters
        ----------
        N    : int
            initial bitstring length
        """
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        return str(self.config)

    def __eq__(self, other):      
        """
        Return true if bitstrings are identical, false otherwise
        
        Parameters
        ----------
        other    : BitString
            Bitstring to compare
        
        Returns
        -------
        eq : boolean
            true if bitstrings are identical
        """  
        return all(self.config == other.config)
    
    def __len__(self):
        """
        Return length of bitstring
        
        Returns
        -------
        len : int
            length of bitstring
        """
        return len(self.config)

    def on(self):
        """
        Return number of bits that are on
        
        Returns
        -------
        sum : int
            number of on bits
        """
        sum = 0
        for i in self.config:
            if i == 1:
                sum += 1
        
        return sum
        

    def off(self):
        """
        Return number of bits that are off
        
        Returns
        -------
        sum : int
            number of off bits
        """
        sum = 0
        for i in self.config:
            if i == 0:
                sum += 1
        
        return sum

    def flip_site(self,i):
        """
        Flip the bit at site i
        
        Parameters
        ----------
        i    : int
            index of site to flip
        """
        self.config[i] = -1*(self.config[i] - 1)
        return self
    
    def integer(self):
        """
        Return the decimal integer corresponding to BitString
        
        Returns
        -------
        sum : int
            integer corresponding to bitstring
        """
        exp = 0
        sum = 0
        for i in range(self.config.size-1, 0, -1):
            sum += self.config[i] * 2**exp
            exp += 1
        
        return sum
 
    def set_config(self, s:list[int]):
        """
        Set the config from a list of integers
        
        Parameters
        ----------
        s    : list[int]
            input list of ints
        """
        self.config = np.array(s)
        return self

    def set_integer_config(self, dec:int):
        """
        convert a decimal integer to binary
    
        Parameters
        ----------
        dec    : int
            input integer
        """
             
        temp = dec
        rem = 0
        arr = []
        
        self.config = np.zeros(len(self.config), dtype=int)
        
        while temp != 0:
            rem = temp % 2
            temp = temp // 2
            arr.append(rem)
        
        arr.reverse()
        
        for i in range(len(arr)):
            self.config[len(self.config)-len(arr)+i] = arr[i]
        
        return self
            
            
class IsingHamiltonian:
    """
    Class to implement ising hamiltonian functions, compute energy, various average values
    """
    
    def __init__(self, G : nx.Graph):
        """
        Set hamiltonian graph as G, initialize mu array with 0s

        Parameters
        ----------
        G   : nx.Graph
            hamiltonian graph
        """
        self.G = G
        self.mu = np.array([0 for i in range(len(G))])
        self.N = len(self.mu)
    
    def energy(self, bs: BitString):
        """
        Compute energy of configuration, `bs`

        Parameters
        ----------
        bs   : Bitstring
            input configuration
            
        Returns
        -------
        sum  : float
            Energy of the input configuration
        """
        sum = 0
        a_array = nx.to_numpy_array(self.G)
            
        for i in range(len(a_array)):
            for j in range(i+1, len(a_array[i])):
                sum += a_array[i][j] * (bs.config[i]*-2 + 1) * (bs.config[j]*-2 + 1)
        
        for i in range(len(self.mu)):
            sum += self.mu[i]*-(bs.config[i]*-2+1)
        
        return sum
      
    def magnetization(self, bs: BitString):
        """
        Return magnetization of bitstring

        Parameters
        ----------
        bs   : Bitstring
            input configuration
            
        Returns
        -------
        magnetization  : int
            Magnetization of the input configuration
        """
        return (bs.on() - bs.off())
    
    def set_mu(self, mus:np.array):
        """
        Set mu array of hamiltonian
    
        Parameters
        ----------
        mus   : np.array
            array of mus
        """
        self.mu = mus
        
        return self

    
    def compute_average_values(self, T: float):
        """
        Compute average energy, magnetization, heat capacity, and
        magnetic susceptibility of hamiltonian at a given temp
    
        Parameters
        ----------
        T   : float
            temperature to compute at
            
        Returns
        -------
        E  : float
            Average energy of the hamiltonian
        M  : float
            Average magnetization of the hamiltonian
        HC  : float
            Average heat capacity of the hamiltonian
        MS  : float
            Average magnetic susceptibility of the hamiltonian
        """
        
        bs = BitString(len(self.G))

        E = 0.0
        M = 0.0
        EE = 0.0
        MM = 0.0 
        Z = 0.0
        
        HC = 0.0
        MS = 0.0
        
        k = 1
        beta = 1/(k*T)
 
        for i in range((2**len(bs))):
            bs.set_integer_config(i)
            curr_e = self.energy(bs)
            
            temp = (np.exp(-beta*curr_e))
            
            Z += temp
            
            E += curr_e*temp
            M += (bs.on() - bs.off())*temp
            EE += (curr_e**2)*temp
            MM += ((bs.on() - bs.off())**2)*temp
            
        E = E / Z
        M = M / Z
        EE = EE / Z
        MM = MM / Z
        
        HC = (EE-E**2)*(T**-2)
        MS = (MM-M**2)*(T**-1)
        
        return E, M, HC, MS
    
    def delta_e(self, bs:BitString, change_index:int):
        """
        Compute difference in energy of changine one index of a bitstring
    
        Parameters
        ----------
        bs  : BitString
            initial bitstring
        change_index : int
            index to change for potential bitstring
            
        Returns
        -------
        delta: float
            difference in energy between the two bitstrings
        """
        
        b = bs.config * -2 + 1
        delta = 0
        
        for neighbor in self.G.neighbors(change_index):
            delta += -2 * b[change_index] * b[neighbor] * self.G.edges[(change_index, neighbor)]["weight"]
        
        delta += -2 * self.mu[change_index] * b[change_index]
        
        return delta
    
class MonteCarlo:
    """
    Class to implement montecarlo functions to compute average values of energy and magnetism
    """
    
    def __init__(self, ham:IsingHamiltonian):
        """
        Assign hamiltonian for use

        Parameters
        ----------
        ham   : IsingHamiltonian
            desired hamiltonian
        """
        self.ham = ham
        
    def run(self, T:float, n_samples:int, n_burn:int):
        """
        Initialize configuration, i 
        Loop over Monte Carlo steps	    
            Loop over sites, n
                Propose new configuration, j, by flipping site, n.
                Compute flipping probability, W(i→j). 
                If  W(i→j) is greater than a randomly chosen number between 0 and 1, 
                    Accept (i = j), 
                else: 
                    Reject 
            Update average values with updated i
        
        Parameters
        ----------
        T   : float
            temperature to compute at
        n_samples: int
            number of samples to run
        n_burn: int
            nuber of samples to burn before saving measurements
            
        Returns
        -------
        E  : list
            list of average energy values
        M  : list
            list of average magnetization values
        """    
        bs = BitString(len(self.ham.G))
        E = []
        M = []
        
        for i in range(n_samples):            
            for j in range(len(bs)):
                delta = self.ham.delta_e(bs, j)
                                
                if delta <= 0 or np.exp(-delta / T) > rand.random():
                    bs.flip_site(j)
            
            if i >= n_burn:
                E.append(self.ham.energy(bs))
                M.append(self.ham.magnetization(bs))
        
        return E, M
                                
    
    def flip_prob(self, T:float, i_en:float, j_en:float):
        """
        Calculate probability of moving from state with i energy to j energy at T temp
        Deprecated, using delta energy for run function now
        
        Parameters
        ----------
        T   : float
            temperature to compute at
        i_en: float
            energy of initial configuration
        j_en: float
            energy of desired configuration
                
        Returns
        -------
        probability  : float
            probability of moving to desired configuration
        """
        if i_en >= j_en:
            return 1.0
        else:
            return np.exp(-(j_en - i_en) / T)
                    
                    

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print("Energy Package")
