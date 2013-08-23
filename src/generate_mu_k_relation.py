"""
This script takes a binding energy landscape, stored as a Python
array, and generates a table relating copy number to chemical
potential.
"""
import sys
from math import exp
from chem_pot_utils import load_array
sys.path.append("lib/utils")
from utils import data2csv
from params import *

mus = range(-5,50)

if __name__ == "__main__":
    arr_fname = sys.argv[1]
    out_fname = sys.argv[2]
    approximation = len(sys.argv) == 4 and sys.argv[3] == "approximation"
    binding_energies = load_array(arr_fname,'f')
    mu_k_dict = {}
    for mu in mus:
        print "mu:",mu
        if approximation:
            Z = sum(exp(-beta*ep) for ep in binding_energies)
            mu_k_dict[mu] = exp(beta*mu)*Z
        else:
            probs = [1/(1+exp(beta*(ep-mu))) for ep in binding_energies]
            copy_number = sum(probs)
            mu_k_dict[mu] = copy_number
    data2csv([[mu,mu_k_dict[mu]] for mu in mus],header=["mu","k"],filename=out_fname)
        
