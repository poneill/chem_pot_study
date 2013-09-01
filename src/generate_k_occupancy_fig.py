"""This script generates k-occupancy figs"""

import sys,os
sys.path.append("data/motifs")
from motifs import Escherichia_coli
sys.path.append("data")
from copy_numbers import copy_numbers
from chem_pot_utils import *
sys.path.append("src/sufficache")
from sufficache import PSSM
sys.path.append("src/utils")
from utils import transpose
from math import exp
from params import *
from matplotlib import pyplot as plt

if __name__ == "__main__":
    (mu_k_file,
     mu_k_approximate_file,
     outfile) = sys.argv[1:]
    folder,filename = os.path.split(mu_k_file)
    tf_name = filename[:filename.index("_")]
    copy_number = copy_numbers[tf_name] if tf_name in copy_numbers else None
    mus_ks = read_mu_k_table(tf_name,approx=False)
    mus,ks = transpose(mus_ks)
    approximate_mus_ks = read_mu_k_table(tf_name,approx=True)
    approx_mus,approx_ks = transpose(approximate_mus_ks)
    site_energies = get_regulon_binding_energies(tf_name)
    occupancies = [sum([1/(1+exp(beta*(ep-mu))) for ep in site_energies])
                   for mu in mus]
    approx_occupancies = [sum([1/(1+exp(beta*(ep-mu))) for ep in site_energies])
                          for mu in approx_mus]
    plt.plot(ks,occupancies,label=r"\mu")
    plt.plot(approx_ks,approx_occupancies,label=r"\hat\mu")
    plt.semilogx()
    plt.xlabel("Copy number")
    plt.ylabel("Regulon Occupancy")
    plt.legend(loc=0)
    plt.savefig(outfile,dpi=400)
    plt.close()
    
