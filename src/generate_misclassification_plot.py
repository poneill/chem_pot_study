"""
Generate a histogram of misclassification rates for all tfs with
known copy number
"""

import sys
sys.path.append("data")
from copy_numbers import copy_numbers
from motifs import Escherichia_coli
from chem_pot_utils import *
from math import exp
from params import *
from matplotlib import pyplot as plt

def misclass_rate(tf_name):
    site_energies = get_regulon_binding_energies(tf_name)
    if not tf_name in copy_numbers:
        return None
    copy_number = copy_numbers[tf_name]
    mus_ks = read_mu_k_table(tf_name)
    approx_mus_ks = read_mu_k_table(tf_name,approx=True)
    mu = closest_mu(mus_ks,copy_number)
    approx_mu = closest_mu(approx_mus_ks,copy_number)
    probs = [1/(1+exp(beta*(ep-mu))) for ep in site_energies]
    approx_probs = [1/(1+exp(beta*(ep-approx_mu))) for ep in site_energies]
    misclass_rate = len([p for (p,app_p) in zip(probs,approx_probs)
                         if abs(p - app_p) > 0.5])/float(len(site_energies))
    return misclass_rate

def make_plot(outfile):
    missclass_rates = filter(lambda x:x,map(misclass_rate,Escherichia_coli.tfs))
    plt.hist(missclass_rates)
    plt.xlabel("Misclassification Rate")
    plt.savefig(outfile,dpi=400)
    plt.close()

if __name__ == "__main__":
    outfile = sys.argv[1]
    make_plot(outfile)

    
        
