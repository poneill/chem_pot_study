"""This script generates mu-mu-hat diagrams"""

import sys,os
sys.path.append("data/motifs")
from motifs import Escherichia_coli
sys.path.append("data")
from copy_numbers import copy_numbers
from chem_pot_utils import *
sys.path.append("src/sufficache")
from sufficache import PSSM
from math import exp
from params import *
from matplotlib import pyplot as plt

def make_mu_mu_hat_plot(mu,approx_mu,tf_name,outfile):
    site_energies = get_regulon_binding_energies(tf_name)
    probs = [1/(1+exp(beta*(ep-mu))) for ep in site_energies]
    approx_probs = [1/(1+exp(beta*(ep-approx_mu))) for ep in site_energies]
    plt.scatter(probs,approx_probs)
    plt.xlabel(r"Occupancy given $\mu$")
    plt.ylabel(r"Occupancy given $\hat\mu$")
    plt.xlim([-0.05,1.05])
    plt.ylim([-0.05,1.05])
    plt.savefig(outfile,dpi=400)
    plt.close()

if __name__ == "__main__":
    (mu_k_file,
     mu_k_approximate_file,
     outfile) = sys.argv[1:]
    folder,filename = os.path.split(mu_k_file)
    tf_name = filename[:filename.index("_")]
    copy_number = copy_numbers[tf_name] if tf_name in copy_numbers else None
    mus_ks = read_mu_k_table(tf_name,approx=False)
    approximate_mus_ks = read_mu_k_table(tf_name,approx=True)
    if not copy_number:
        print "Found no copy number for",tf_name,"."
        print "Quitting mu_mu_hat plot."
        pass
    else:
        mu = closest_mu(mus_ks,copy_number)
        approx_mu = closest_mu(approximate_mus_ks,copy_number)
        make_mu_mu_hat_plot(mu, approx_mu, tf_name, outfile)
