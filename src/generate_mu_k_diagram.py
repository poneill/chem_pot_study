"""This script generates mu-k diagrams."""

import sys,csv,os
from matplotlib import pyplot as plt
sys.path.append("lib/utils")
from utils import transpose
from chem_pot_utils import read_mu_k_table
sys.path.append("data")
from copy_numbers import copy_numbers

def make_plot(mus_ks,
              control_mus_ks,
              approximate_mus_ks,
              control_approximate_mus_ks,
              copy_number,
              outfile):
    plt.plot(*transpose([(k,mu) for (mu,k) in mus_ks]),label=r"$\mu$")
    plt.plot(*transpose([(k,mu) for (mu,k) in control_mus_ks]),label=r"Control $\mu$")
    plt.plot(*transpose([(k,mu) for (mu,k) in approximate_mus_ks]),label=r"$\hat\mu$")
    plt.plot(*transpose([(k,mu) for (mu,k) in control_approximate_mus_ks]),
              label=r"Control $\hat\mu$")
    if copy_number:
        plt.plot([copy_number,copy_number],[0,50],label="Copy number",linestyle="--")
    plt.xlabel("Copy number")
    plt.ylabel("Chemical potential + Const. (kBT)")
    plt.semilogx()
    plt.xlim(1,10**6)
    plt.title("Copy number vs. Chemical Potential")
    plt.legend(loc='upper left')
    plt.savefig(outfile,dpi=400)
    plt.close()
    
if __name__ == "__main__":
    print len(sys.argv)
    print sys.argv
    (mu_k_file,
     control_mu_k_file,
     mu_k_approximate_file,
     mu_k_control_approximate_file,
     outfile) = sys.argv[1:]
    data_files = sys.argv[1:5]
    folder,filename = os.path.split(mu_k_file)
    tf_name = filename[:filename.index("_")]
    copy_number = copy_numbers[tf_name] if tf_name in copy_numbers else None
    mus_ks                     = read_mu_k_table(tf_name,control=False,approx=False)
    control_mus_ks             = read_mu_k_table(tf_name,control=True, approx=False)
    approximate_mus_ks         = read_mu_k_table(tf_name,control=False,approx=True)
    control_approximate_mus_ks = read_mu_k_table(tf_name,control=True, approx=True)
    make_plot(mus_ks,
              control_mus_ks,
              approximate_mus_ks,
              control_approximate_mus_ks,
              copy_number,
              outfile)
