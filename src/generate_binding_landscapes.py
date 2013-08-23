"""This script generates a binding landscape for each TF"""

import sys
sys.path.append("src/sufficache")
from sufficache import PSSM
sys.path.append("data/motifs")
from motifs import *
from chem_pot_utils import get_ecoli_genome
from array import array

if __name__ == "__main__":
    "usage: generate_binding_landscapys.py tf_name [control]"
    tf_name = sys.argv[1]
    control = len(sys.argv) == 3 and sys.argv[2] == "control"
    genome = get_ecoli_genome() if not control else random_site(len(get_ecoli_genome()))
    print "Generating %s landscape for %s " % ("Control" * control,tf_name)
    tf = getattr(Escherichia_coli,tf_name)
    pssm = PSSM(tf)
    binding_energies = pssm.slide_trap(genome)
    arr = array('f')
    arr.extend(binding_energies)
    if not control:
        fname = "results/binding_landscapes/%s_genome_binding_landscape.dat" % tf_name
    else:
        fname = "results/binding_landscapes/%s_control_binding_landscape.dat" % tf_name
    with open(fname,'w') as f:
        arr.tofile(f)

