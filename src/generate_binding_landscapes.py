"""This script generates a binding landscape for each TF"""

import sys
sys.path.append("sufficache")
from sufficache import PSSM
sys.path.append("../data/motifs")
from motifs import *
from chem_pot_utils import get_ecoli_genome
from array import array

genome = get_ecoli_genome()

print "Generating binding landscapes"

for tf_name in Escherichia_coli.tfs:
    print "Generating landscape for ",tf_name
    tf = getattr(Escherichia_coli,tf_name)
    pssm = PSSM(tf)
    binding_energies = PSSM.slide_trap(genome)
    arr = array('f')
    arr.extend(binding_energies)
    fname = "../results/binding_landscapes/%s_binding_landscape.dat" % tf_name
    arr.tofile(fname)
    
print "Finished generating binding landscapes"
