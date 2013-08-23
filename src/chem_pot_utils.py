from array import array
import csv,sys
sys.path.append("data/motifs")
from motifs import Escherichia_coli
sys.path.append("src/sufficache")
print sys.path
from sufficache import PSSM


def get_ecoli_genome():
    genome_file = "data/NC_000913.fna"
    with open(genome_file) as f:
        raw_genome = "".join([line.strip() for line in f.readlines()[1:]])
    return "".join(b for b in raw_genome if b in "ATGC") # may contain
                                                         # other iupac
                                                         # symbols

def load_array(filename,typ):
    n = 5000000 # this constant must be hard-coded, but in practice we
                # know in advance that the length of the array will
                # always be less than this.
    arr = array(typ)
    with open(filename) as f:
        try:
            arr.fromfile(f,n)
        except:
            pass
    return list(arr)

def read_mu_k_table(tf_name,control=False,approx=False):
    data_file = ("results/mu_k_tables/%s_%s_%s_mu_k_table.csv" %
                 (tf_name,
                  "control" if control else "genome",
                  "approximation" if approx else "exact"))
    with open(data_file) as f:
        r = csv.reader(f)
        table = [(float(mu),float(k)) for (mu,k) in list(r)[1:]] # skip header
    return table

def closest_mu(mus_ks,copy_number):
    """take as mu that which gives us the closest approximation to
    physiological copy number"""
    mu,k = sorted(mus_ks,key=lambda(m,k):(copy_number - k)**2)[0]
    return mu

def get_regulon_binding_energies(tf_name):
    sites = getattr(Escherichia_coli,tf_name)
    pssm = PSSM(sites)
    site_energies = map(pssm.trap,sites)
    return site_energies
