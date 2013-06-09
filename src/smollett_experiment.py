# The purpose of this script is to compare Mycobacterium tuberculosis
# LexA chip-seq data with the model.

from utils import *
from smollett_data import smollett_data,smollett_sites
from sufficache import PSSM
from chem_pot_experiment import fd_probs,approximate_copy_number_from_mu
from scipy.stats import spearmanr,pearsonr

mus = range(0,50)
def get_mycobacterium_genome():
    with open("NC_000962.fna") as f:
        return "".join([line.strip() for line in f.readlines()[1:]])
        
genome = get_mycobacterium_genome()

#Site1 in smollett_data does not appear in Myco genome.  (We use
#strain H37Rv; they use strain 1424 which is derived from the former.)
#For this reason, we need to revise the genome in order to stitch the site in.

start_coordinate = (3811492 #start position of region listed in Table 1
                    +158) # position of strongest site, relative to
                          # start position

site1 = smollett_sites[1]
revised_genome = subst(genome,site1,start_coordinate)
    
model = PSSM(smollett_sites.values())

traps = model.slide_trap(revised_genome)
exact_copies = [sum(fd_probs(traps,mu,beta)) for mu in verbose_gen(mus)]
z = sum(exp(-beta*ep) for ep in traps)
approx_copies = [approximate_copy_number_from_mu(traps,mu,z)
                     for mu in verbose_gen(mus)]
absolute_ns_energy = -8 #kBT = -5 kca/mol
width = len(site1)
ep_ns = 2*width + absolute_ns_energy #Assume binding energy is -2kbt/match
offset = lambda ep:log(exp(-beta*ep) + exp(-beta*ep_ns))/-beta
ns_traps = map(offset,traps)

coordinates = [smollett_data[i][0] for i in range(1,25+1)]
scores = [smollett_data[i][1] for i in range(1,25+1)]
regions = [genome[start_pos:end_pos+18] for (start_pos,end_pos) in coordinates]
normalized_scores = [score/len(region) for score,region in zip(scores,regions)]
select_traps = [traps[start_pos:end_pos] for (start_pos,end_pos) in coordinates]
ns_select_traps = [ns_traps[start_pos:end_pos] for (start_pos,end_pos) in coordinates]
min_traps = map(min,select_traps)
mean_traps = map(mean,select_traps)

def correlate(mu,method=pearsonr,ns=False):
    selections = select_traps if not ns else ns_select_traps
    predicted_counts = [max([fermi_dirac(trap,mu,beta) for trap in region])/len(region)
                        for region in selections]
    return method(predicted_counts,scores)

def max_probs(mu):
    return [fermi_dirac(mt,mu=mu,beta=beta)
            for mt in min_traps]

def normalized_max_probs(mu):
    return [fermi_dirac(mt,mu=mu,beta=beta)/len(select_trap)
                      for mt,select_trap in zip(min_traps,select_traps)]
    
def sum_probs(mu):
    return [sum(fermi_dirac(t,mu=mu,beta=beta)
                for t in select_trap)
            for select_trap in select_traps]

def normalized_sum_probs(mu):
    return [sum(fermi_dirac(t,mu=mu,beta=beta)
                for t in select_trap)/len(select_trap)
            for select_trap in select_traps]

def plot():
    correlations = [pearsonr(scores,max_probs(mu))[0] for mu in mus]
    correlations2 = [pearsonr(scores,sum_probs(mu))[0]
                     for mu in mus]
    correlations3 = [pearsonr(scores,normalized_sum_probs(mu))[0]
                     for mu in mus]
    correlations4 = [pearsonr(scores,normalized_max_probs(mu))[0]
                     for mu,select_trap in zip(mus,select_traps)]
    exact_indices = [i for i,v in enumerate(exact_copies) if 1 <= v < 10**15]
    approx_indices = [i for i,v in enumerate(approx_copies) if 1 <= v < 10**15]
    plt.plot(rslice(exact_copies,exact_indices),rslice(correlations,exact_indices),
             label="$\mu$")
    plt.plot(rslice(approx_copies,approx_indices),rslice(correlations,approx_indices),
             label="$\hat\mu$")
    plt.xlabel("Copy Number")
    plt.ylabel(r"Pearson $\rho$")
    plt.semilogx()
    plt.legend(loc=0)
    plt.title("Correlation between predicted occupancy and Chip-Seq read density")
    plt.savefig("Smollett_experiment.png",dpi=400)
    

    
