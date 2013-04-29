#There is a linear relationship between score and log(probability)
#Intercept depends on TF copy number; slope on temp 
import sys
#pypy imports
at_lab = True
sys.path.append("/home/poneill/python_utils")
sys.path.append("/home/poneill/motifs")
sys.path.append("/home/poneill/sufficache")
from utils import *
from motifs import *
#import numpypy as np
#import numpy as np
import sufficache
import random
from math import *
from array import array
from itertools import combinations
from collections import Counter,defaultdict
from copy_dict import copy_dict,ns_copy_dict
from copy_numbers import copy_numbers
if not sys.executable == "/usr/local/bin/pypy":
    """Load these imports only if not using pypy"""
    from scipy.optimize import fmin
    from matplotlib import pyplot as plt
#from chem_pot_data import *
#from mpmath import mpf,exp
#kB = 1.3806503*10**-23
kB  = 0.0019872041 #kcal/mol (!)
temp_C = 37
C2K = 273.15
temp = temp_C + C2K 
beta = 1/(kB*temp) 
R = 1.9858775*10**-3#ideal gas constant
mus = range(-5,50)

def pairs(xs):
    return zip(xs[:-1],xs[1:])

def choose(n,k):
    return factorial(n)/(factorial(k) * factorial(n-k))

def choose_large(n,k):
    return product(range(n-k+1,n+1))/factorial(k)

def product(xs):
    return reduce(lambda x,y:x*y,xs,1)
    
def transpose(xxs):
    return zip(*xxs)

def fast_partials(probs):
    partials = []
    total = 0
    for prob in probs:
        total += prob
        partials.append(total)
    return partials


def delta_g_from_kd(kd,temp):
    return -R*temp*log(kd)

def index(xs, p):
    """Return index of first x satisfying p, or None if none do"""
    winners = filter(lambda (i, x): p(x), zip(range(len(xs)), xs))
    return winners[0][0] if winners else None

def verbose_gen(xs,n=1):
    for (i,x) in enumerate(xs):
        if i % n == 0:
            print i
        yield x
        
def fast_index(xs,p):
    """Return index of first x satisfying p, or None if none do, using
    binary search.  Assumes that xs are sorted according to p, e.g. xs
    is a sorted numeric list, and p is of the form lambda x: x > k"""
    lo = 0
    hi = len(xs)-1
    while lo != hi:
        if hi - lo == 1:
            if p(xs[lo]):
                return lo
            elif p(xs[hi]):
                return hi
            else:
                return None
        guess = int((lo + hi)/2 + random.random())
        if p(xs[guess]):
            hi = guess
        else:
            lo = guess
    if p(xs[guess]):
        return guess
    else:
        return None
    

def normalize(xs):
    total = float(sum(xs))
    return [x/total for x in xs]

def normalize_dict(counts):
    return dict(zip(counts.keys(),normalize(counts.values())))


def motif_ic(motif):
    return sum(columnwise_ic(motif))


def columnwise_ic(seqs,hg=2):
    return [hg - col_h for col_h in columnwise_h(seqs)]


def columnwise_h(seqs):
    return [h(col) for col in zip(*seqs)]

def mean(xs):
    return sum(xs)/float(len(xs))

def gini(xs):
    print "sorting,normalizing"
    xs = normalize(sorted(xs))
    print "constructin ys"
    y = mean(xs)
    ys = [y for x in xs]
    #xs_partials = [sum(xs[:i+1]) for i in range(len(xs))]
    #ys_partials = [sum(ys[:i+1]) for i in range(len(ys))]
    print "partial xs"
    xs_partials = (fast_partials(xs))
    print "partial ys"
    ys_partials = (fast_partials(ys))
    print "summing partial xs"
    b = sum(xs_partials)
    print "summing partial ys"
    a = sum(ys_partials) - b
    return a/(a + b)

def motif_gini_coefficient(motif):
    ics = sorted(columnwise_ic(motif))
    return gini(ics)

def freq(b,seq):
    return len([c for c in seq if c ==b])/float(len(seq))

def safelog2(x):
    epsilon = 10**-10
    return log(x+epsilon,2)

def h(seq):
    return -1*sum([freq(b,seq)*safelog2(freq(b,seq)) for b in set(seq)]) 

def sample(xs,probs):
    partials = fast_partials(probs)
    r = random.random()
    i = fast_index(partials,lambda x: x > r)
    return xs[i]

def sample_probs(probs):
    partials = fast_partials(probs)
    r = random.random()
    i = fast_index(partials,lambda x: x > r)
    return i

def fast_sample_probs(probs):
    r = random.random()
    total = 0
    for i in range(len(probs)):
        total += probs[i]
        if total > r:
            return i

        
def log_sample(xs,log_props):
    """sample from xs, given the log_propensities"""
    z = sum(log_props)
    pre_partials = fast_partials([lp for lp in log_props])
    def f(x):
        x_0 = pre_partials[0]
        x_n = pre_partials[-1]
        m = 1/(x_n - x_0)
        return m*(x - x_0)
    partials = map(f,pre_partials)
    r = random.random()
    log_r = log(r)
    i = fast_index(partials,lambda x: x > log_r)
    return xs[i]
    
    
def importance_sample(xs,probs):
    """Sample z from Q, u from U(0,1).  If u < P(z)/kQ(z), accept.  k
    is the scaling factor for Q to ensure that Q is an envelope for P."""
    n = len(probs)
    proposed = random.randrange(n)

    
def tree_sample(xs,probs,r=random.random()):
    """Inverse CDF sample xs using probs with a binary search tree """
    if len(xs) == 1:
        return xs[0]
    else:
        pass

def mutate_site(site):
    pos = random.randrange(len(site))
    return site[:pos] + random.choice("ACGT") + site[pos+1:]

def chem_pot_experiment(n=10,g=100000,iterations=10):
    """In this experiment, we simulate the behavior of multiple
    transcription factors binding in a genome.  Assuming Fermi-Dirac
    statistics, the probability of a site being bound by any of n TFs
    is given by: P(bound|e_i) = 1/(e^[(e_i-mu)*n/(kBT)] + 1).  We wish
    to test the goodness of fit of the approximation mu = F_0/(kBT) +
    log(n) for the chemical potential, where F_0 = -kBT*log(Z_0) is
    the free energy of a random genome, Z_0 is the partition function
    of the same."""
    #motif = [sufficache.utils.random_site(tf_width) for i in range(10)]
    motif = [mutate_site("a"*tf_width) for i in range(20)]
    tf = sufficache.PSSM(motif)
    tf_width = len(tf.motif[0])
    genome = sufficache.utils.random_site(g + tf_width)
    scores = [tf.score(genome[i:i+tf_width]) for i in range(100000)]
    scale_function = biochem_scale_output(tf)
    energies = [scale_function(score) for score in scores]
    props = [exp(-beta*energy) for energy in energies]
    #log_props = [(-beta*energy) for energy in energies]
    bound_positions = []
    open_positions = range(g)
    bound_counts = [0] * g
    for iteration in xrange(iterations):
        print iteration
        for i in range(n):
            print "tf:",i
            open_positions = [site for site in range(g)
                              if not site in bound_positions]            
            open_props = [props[j] for j in open_positions]
            #open_log_props = [log_props[j] for j in open_positions]
            z = sum(open_props)
            #log_z = sum(log_open_props)
            open_probs = [op/z for op in open_props]
            bound_site = sample(open_positions,open_probs)
            #log_bound_site = log_sample(open_positions,open_log_props)
            bound_counts[bound_site] += 1
            bound_positions.append(bound_site)
        print bound_positions
        bound_positions = []
    return bound_counts

def compute_gs(energies):
    """compute the multiplicity of the energy levels: for each site in
    scores, find the number of sites sharing that score.  See Sec 2 of
    chemical_potential.pdf"""
    counts = Counter(energies)
    return [counts[energy] for energy in energies]

def compute_cs(energies,beta):
    return [exp(-beta * energy) for energy in energies]

def root_find_b(energies,num_tfs,beta,mu_tolerance):
    # n = num_tfs
    # print "computing gs"
    # gs = compute_gs(energies)
    # print "computing cs"
    # cs = compute_cs(energies,beta) #even though this takes 11 sec, worth it
    # def f(b):
    #     return sum(gs[i]/(b*cs[i] + 1)
    #                for i in xrange(len(gs))) - n
    print "computing counts"
    counts = Counter(energies)
    print "computing cs"
    cs = {count:exp(-beta*count) for count in counts}
    f = lambda(b):safe_f(b,counts,cs,num_tfs)
    min_mu = -50 #magic numbers come from biochemical savvy
    max_mu = 50
    print "computing b"
    b = bisect(f,min_mu,max_mu,beta=beta, mu_tolerance=mu_tolerance)
    return b

def safe_f(b,counts,cs,n):
    """Compute f by summing over energy equivalence classes, not sites"""
    return sum(counts[count]**2/(b*cs[count] + 1) for count in counts) - n
    
def verbose_sum(xs):
    total = 0
    for x in xs:
        total += x
        print total
    return total

def compute_mu(b,beta):
    """compute mu from: b = exp(-beta*mu)"""
    return log(b)/-beta

def compute_b(mu,beta):
    """compute b from: b = exp(-beta*mu)"""
    return exp(-beta*mu)

def fermi_dirac(e,mu,beta=beta):
    return 1/(exp((e-mu)*beta) + 1)

def fd(e_i,beta,mu):
    return 1/(exp(beta*(e_i - mu)) + 1)

def bisect(f,mu_lo,mu_hi,beta=beta,mu_tolerance=1e-2,last_mu_diff=None,verbose=False):
    """find x such that f(x) == 0, x in [lo,high]"""
    current_mu_diff = mu_hi - mu_lo
    mu_guess = (mu_lo + mu_hi)/2
    b_guess = compute_b(mu_guess,beta)
    y_guess = f(b_guess)
    if abs(mu_hi - mu_lo) < mu_tolerance:
        return b_guess
    elif current_mu_diff == last_mu_diff:
        if verbose:
            print "warning: could not compute mu to desired tolerance"
            print "mu within tolerance:",current_mu_diff
        return b_guess
    else:
        if verbose:
            print "mu_lo",mu_lo,"mu_guess",mu_guess, "mu_hi",mu_hi,"diff",mu_hi-mu_lo
        """zero lies to right of guess"""
        if y_guess < 0:
            if verbose:
                print "y_guess:",y_guess, "raising left endpoint"
            new_mu_lo = mu_guess
            new_mu_hi = mu_hi
        else:
            if verbose:
                print "y_guess:",y_guess, "lowering right endpoint"
            new_mu_lo = mu_lo
            new_mu_hi = mu_guess
        return bisect(f,new_mu_lo,new_mu_hi,
                      mu_tolerance=mu_tolerance,beta=beta,
                      last_mu_diff=current_mu_diff)
        
def secant_method(f,mu_lo,mu_hi,y_lo,y_hi,
                  beta=beta,mu_tolerance=1e-2):
    """find x such that f(x) == 0, x in [lo,high]. DEPRECATED --
    bisection faster in practice"""
    
    m = ((y_hi) - (y_lo))/(mu_hi - mu_lo)
    mu_guess = -(y_lo)/m + mu_lo
    b_guess = compute_b(mu_guess,beta)
    if abs(mu_hi - mu_lo) < mu_tolerance:
        return b_guess
    else:
        print "range:",lo,hi
        y_guess = f(b_guess)
        if y_lo * y_guess > 0:
            """zero lies to right of guess"""
            print "y_guess:",y_guess, "raising left endpoint"
            print "mus:",mu_lo,mu_guess,mu_hi
            return secant_method(f,mu_guess,mu_hi,y_guess,y_hi,beta,mu_tolerance)
        else:
            print "y_guess:",y_guess, "lowering right endpoint"
            print "mus:",mu_lo,mu_guess,mu_hi
            return secant_method(f,mu_lo,mu_guess,y_lo,y_guess,beta,mu_tolerance)
    

def biochem_scale_output(tf,y_min=-.05,energy_per_nt = -2):
    """Return a function that scales the output of a tf (qua PSSM
    object) so that its minimum output is -.05, and its maximum is
    -17.  I.e. return minimum and maximum energy levels consistent
    with empirical data."""
    x_min = sum([min(column) for column in tf.columns])
    x_max = sum([max(column) for column in tf.columns])
    #y_max = -17
    y_max = energy_per_nt * len(tf.motif[0])
    m = (y_max - y_min)/(x_max - x_min)
    b = y_max - m * x_max
    def f(x):
        return m*x + b
    return f
    
def fermi_dirac_exp():
    #motif = [mutate_site("a"*tf_width) for i in range(20)]
    print "setting up"
    motif = get_crp_motif()
    tf = sufficache.PSSM(motif)
    tf_width = len(tf.motif[0])
    scale_function = biochem_scale_output(tf)
    #genome = sufficache.utils.random_site(g + tf_width)
    genome = get_ecoli_genome()
    num_tfs = 1
    #temp = 300
    print "scoring"
    scores = [tf.score(genome[i:i+tf_width])
              for i in verbose_gen(range(len(genome)-tf_width + 1),1000)]
    print "energies"
    energies = [scale_function(score) for score in scores]
    del(scores)
    print "computing b"
    mu_tolerance = 1e-2
    b = root_find_b(energies,num_tfs,beta,mu_tolerance=mu_tolerance)
    mu = compute_mu(b,beta)
    fd_probs = [fermi_dirac(e,mu,beta) for e in energies]
    return fd_probs

def het_index_exp():
    #motif = [mutate_site("a"*tf_width) for i in range(20)]
    print "setting up"
    lexa = get_lexa_motif()
    tf = PSSM(lexa)
    tf_width = len(tf.motif[0])
    #scale_function = biochem_scale_output(tf)
    #genome = sufficache.utils.random_site(g + tf_width)
    genome = get_ecoli_genome()
    num_tfs = 1
    #temp = 300
    print "scoring"
    het_indices = [tf.het_index(genome[i:i+tf_width])
                   for i in verbose_gen(range(len(genome)-tf_width + 1),1000)]
    print "energies"
    k_0 = -1 #see Berg 1988
    energies = [k_0 *exp(-hi) for hi in het_indices]
    print "computing b"
    mu_tolerance = 1e-10
    b = root_find_b(energies,num_tfs,beta,mu_tolerance=mu_tolerance)
    mu = compute_mu(b,beta)
    fd_probs = [fermi_dirac(e,mu,beta) for e in energies]
    return fd_probs


def mb_probs_from_energies(energies,beta=beta):
    props = [exp(-beta*e) for e in energies]
    z = sum(props)
    return [p/z for p in props]

def fd_probs_from_energies(energies,num_tfs,beta=beta,mu_tolerance = 1e-10):
    mu = compute_mu_from_energies(energies,num_tfs,mu_tolerance,method=fd_bisect)
    print "mu:",mu
    fd_probs = [fermi_dirac(e,mu,beta) for e in energies]
    return fd_probs

def fd_probs_from_counts(counts,num_tfs,beta=beta):
    mu_tolerance = 1e-10
    mu = compute_mu_from_energies(counts,num_tfs,mu_tolerance,method=fd_bisect)
    print "fd_probs"
    fd_probs = concat([[fermi_dirac(e,mu,beta)] * count
                for (e,count) in counts.iteritems()])
    return fd_probs

def ghost_probs(energies,num_tfs,beta=beta):
    mb_probs = mb_probs_from_energies(energies,beta)
    return [1-(1-prob)**num_tfs for prob in mb_probs]

def fd_probs(energies,mu,beta=beta):
    if type(energies) is list:
        return [fermi_dirac(e,mu,beta) for e in energies]
    elif isinstance(energies,dict):
        return {e:fermi_dirac(e,mu,beta) for e in energies}
    else:
        assert(False)

def get_copy_number(energies,mu,beta=beta):
    k = 0
    for energy in energies:
        k += fermi_dirac(energy,mu,beta)
    return k

def fd_bisect(lo,hi,n,relative_error,total_prob):
    guess = (lo + hi)/2.0
    total = total_prob(guess)
    if abs(total - 1) < relative_error: #nb: total is already
                                        #normalized by copy number
        return guess
    else:
        print lo,hi,guess,total
        if total > 1:
            """chemical potential should be lower"""
            new_lo = lo
            new_hi = guess
        else:
            """chemical potential should be higher"""
            new_lo = guess
            new_hi = hi
        return fd_bisect(new_lo,new_hi,n,relative_error,total_prob)


def fd_secant(lo,hi,n,relative_error,total_prob,y_lo=None,y_hi=None):
    #NB: do not rely on this code!  solves for n = 1!
    f = lambda mu: log(total_prob(mu))
    if not y_lo:
        y_lo = f(lo)
    if not y_hi:
        y_hi = f(hi)
    m = ((y_hi) - (y_lo))/(hi - lo)
    guess = -(y_lo)/m + lo
    y_guess = (f((guess)))
    if abs(y_guess) < relative_error:
        return guess
    else:
        print lo,y_lo,hi,y_hi,guess,y_guess
        if y_guess > 0:
            """chemical potential should be lower"""
            new_lo = lo
            new_y_lo = y_lo
            new_hi = guess
            new_y_hi = y_guess
        else:
            """chemical potential should be higher"""
            new_lo = guess
            new_y_lo = y_guess
            new_hi = hi
            new_y_hi = y_hi
        return fd_secant(new_lo,new_hi,n,relative_error,total_prob,new_y_lo,new_y_hi)


def compute_fermi_dirac_probs(motif,genome,num_tfs,beta=beta):
    tf = sufficache.PSSM(motif)
    tf_width = len(tf.motif[0])
    scale_function = biochem_scale_output(tf)
    print "scores"
    scores = [tf.score(genome[i:i+tf_width])
              for i in verbose_gen(range(len(genome)-tf_width + 1),100000)]
    print "energies"
    energies = [scale_function(score) for score in scores]
    print "computing b"
    mu_tolerance = 1e-10
    b = root_find_b(energies,num_tfs,beta,mu_tolerance=mu_tolerance)
    mu = compute_mu(b,beta)
    print "fd_probs"
    fd_probs = [fermi_dirac(e,mu,beta) for e in energies]
    return fd_probs
    
def x(beta,mu,e_i):
    return exp(beta*(e_i-mu))

def get_ecoli_genome(apec=False):
    if apec:
        lab_file = "/home/poneill/ecoli/NC_008563.fna"
    else:
        lab_file = "/home/poneill/ecoli/NC_000913.fna"
    home_file = "/home/pat/Dropbox/entropy/NC_008563.fna"
    with open(lab_file if at_lab else home_file) as f:
        genome = "".join([line.strip() for line in f.readlines()[1:]])
    return "".join(g for g in genome if g in "ATGC") # contains other iupac symbols

def get_crp_motif():
    lab_file = "/home/poneill/euksites/crp.csv"
    home_file = "/home/pat/Dropbox/entropy/crp.txt"
    with open(lab_file if at_lab else home_file) as f:
        lines = [line.strip() for line in f.readlines()[1:]]
    return [line.replace(",","") for line in lines]


def get_lexa_e_coli_motif():
    home_file = "/home/pat/Dropbox/lexa.csv"
    lab_file = "/home/poneill/lexa/LexA_E_coli_120.csv"
    with open(lab_file if at_lab else home_file) as f:
        lines = [line.strip() for line in f.readlines()[1:]]
    return [line.replace(",","")[50:70] for line in lines]

def get_lexa_b_subtilis_motif():
    home_file = "/home/pat/Dropbox/lexa.csv"
    lab_file = "/home/poneill/lexa/LexA_B_subtilis_116.csv"
    with open(lab_file if at_lab else home_file) as f:
        lines = [line.strip() for line in f.readlines()[1:]]
    return [line.replace(",","")[51:65] for line in lines]


def entropy(probs):
    return -sum(p*log(p,2) for p in probs)

def entropy_from_counts(counts,num_tfs):
    probs = fd_probs_from_counts(counts,num_tfs)
    return entropy(normalizeprobs)
    
def counts_entropy(counts):
    return -sum(p*log(p,2) for p in counts.values())

def rFreq(probs):
    return log(len(probs),2) - entropy(probs)

def numerical_exp():
    print "sorting"
    top = sorted(energies,reverse=True)[:100000]
    print "finished sorting"
    for i in [1,2,5,10,20,50,100,500,1000,10000]:
 	f = lambda(b):safe_f(b,counts,cs,i)
 	b = bisect(f,-50,-.05,beta=beta, mu_tolerance=1e-2)
 	mu = compute_mu(b,beta)
 	print sum([fermi_dirac(e,mu,beta) for e in top])

print("loaded chem_pot_experiment")

def countify(energies_or_counts):
    if type(energies_or_counts) is list:
        counts = Counter(energies_or_counts)
    else:
        counts = energies_or_counts
    return counts

def uncountify(counts):
    return sum([[k] * counts[k] for k in verbose_gen(counts,1e5)],[])
    
def compute_mu_from_energies(energies_or_counts,num_tfs,relative_error,
                             method=fd_bisect,lb=-200,ub=200):
    counts = countify(energies_or_counts)
    def total_prob(mu):
        """sum probability over all sites; should sum to num_tfs;divide by n"""
        return sum(fermi_dirac(energy,mu,beta) * counts[energy]
                   for energy in counts)/float(num_tfs)
    mu = method(lb,ub,num_tfs,relative_error,total_prob)
    return mu
    
def fermi_energy(energies_or_counts,num_tfs):
    if isinstance(energies_or_counts,dict):
        energies = expand_counts(energies_or_counts)
    else:
        energies = energies_or_counts
    return sorted(energies)[num_tfs-1]

def unnormalized_entropy_exp(n):
    for i in range(n):
        print i
        energies = ([random.random() for i in range(100)])
        #plt.plot([entropy(fermi_dirac_probs(energies,j)) for j in range(1,100+1)])
        plt.plot([entropy(map(lambda x:x*j,normalize(energies))) for j in range(1,100+1)])
    plt.show()

def log2(x):
    return log(x,2)

def log10(x):
    return log(x,10)

comp = {"a":"t",
        "t":"a",
        "g":"c",
        "c":"g",
        "A":"T",
        "T":"A",
        "G":"C",
        "C":"G"}

def wc(dna):
    return "".join([comp[b] for b in dna])[::-1]

def auto_reg_exp(n=1,basal_rate=1,binding_factor=10000,iterates=10,counts=None):
    crp = sufficache.PSSM(get_crp_motif())
    crp_binding_site = min([(site,crp.trap(site))
                            for site in crp.motif],key=lambda(x,y):y)[0]
    e = crp.trap(crp_binding_site)
    if counts is None:
        counts = Counter(crp.slide_trap(get_ecoli_genome()))
        if counts[e] == 0: #make sure the binding site itself is in the genome!
            counts[e] = 1
    def copy_num_from_prob(p):
        return p * binding_factor + basal_rate #don't interpret; just
                                             #assume it's some affine
                                             #function
    ps = []
    ns = []
    mus = []
    for i in range(iterates):
        mu = compute_mu_from_energies(counts,n,1e-10)
        p = fermi_dirac(e,mu)
        n = copy_num_from_prob(p)
        ps.append(p)
        ns.append(n)
        mus.append(mu)
        print p,n,mu
    return (ps,ns,mus)
    
def main():
    crp = sufficache.PSSM([site.upper() for site in get_crp_motif()])
    genome = get_ecoli_genome().upper()
    esa = sufficache.ESA(genome)
    # crp_scale_function = biochem_scale_output(crp)
    # crp_energies = sorted(map(crp_scale_function,crp.slide_score(get_ecoli_genome())))
    # crp_traps = crp.slide_trap(get_ecoli_genome())

def diffs(xs):
    return map(lambda (x,y): x-y,zip(xs+[None],[None]+xs)[1:len(xs)])

def genomic_nmers(genome,n):
    g = len(genome)
    for i in xrange(g - n + 1):
        yield genome[i:i+n]

def sample_marginal(energies,bound_states):
    """Given a list of energies and a list of bound states, sample
    from the marginal M-B distribution on unbound states"""
    props = [exp(-beta*e) if not i in bound_states else 0
             for (i,e) in enumerate(energies)]
    z = sum(props)
    return fast_sample_probs([p/z for p in props])

def gibbs_double_stranded(forwards,backs,n,iterations,beta=beta):
    """given a list of scores corresponding to the 5'-3' and 3'-5'
    scores for a sequence, simulate the binding behavior of n tfs
    given that binding on the forward strand excludes binding on the
    backward strand."""
    counts = defaultdict(int)
    sampled = defaultdict(int)
    fb_sequence = forwards+backs
    bound_states = []
    length = len(forwards)
    #seed tfs
    def both_strands(i):
        if i >= length:
            return [i - length,i]
        else:
            return [i,i + length]
    for i in range(n):
        bound_state = sample_marginal(fb_sequence,bound_states)
        bound_states.extend(both_strands(bound_state))
        print bound_states
    for i in verbose_gen(xrange(iterations),10000):
        remove_indices = both_strands(random.choice(bound_states))
        #print "remove indices:", remove_indices
        #print "bound_states:", bound_states
        bound_states.remove(remove_indices[0])
        bound_states.remove(remove_indices[1])
        sampled_index = sample_marginal(fb_sequence,bound_states)
        sampled[sampled_index] += 1
        bound_states.extend(both_strands(sampled_index))
        for bs in bound_states:
            counts[bs] += 1
    return [counts[i]/(float(iterations)) if i in counts else 0
            for i in range(len(fb_sequence))], sampled
    
    
def sample_marginal_old(energies,bound_states):
    length = len(energies)
    unbound_states = [i for i in range(length) if not i in bound_states]
    unbound_energies = [e for (i,e) in enumerate(energies)
                        if i in unbound_states]
    qs = mb_probs_from_energies(unbound_energies)
    return unbound_states[sample(range(len(qs)),qs)]

def get_factoradic(factoradic_rep,xs):
    ys = xs[:]
    zs = []
    for i in factoradic_rep:
        x = ys[i]
        ys.remove(x)
        zs.append(x)
    return zs

def increment_fac_rep(fac_rep):
    rf = fac_rep[::-1]
    carry = 1
    for i in range(len(rf)):
        rf[i] += carry
        if rf[i] > i:
            rf[i] = 0
            carry = 1
        else:
            carry = 0
    return rf[::-1]        

def config_energy(config,energies):
    #print config
    return sum([energies[i] for i in config]) #this (explicit list )beats a
                                                  #generator
                                                  #expression by a factor of two!
def subset_sum_exp(energies,k,r,s):
    """Is there a k-combination of energies that sums to between r and s?"""
    energies = sorted(energies)
    memo = {}
    def q(i,r,s):
        if (i,r,s) in memo:
            return memo[(i,r,s)]
        elif i == 0:
            memo[(i,r,s)] = (r <= energies[i] <= s)
        else:
            cur = energies[i]
            memo[(i,r,s)] = (q(i-1,r,s) or
                             q(i-1,r - cur,s - cur) or
                             r <= cur <= s)
        print (i,r,s),"=>",memo[(i,r,s)]
        return memo[(i,r,s)]
    
    def q2(es,r,s):
        if (es,r,s) in memo:
            return memo[(i,r,s)]
        elif len(es) == 0:
            return False
        else:
            for cur in energies:
                memo[(i,r,s)] = (q2(i-1,r,s) or
                                 q2(i-1,r - cur,s - cur) or
                                 r <= cur <= s)
        print (i,r,s),"=>",memo[(i,r,s)]
        return memo[(i,r,s)]

def prefix_sum(config,i):
    return sum(config[:i])

def ith_remainder(cur_config,prop_config,i):
    return sum(cur_config) - prefix_sum(prop_config,i)

def index_by(xs,indices):
    return [xs[i] for i in indices]

def too_high_or_low_test(cur_energy,mng,prop_config,energies,k,n):
    print "starting thotl test w/",prop_config
    for i in range(1,k+1):
        prefix = prop_config[:i]
        prefix_energy = config_energy(prefix,energies)
        print "prefix:",prefix
        print "prefix energy:",prefix_energy
        if prefix_energy > mng:
            print "rejected:",prop_config,"as too high at prefix",i
            return ("high",i)
        r = cur_energy - prefix_energy
        #max_possible_remaining_energy = (k - i) * energies[prop_config[i]]
        print "max_possible indices:",n-(k-1),n
        max_possible_remaining_energy = sum(energies[n-(k - i):n])
        print "r:",r
        print "max possible remaining energy:",max_possible_remaining_energy
        if r > max_possible_remaining_energy:
            print "rejected:",prop_config,"as too low at prefix",i
            return ("low",i)
    return ("passed",0)

def update_indices_deprecated(indices,k,n,i=None):
    inds = list(indices)
    if not i is None:
        inds[i-1] += 1
        for offset,j in enumerate(range(i-1,k)):
            inds[j] = inds[i-1] + offset
    else:
        i = k - 1
        inds[i] += 1
        carry = False
        while inds[i] > n - (k - i + -1):
            i -= 1
            inds[i] += 1
            carry = True
        if carry:
            print i,k
            for offset,j in enumerate(range(i,k)):
                inds[j] = inds[i] + offset
    return tuple(inds)

def update_indices(indices,n,i=None):
    """Return the next k-combination of n, optionally the next one not
    sharing an i-prefix with the current combination"""
    print "updating indices:",indices,n,i
    k = len(indices)
    if i == None:
        i = k
    inds = list(indices)
    i -= 1
    inds[i] += 1
    carry = False
    while inds[i] > n - (k - i):
        i -= 1
        inds[i] += 1
        carry = True
    if carry:
        print i,k
        for offset,j in enumerate(range(i,k)):
            print "ofsetting",offset,j
            inds[j] = inds[i] + offset
    return tuple(inds)

def update_indices2(indices,n,i=None):
    print "entering update_indices2:",indices,n,i
    k = len(indices)
    if i == None:
        i = k
    prefix = indices[:i]
    updated_prefix = update_indices(prefix,n-(k-i))
    last = updated_prefix[-1]
    return updated_prefix + tuple(range(last + 1,last + 1 + k - i))

def update_indices_test():
    n = 10
    k = 3
    init_config = tuple(range(k))
    configs = combinations(range(n),k)
    configs_by_update = iterate_list(lambda c:update_indices(c,n),
                                     init_config,choose(n,k)-1)

def lexicographic_cmp(indices1,indices2):
    comps = filter(lambda x: x!=0,zipWith(cmp,indices1,indices2))
    return comps[0] if comps else 0
    
def next_config_by_score(energies,cur_config):
    print "entering next config by score:",cur_config
    n = len(energies)
    k = len(cur_config)
    prop_config = range(k)
    ultimate_config = tuple(range(n-k,n))
    if cur_config == ultimate_config:
        return ultimate_config
    cur_energy = config_energy(cur_config,energies)
    neighborhood = neighbors(cur_config,n)
    best_energy = min(map(lambda c:config_energy(c,energies),neighborhood))
    best_configs = filter(lambda c: (config_energy(c,energies) == best_energy and
                                     lexicographic_cmp(cur_config,c) == -1),
                          neighborhood)
    best_config = sorted(best_configs,cmp = lexicographic_cmp)[0]
    print "best config, best energy:",best_config,best_energy
    #test to see if energies too high
    while prop_config != ultimate_config:
        print "cur_config:",cur_config,cur_energy
        print "prop_config:",prop_config,config_energy(prop_config,energies)
        print "best_config:",best_config,best_energy
        (status,i) = too_high_or_low_test(cur_energy,best_energy,
                                          prop_config,energies,k,n)
        if status == "high":
            print "iffing"
            print "config before updating:",prop_config
            if i == 0:
                break
            prop_config = update_indices2(prop_config,n,i)
            print "config after updating:",prop_config
        elif status == "low":
            print "eliffing"
            print "config before updating:",prop_config
            if i == 0:
                break
            prop_config = update_indices2(prop_config,n,i)
            print "config after updating:",prop_config
        else:
            print "elsing"
            prop_energy = config_energy(prop_config,energies)
            print "comparing",cur_energy,prop_energy,best_energy
            print "comparing configs:",cur_config,prop_config
            if (cur_energy < prop_energy < best_energy or
                (cur_energy <= prop_energy <= best_energy and
                 lexicographic_cmp(cur_config,prop_config) == -1 or
                 lexicographic_cmp(prop_config,best_config) == 1)):
                print "updating best energy"
                best_energy = prop_energy
                best_config = prop_config
                if cur_energy == prop_energy:
                    print "returning best config:",best_config
                    return best_config
            print "config before updating:",prop_config
            prop_config = update_indices2(prop_config,n)
            print "config after updating:",prop_config
    print "returning best config:",best_config
    return best_config

def next_config_by_score_test():
    n = 10
    energies = range(n)
    k = 3
    configs = list(combinations(range(n),k))
    configs_by_sort = sorted(configs,
                             key = lambda c:config_energy(c,energies))
    configs_by_score = iterate_list(lambda c: next_config_by_score(energies,c),
                                  tuple(range(k)),choose(n,k)-1)
    
    
def lattice_exp():
    """enumerate k-combinations of energies by sum,using the partial
    order induced by strict domination of one combination's energy
    levels by another"""
    def dominates(x1,x2):
        """x1 dominates x2: x1's energies are all less then or equal
        to x2s'.  if x1 dominates x2 return 1.  If x2 dominates x1
        return -1.  If incomparable, return 0"""
        comps = (zipWith(lambda c1,c2: c1<=c2),x1,x2)
        if all([c == True for c in comps]):
            return 1
        elif all([c == False for c in comps]):
            return -1
        else:
            return 0
        
    domination_graph = {}
    energies = range(5)
    n = len(energies)
    k = 2
    configs = lambda:(combinations(range(n),k))
    cur_config = configs().next()
    for i in xrange(choose_large(n,k)):
        yield cur_config
        mng = max_neighbor_gain(cur_config,energies)
        

def is_valid(config,n):
    return (len(set(config)) == len(config) and
            tuple(sorted(config)) == config
            and max(config) < n)
    
def neighbors(config,n):
    k = len(config)
    neighbor_vectors = [(0,)*(i) + (1,) + (0,)*(k-i-1) for i in range(k)]
    candidates = map(lambda v:tuple(zipWith(lambda x,y:x+y,config,v)),
                     neighbor_vectors)
    return filter(lambda c: is_valid(c,n),candidates)

def max_neighbor_gain(config,energies):
    neighbor_energies = map(lambda c: config_energy(c,energies),neighbors)
    return min(neighbor_energies)
    
def enumerate_configs_by_energy(energies,k):
    neighbor_vectors = [(0,)*(i) + (1,) + (0,)*(k-i-1) for i in range(k)]
    is_valid = lambda config:(len(set(config)) == k and
                              tuple(sorted(config)) == config
                              and max(config) < n)
    energies = sorted(energies)
    n = len(energies)
    configs = (combinations(range(len(energies)),k))
    cur_config = tuple(range(k))
    cur_energy = config_energy(cur_config,energies)
    neighbors = filter(lambda config: is_valid(config),
                       map(lambda v:tuple(zipWith(lambda x,y:x+y,cur_config,v)),
                         neighbor_vectors))
    neighbor_energies = map(lambda c: config_energy(c,energies),neighbors)
    max_gain = min(neighbor_energies)
def config_mass_exp():
    energies = range(10)
    n = 2
    configs = list(combinations(range(len(energies)),n))
    for config in configs:
        print config,config_energy(config,energies)
    print
    sorted_configs = sorted(configs,key=lambda c:config_energy(c,energies))
    for sorted_config in sorted_configs:
        print sorted_config,config_energy(sorted_config,energies)
    
def gibbs_matrix_analytic(energies,n):
    configs = list(combinations(range(len(energies)),n))
    def energy(config):
        return sum(energies[i] for i in config)
    def relates(config1,config2):
        return len(set(config1).intersection(set(config2))) >= n - 1
    matrix = []
    for config1 in configs:
        row_probs = [exp(-beta*energy(config2)) if relates(config1,config2) else 0
                     for config2 in configs]
        z = sum(row_probs)
        matrix.append([p/z for p in row_probs])
    return matrix

def gibbs_matrix_analytic2(energies,n):
    configs = list(combinations(range(len(energies)),n))
    def energy(config):
        return sum(energies[i] for i in config)
    def relates(config1,config2,i):
        return (all(c1 in config2 for c1 in config1 if c1 != i) or
                all(c2 in config1 for c2 in config2 if c2 != i))
    def relating_set(config,i):
        return filter(lambda c:relates(c,config,i),configs)
    matrix = []
    for config1 in verbose_gen(configs):
        row_probs = map(mean,zip(*[[exp(-beta*energy(config2)) * relates(config1,config2,i) / sum([exp(-beta*energy(c)) for c in relating_set(config1,i)])
                                    for config2 in configs]
                                   for i in config1]))
        matrix.append(row_probs)
    return matrix


def gibbs_matrix_analytic3(energies,n):
    configs = list(combinations(range(len(energies)),n))
    def energy(config):
        return sum(energies[i] for i in config)
    def relates(config1,config2,i):
        return (all(c1 in config2 for c1 in config1 if c1 != i) or
                all(c2 in config1 for c2 in config2 if c2 != i))
    def relating_set(config,i):
        return filter(lambda c:relates(c,config,i),configs)
    matrix = []
    for config1 in verbose_gen(configs):
        rs = {i:relating_set(config1,i) for i in config1}
        row_probs = map(mean,zip(*[[exp(-beta*energy(config2)) * relates(config1,config2,i) / sum([exp(-beta*energy(c)) for c in rs[i]])
                                    for config2 in configs]
                                   for i in config1]))
        matrix.append(row_probs)
    return matrix


def gibbs_analytic(energies,n,beta=beta):
    print "computing matrix"
    matrix = gibbs_matrix_analytic3(energies,n)
    configs = list(combinations(range(len(energies)),n))
    init_config_vector = [[1] + ([0]*(len(configs)-1))]
    print "computing config vector"
    config_vector = converge2(lambda state: matrix_mult(state,matrix),
                             init_config_vector,verbose=True)[0]
    print "projecting etas"
    etas = project_etas(configs,config_vector)
    return etas

def gibbs_more_analytic(energies,k,beta=beta):
    configs = lambda:(combinations(range(len(energies)),k))
    print "num configs:",choose(len(energies),k)
    z = sum(exp(-beta*config_energy(config,energies)) for config in configs())
    print "computing kappas"
    kappas_final = [exp(-beta*config_energy(config,energies))/z
                    for config in configs()]
    print "projecting etas"
    return project_etas(list(configs()),kappas_final)

def my_exp(x):
    return 1 + x + x**2/2 + x**6/6

def gibbs_more_analytic_opt(energies,k,beta=beta):
    positions = range(len(energies))
    configs = lambda:(combinations(positions,k))
    #print "num configs:",choose(len(energies),k)
    #z = sum(exp(beta*config_energy(config,energies)) for config in configs())
    #print "computing kappas"
    #print "projecting etas"
    etas = [0 for i in positions]
    z = 0
    for config in configs():
        energy = exp(-beta*config_energy(config,energies))
        for i in config:
            etas[i] += energy
        z += energy
    return map(lambda eta:eta/z,etas)

def cutoff_gen(xs,n):
    for i,x in zip(range(n),xs):
        yield x
    
def gibbs_more_analytic_opt_cutoff(energies,k,cutoff,beta=beta):
    positions = range(len(energies))
    configs = lambda:cutoff_gen(combinations(positions,k),cutoff)
    #print "num configs:",choose(len(energies),k)
    #z = sum(exp(-beta*config_energy(config,energies)) for config in configs())
    #print "computing kappas"
    #print "projecting etas"
    etas = [0 for i in positions]
    z = 0
    for config in configs():
        energy = exp(-beta*config_energy(config,energies))
        for i in config:
            etas[i] += energy
        z += energy
    return map(lambda eta:eta/z,etas)

def gibbs_most_analytic_opt_cutoff(energies,k,cutoff_percentile=0.99,beta=beta):
    positions = range(len(energies))
    configs = lambda:combinations(positions,k)
    num_configs = choose_large(len(energies),k)
    remaining_configs = num_configs
    #print "num configs:",choose(len(energies),k)
    #z = sum(exp(-beta*config_energy(config,energies)) for config in configs())
    #print "computing kappas"
    #print "projecting etas"
    etas = [0 for i in positions]
    z = 0
    last_mass = 0
    for config in configs():
        mass = exp(-beta*config_energy(config,energies))
        for i in config:
            etas[i] += mass
        z += mass
        remaining_configs -= 1
        remaining_mass_ub = mass * remaining_configs
        print z,remaining_mass_ub,mass,config,("*" if mass > last_mass else "")
        last_mass = mass
    return map(lambda eta:eta/z,etas)


def recover_matrix(xs,ys):
    """Given a system of matrix equations of the form:

    y = Ax

    solve for A
    """
    dim = len(xs)
    X = np.matrix(xs).transpose()
    Y = np.matrix(ys).transpose()
    # A = np.matrix(([np.linalg.solve(X,np.transpose(ys[i]))
    #                 for i in range(dim)])).transpose()
    return Y*(X.I)

def recover_matrix2(xs,ys):
    """Given a system of matrix equations of the form:

    y = Ax

    solve for A
    """
    dim = len(xs)
    X = np.matrix(xs)
    Y = np.matrix(ys)
    # A = np.matrix(([np.linalg.solve(X,np.transpose(ys[i]))
    #                 for i in range(dim)])).transpose()
    A = np.linalg.solve(X,Y).transpose()
    return A

def test_recover_matrix():
    m = [[random.randrange(0,10) for i in range(10)] for j in range(10)]
    xs = [[random.randrange(0,10) for i in range(10)] for j in range(10)]
    ys = [concat(matrix_mult(m,transpose([x]))) for x in xs]
    m_test = recover_matrix(xs,ys)
    m_test2 = recover_matrix2(xs,ys)

def linearization_recovery_exp():
    """Can we learn a simple, linear relationship between energies and
    etas?"""
    k = 4
    epsilonses = [sorted([random.randrange(10) for i in range(10)]) for j in range(10)]
    etases = [gibbs_more_analytic_opt(epsilons,k) for epsilons in epsilonses]
    new_epsilonses = [[random.randrange(10) for i in range(10)] for j in range(10)]
    new_etases = [gibbs_more_analytic_opt(new_epsilons,k)
                  for new_epsilons in new_epsilonses]
    M1 = recover_matrix(epsilonses,etases)
    M2 = recover_matrix2(epsilonses,etases)
    
def eta_function_recovery_exp():
    energies = range(5)
    k = 2
    iterations = 100
    configs = list(combinations(range(len(energies)),k))
    z = sum(exp(-beta*config_energy(config,energies)) for config in configs)
    kappas_final = [exp(-beta*config_energy(config,energies))/z
                    for config in configs]
    pi = lambda kappas: project_etas(configs,kappas[0])
    matrix = gibbs_matrix_analytic3(energies,k)
    kappas0 = [([0]*(len(configs)-1)) + [1]]
    kappa_traj = iterate_list(lambda kappas:matrix_mult(kappas,matrix),
                              kappas0,iterations)
    eta_traj = map(pi,kappa_traj)
    deta_traj = map(lambda (x,y):zipWith(lambda u,v: v-u,x,y),pairs(eta_traj))
    dfs = [dfdt(energies,etas,k) for etas in eta_traj]
    dgs = [dgdt(energies,etas,k) for etas in eta_traj]
    dhs = [dhdt(energies,etas,k) for etas in eta_traj]
    
def gibbs(energies,n,iterations,beta=beta):
    #initialize
    counts = defaultdict(int)
    bound_states = []
    for i in range(n):
        if i % 10000 == 0:
            print "appending ",i
        bound_states.append(sample_marginal(energies,bound_states))
    #print("initial states:",bound_states)
    #sample
    for i in xrange(iterations):
        if i % 1000 == 0:
            print i
        remove_index = random.choice(bound_states)
        bound_states.remove(remove_index)
        bound_states.append(sample_marginal(energies,bound_states))
        for bs in bound_states:
            counts[bs] += 1
    return [counts[i]/(float(iterations)) if i in counts else 0
            for i in range(len(energies))]
    #return {state:counts[state]/float(iterations) for state in counts}

def project_etas(configs,config_vector):
    n = max(map(max,configs)) + 1
    return [sum(config_prob * (i in config)
                for (config,config_prob) in zip(configs,config_vector))
            for i in range(n)]

def gibbs_trajectory(energies,n,iterations,beta=beta):
    #initialize
    bound_states = []
    history = []
    for i in range(n):
        #if i % 10000 == 0:
            #print "appending ",i
        bound_states.append(sample_marginal(energies,bound_states))
    #print("initial states:",bound_states)
    #sample
    for i in verbose_gen(xrange(iterations),iterations/100):
        if i % 10000 == 0:
            pass #print i
        remove_index = random.choice(bound_states)
        bound_states.remove(remove_index)
        bound_states.append(sample_marginal(energies,bound_states))
        history.append(bound_states[:])
    return history

def gibbs_matrix(energies,n,iterations,beta=beta):
    """Return (some approximation to) the transition matrix given by a
    trajectory"""
    raw_trajectory = gibbs_trajectory(energies,n,iterations,beta=beta)
    traj = map(lambda state: tuple(sorted(state)),raw_trajectory)
    configs = sorted(set(traj))
    d = {config:defaultdict(int) for config in configs}
    for (x,y) in pairs(traj):
        d[x][y] += 1
    matrix = [[d[source_config][dest_config]/float(sum(d[source_config].values()))
               for dest_config in configs]
              for source_config in configs]
    return configs,matrix
            
        

def pick_tf_for_removal(bound_states):
    return random.choice(sum([[e]*bound_states[e] for e in bound_states],[]))

def mb_sample_counts(counts):
    """given a probability distribution expressed as a dictionary of the form
    [energy:degeneracy], sample an energy level"""
    props = {e:exp(-beta*e)*counts[e] for e in counts}
    z = sum(props.values())
    mbs = {e:props[e]/z for e in counts}
    r = random.random()
    total = 0
    for e in mbs:
        total += mbs[e]
        if total > r:
            return e
        
    
def sample_marginals_counts(counts,bound_states):
    unbound_energy_counts = {key: (counts[key] - bound_states[key])
                             for key in counts}
    return mb_sample_counts(unbound_energy_counts)
    
    
def gibbs_counts(counts,n,iterations,beta=beta):
    "What is the probability that the ith energy level is bound, for all i?"
    history = defaultdict(int)
    bound_states = {e:0 for e in counts}
    for i in range(n):
        #print "appending",i
        bound_states[sample_marginals_counts(counts,bound_states)] += 1
    for i in xrange(iterations):
        #print i
        remove_index = pick_tf_for_removal(bound_states)
        #print "chose index for removal:",remove_index
        bound_states[remove_index] -= 1
        add_index = sample_marginals_counts(counts,bound_states)
        bound_states[add_index] += 1
        for bs in bound_states:
            history[bs] += bound_states[bs]
    return {h:history[h]/(counts[h]*float(iterations)) for h in history}

def plot_gibbs_comparison(n,G,iterates = 100000):
    counts = Counter([random.randrange(10) for i in range(G)])
    gibbs_results = gibbs_counts(counts,n,iterates)
    mu = compute_mu_from_energies(counts,n,0.0001)
    plt.plot(gibbs_results.keys(),gibbs_results.values(),linestyle='',marker='o')
    plt.plot(sorted(gibbs_results.keys()),
             fd_probs(sorted(gibbs_results.keys()),mu))
    
def gibbs_counts_framework(ns,Gs,iterates=100000):
    errors = []
    for n in ns:
        for G in Gs:
            if n >= G:
                errors.append(None)
                continue
            counts = Counter([random.randrange(10) for i in range(G)])
            
            gibbs_results = gibbs_counts(counts,n,iterates*n)
            mu = compute_mu_from_energies(counts,n,0.0001)
            error = sum(abs(gibbs_results[key] -
                             fermi_dirac(key,mu)) * counts[key] 
                        for key in gibbs_results)/float(G)
            print "TF:",n,"genome:",G,"error:",error
            errors.append(error)
    errors = np.array(errors).reshape(len(ns),len(Gs))
    return errors

def test_mu_approximation():
    """Test the approximation mu = log(n/Z_0), where Z_0 is the
    partition function of a random genome"""
    pass

def entropy_hunting_exp():
    """Insert a fake genome into E. coli, fish it out via the entropy
    bump method."""
    pass

def expand_counts(counts):
    """Given a dictionary of counts, expand it into a list"""
    return [k for k in counts for i in range(counts[k])]

def zipWith(f,xs,ys):
    return map(lambda (x,y):f(x,y),zip(xs,ys))

def correct (x,y):
    return log(exp(-beta*x) + exp(-beta*y))/-beta

def half_site_exp():
    plt.close()
    tfs = 5
    G = 20
    iterates = int(1e3)
    #fds = sorted([random.random() for i in range(G)])
    fds = [1]*(G-10) + [0]*10
    #bks = sorted([random.random() for i in range(G)])
    bks = [1]*(G-9) + [0]*9
    #bks = [100] * G
    corrects = map(correct,zip(fds,bks))
    # corrects2 = map(lambda (x,y): log((exp(-beta*x) + exp(-beta*y))/2)/-beta,zip(fds,bks))
    # corrects3 = map(lambda (x,y): log((exp(-beta*x) + exp(-beta*y)))/(2*-beta),zip(fds,bks))
    sums = map(lambda (x,y): (x+y),zip(fds,bks))
    avgs = map(lambda (x,y): (x+y)/2.0,zip(fds,bks))
    mins = map(lambda (x,y): min(x,y),zip(fds,bks))
    doubles = gibbs_double_stranded(fds,bks,tfs,iterates)[0][:G]
    singles = gibbs(fds,tfs,iterates)
    exact_correct_probs = (fd_probs_from_energies(corrects,tfs))
    # exact_correct_probs2 = (fd_probs_from_energies(corrects2,tfs))
    # exact_correct_probs3 = (fd_probs_from_energies(corrects3,tfs))
    exact_min_probs = (fd_probs_from_energies(mins,tfs))
    exact_avg_probs = (fd_probs_from_energies(avgs,tfs))
    exact_sum_probs = (fd_probs_from_energies(sums,tfs))
    exact_fd_probs = (fd_probs_from_energies(fds,tfs))
    exact_bk_probs = (fd_probs_from_energies(bks,tfs))
    all_roider_probs = fd_probs_from_energies(fds+bks,tfs)
    fd_roider_probs = all_roider_probs[:G]
    bk_roider_probs = all_roider_probs[G:]
    exact_roider_probs = map(lambda(f,d):f+d-f*d,zip(fd_roider_probs,bk_roider_probs))
    charitable_roider_probs = map(lambda(f,d):f+d-f*d,zip(exact_fd_probs,
                                                          exact_bk_probs))
    linestyles = ["-","--","-.",":"]
    plt.plot(doubles,label="gibbs double_stranded",linestyle=random.choice(linestyles))
    plt.plot(singles,label="gibbs single_stranded",linestyle=random.choice(linestyles))
    plt.plot(exact_correct_probs,label="inverse",linestyle=random.choice(linestyles))
    # plt.plot(exact_correct_probs2,label="inverse2")
    # #plt.plot(exact_correct_probs3,label="inverse3")
    plt.plot(exact_fd_probs,label="exact fd on forward strand",linestyle=random.choice(linestyles))
    plt.plot(exact_min_probs,label="min",linestyle=random.choice(linestyles))
    plt.plot(exact_sum_probs,label="sum",linestyle=random.choice(linestyles))
    plt.plot(exact_avg_probs,label="avg",linestyle=random.choice(linestyles))
    plt.plot(exact_roider_probs,label="roider",linestyle=random.choice(linestyles))
    plt.plot(charitable_roider_probs,label="charitable roider",linestyle=random.choice(linestyles))
    #plt.semilogy()
    fd_vs_single_residuals = zipWith(lambda(x,y):x-y,singles,exact_fd_probs)
    min_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,exact_min_probs)
    avg_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,exact_avg_probs)
    sum_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,exact_sum_probs)
    exact_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,exact_correct_probs)
    roider_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,exact_roider_probs)
    charitable_roider_vs_double_residuals = zipWith(lambda(x,y):x-y,doubles,
                                                    charitable_roider_probs)
    for resid in ["fd_vs_single_residuals",
                  "min_vs_double_residuals",
                  "avg_vs_double_residuals",
                  "sum_vs_double_residuals",
                  "exact_vs_double_residuals",
                  "roider_vs_double_residuals",
                  "charitable_roider_vs_double_residuals"]:
        print "%s\t\t\t mean: %e\t absmean: %e max:%e" %(resid.replace("_residuals",""),
                                                      mean(eval(resid)),
                                                      mean(map(abs,eval(resid))),
                                                      max(eval(resid)))
    plt.legend()
    plt.title("[TF]=%s"%tfs)
    plt.show()

def how_many_geq_than(threshold,rev_sorted_xs):
    count = 0
    for x in rev_sorted_xs:
        if x >= threshold:
            count+=1
        else: return count
    return count

def mask_motif_scores(genome_scores,motif_scores):
    for motif_score in verbose_gen(motif_scores):
        genome_scores.remove(motif_score)
    return genome_scores

def tpr(sorted_motif_scores,threshold):
    tp = how_many_geq_than(threshold,sorted_motif_scores)
    fn = len(sorted_motif_scores) - tp
    return tp / float(tp + fn)

def fpr(sorted_motif_scores,sorted_genome_scores,threshold):
    tp = how_many_geq_than(threshold,sorted_motif_scores)
    fp = how_many_geq_than(threshold,sorted_genome_scores) - tp
    tn = len(sorted_motif_scores) - tp
    return fp / float(fp + tn)

def fpr_and_tpr(sorted_motif_scores,sorted_genome_scores):
    m = []
    for sms in sorted_motif_scores:
        m.append((fpr(sorted_motif_scores,sorted_genome_scores,sms),
                  tpr(sorted_motif_scores,sms)))
        print sms
    return transpose(m)

    
def roc():
    motif = get_crp_motif()
    genome = get_ecoli_genome()
    true_motif = [site for site in verbose_gen(motif)
                  if site in genome or wc(site) in genome]
    fd_motif = [site for site in verbose_gen(true_motif) if site in genome]
    crp = sufficache.PSSM(true_motif)
    fd_crp = sufficache.PSSM(fd_motif)
    naive_motif_scores = sorted([fd_crp.score(site) for site in fd_crp.motif],
                                reverse=True)
    max_motif_scores = sorted([max(crp.score(site),crp.score(wc(site)))
                               for site in crp.motif],reverse=True)
    correct_motif_scores = sorted([correct(crp.score(site),crp.score(wc(site)))
                                   for site in crp.motif],reverse=True)
    fw_scores = crp.slide_score(genome)
    bk_scores = list(reversed(crp.slide_score(wc(genome))))
    max_scores = sorted(zipWith(max,fw_scores,bk_scores), reverse=True)
          
    correct_scores = sorted(zipWith(correct,fw_scores,bk_scores), reverse=True) 
    naive_scores = sorted(fw_scores)
    del fw_scores,bk_scores
    
def r_freq_exp(g,G,epsilon,epsilon_min,temp=37):
    """Schneider's result assumes a perfect recognizer which can avoid
    binding to the genomic background entirely.  How does this result
    vary with the imperfection of the recognizer? """
    T = 273.15 + temp #37C
    G = int(G) #in case G is entered in float notation, i.e. 5e6
    Z = sum([exp(-beta*epsilon) for i in range(g)] +
            [exp(-beta*epsilon_min) for i in range(G)])
    fg_prob = exp(-beta*epsilon)/Z
    bg_prob = exp(-beta*epsilon_min)/Z
    H = (-sum(fg_prob * log2(fg_prob) for i in xrange(g)) +
          -sum(bg_prob * log2(bg_prob) for i in xrange(G)))
    return H

def background_binding(w,G):
    [choose(w,i)*(1/4.0)**i*(3/4.0)**(w-i)*exp(2*i/T) for i in range(w+1)]

def iter_solve(epsilons,etas):
    return[1/(sum(exp(-beta*(ei-ej))*(1-etaj) for (ej,etaj) in zip(epsilons,etas))/
              sum(exp(-beta*(ej-ei))*(etaj) for (ej,etaj) in zip(epsilons,etas)) + 1)
           for ei in etas]

def dfdt(epsilons,etas,k,beta=beta):
    return [sum(exp(-beta*(ei-ej)) * etaj * (1-etai) -
                exp(-beta*(ej-ei)) * etai * (1-etaj)
                for (etaj,ej) in zip(etas,epsilons))
            for (etai,ei) in zip(etas,epsilons)]

def dgdt(epsilons,etas,k,beta=beta):
    z = sum(exp(-beta*ek) * (1-etak) for (etak,ek) in zip(etas,epsilons))
    return [1/z * sum(exp(-beta*ei) * (1-etai) * etaj -
                      exp(-beta*ej) * (1-etaj) * etai
                for (etaj,ej) in zip(etas,epsilons))
            for (etai,ei) in zip(etas,epsilons)]

def dhdt(epsilons,etas,k,beta=beta):
    def flow(ei,ej,etai,etaj):
        adjust = lambda eta: eta*(k-1)/(k-etai)
        z = sum(exp(-beta*ek) * (1-adjust(etak))
                for (etak,ek) in zip(etas,epsilons))
        zp = (z - exp(-beta*ei) * (1-adjust(etai))
                - exp(-beta*ej) * (1-adjust(etaj)))
        zpp = zp + exp(-beta*ej) + exp(-beta*ei)
        i_is_occ  = etai
        j_is_free = 1-adjust(etaj)
        i_goes_to_j = exp(-beta*ej)/zpp
        return i_is_occ*j_is_free * i_goes_to_j
    return [sum(flow(ej,ei,etaj,etai) - flow(ei,ej,etai,etaj)
                for (ej,etaj) in zip(epsilons,etas))
            for (ei,etai) in zip(epsilons,etas)]

def didt(epsilons,etas,k,beta=beta):
    def flow(ei,ej,etai,etaj):
        adjust = lambda eta: eta - (1-eta)*(k-1)/(k-etai)
        z = sum(exp(-beta*ek) * (1-adjust(etak))
                for (etak,ek) in zip(etas,epsilons))
        zp = (z - exp(-beta*ei) * (1-adjust(etai))
                - exp(-beta*ej) * (1-adjust(etaj)))
        zpp = zp + exp(-beta*ej) + exp(-beta*ei)
        i_is_occ  = etai
        j_is_free = 1-adjust(etaj)
        i_goes_to_j = exp(-beta*ej)/zpp
        return i_is_occ*j_is_free * i_goes_to_j
    return [sum(flow(ej,ei,etaj,etai) - flow(ei,ej,etai,etaj)
                for (ej,etaj) in zip(epsilons,etas))
            for (ei,etai) in zip(epsilons,etas)]

def djdt(epsilons,etas,k,beta=beta):
    def flow(ei,ej,etai,etaj):
        adjust = lambda eta: eta
        z = sum(exp(-beta*ek) * (1-adjust(etak))
                for (etak,ek) in zip(etas,epsilons))
        zp = (z - exp(-beta*ei) * (1-adjust(etai))
                - exp(-beta*ej) * (1-adjust(etaj)))
        zpp = zp + exp(-beta*ej) + exp(-beta*ei)
        i_is_occ  = etai
        j_is_free = 1-adjust(etaj)
        i_goes_to_j = exp(-beta*ej)/zpp
        return i_is_occ*j_is_free * i_goes_to_j
    return [sum(flow(ej,ei,etaj,etai) - flow(ei,ej,etai,etaj)
                for (ej,etaj) in zip(epsilons,etas))
            for (ei,etai) in zip(epsilons,etas)]

def dkdt(epsilons,etas,k,beta=beta):
    def flow(ei,ej,etai,etaj):
        "flow from i to j"
        i_occ = etai
        j_free = 1-etaj
        flow_factor = exp(-beta*(ej-ei))
        return i_occ * j_free * flow_factor
        
    return [sum(flow(ej,ei,etaj,etai) - flow(ei,ej,etai,etaj)
                for (ej,etaj) in zip(epsilons,etas))
            for (ei,etai) in zip(epsilons,etas)]

def integrate(epsilons,etas,method,dt,tolerance,beta=beta):
    num_tfs = sum(etas)
    i = 0
    acc = 0
    while True:
        detas = [de * dt for de in method(epsilons,etas,k=num_tfs,beta=beta)]
        old_etas = etas
        etas = zipWith(lambda x,y: x+y,etas,detas)
        l2_diff = sum(zipWith(lambda x,y:(x-y)**2,etas,old_etas))
        i += 1
        acc += l2_diff
        print i, l2_diff, acc
        if l2_diff < tolerance:
            break
    return etas

def print_motif(motif):
    for site in motif:
        print site

def iterate_motif(motif):
    num_sites = len(motif)
    width = len(motif[0])
    pssm = sufficache.PSSM(motif)
    scores = pssm.slide_score(genome)
    indices = sorted(enumerate(scores),key=lambda (x,y):y,reverse=True)[:num_sites]
    new_motif = [genome[i:i+width] for (i,score) in indices]
    return new_motif

def random_genomic_motif(num_sites,width):
    return [(lambda r:genome[r:r+width])(random.randint(0,len(genome)))
                 for i in range(num_sites)]

def motif_convergence_exp(motif=None):
    num_sites = 20
    width = 10
    if motif is None:
        motif = random_genomic_motif(num_sites,width)
    history = {}
    iteration = 0
    while True:
        if motif in history.values():
            past_iteration = [k for k in history if history[k] == motif][0]
            cycle_length = iteration - past_iteration
            print "found cycle of length: ",cycle_length
            return (iteration,cycle_length,motif)
        else:
            history[iteration] = motif
            motif = iterate_motif(motif)
            print_motif(motif)
            print "iteration:",iteration
            iteration += 1

def sliding_window(genome,w):
    for i in range(len(genome)-w+1):
        yield genome[i:i+w]
    
def count_words_of_len(n,genome):
    counts = defaultdict(int)
    for word in sliding_window(genome,n):
        counts[word] += 1
    return counts

def binomial_prob(p,n,k):
    """Return the probability that we obtain k successes out of n
    trials where success occurs w/ prob. p"""
    return choose_large(n,k)*(p**k)*((1-p)**(n-k))

def effective_parameter_exp():
    n = 20
    k = 7
    energies = range(n)
    gibbs_result = gibbs_more_analytic_opt(energies,k)
    mu = compute_mu_from_energies(energies,k,1e-10)
    diffs = [[l2(gibbs_result,fd_probs(energies,mu=mu*m,beta=beta*b))
              for m in [m_raw/100.0 for m_raw in range(200)]] # xaxis
             for b in [b_raw/100.0 for b_raw in range(200)]]# yaxis
    return diffs

def l2_diff(n,k,m,b):
    energies = range(n)
    gibbs_result = gibbs_more_analytic_opt(energies,k)
    mu = compute_mu_from_energies(energies,k,1e-10)
    return l2(gibbs_result,fd_probs(energies,mu=mu*m,beta=beta*b))
    
def min_element_of_matrix(xxs):
    min_diff = None
    for i in range(len(xxs)):
 	for j in range(len(xxs[0])):
 		if min_diff is None or xxs[i][j] < min_diff:
 			min_diff = xxs[i][j]
 			min_i = i
 			min_j = j
 			print "revising",min_i,min_j
    return min_i,min_j
 
def gradient_ascent(f,min_x,max_x,min_y,max_y,x=None,y=None,epsilon = 1e-2):
    if x is None or y is None:
        x = random.random() * (max_x - min_x) + min_x
        y = random.random() * (max_y - min_y) + min_y
    direction = None
    stay = (0,0)
    i = 0
    def cond_print(text):
        if i % 1000 == 0:
            print text
            
    north = (0,epsilon)
    south = (0,-epsilon)
    east = (-epsilon,0)
    west = (epsilon,0)
    while not direction == stay:
        z_cur = f(x,y)
        z_dict = {
            f(*zipWith(lambda x,y:x+y,(x,y),north)) : north,
            f(*zipWith(lambda x,y:x+y,(x,y),south)) : south,
            f(*zipWith(lambda x,y:x+y,(x,y),east)) : east,
            f(*zipWith(lambda x,y:x+y,(x,y),west)) : west,
            z_cur : stay
            }
        z_opt = max(z_dict.keys())
        direction = z_dict[z_opt]
        if direction == stay:
            cond_print(i)
            return (x,y)
        else:
            x,y = zipWith(lambda x,y: x + y,(x,y),direction)
            z_last = z_dict
            df = z_opt - z_cur
            cond_print("df:%s"%df + "z_opt:%s"%z_opt + "(%s,%s) %s" % (x,y,i))
            i += 1

def boltzmann_dist(xs):
    """return the boltzmann distribution on xs"""
    ys = [exp(-beta*x) for x in xs]
    z = sum(ys)
    return [y/z for y in ys]

def delta_g_to_kd(energy):
        return exp(energy/(R*temp)) #km,

def etas_from_energies(energies,p):
    kds = [delta_g_to_kd(energy) for energy in energies]
    return [p/(kd + p) for kd in kds]

def find_free_protein_experimental(energies,copy_number):
    kds = [delta_g_to_kd(energy) for energy in energies]
    n = len(energies)
    k_mean = exp(mean(map(log,kds)))
    
def find_free_protein(energies,copy_number,fp_min=0,fp_max=None,tolerance=1e-10):
    if fp_max == None:
        fp_max = copy_number
    total_protein = lambda fp:sum(etas_from_energies(energies,fp)) + fp - copy_number
    return bisect_interval(total_protein,fp_min,fp_max,tolerance=tolerance)

def find_free_protein_opt(energies,copy_number,fp_min=0,fp_max=None,tolerance=1e-10):
    if fp_max == None:
        fp_max = copy_number
    kds = [delta_g_to_kd(energy) for energy in energies]
    total_protein = lambda fp:sum(fp/(fp + kd) for kd in kds) + fp - copy_number
    return bisect_interval(total_protein,fp_min,fp_max,tolerance=tolerance)


def find_free_protein_opter(energies,copy_number,fp_min=0,fp_max=None,tolerance=1e-10):
    if fp_max == None:
        fp_max = copy_number
    kds = [delta_g_to_kd(energy) for energy in energies]
    total_protein = lambda fp:sum(fp/(fp + kd) for kd in kds) + fp - copy_number
    return secant_interval(total_protein,fp_min,fp_max,tolerance=tolerance)

def find_free_protein_optest(energies,copy_number,fp_min=0,fp_max=None,tolerance=1e-10,p=0.1):
    if fp_max == None:
        fp_max = copy_number
    kds = [delta_g_to_kd(energy) for energy in energies]
    total_protein = lambda fp:sum(fp/(fp + kd) for kd in kds) + fp - copy_number
    return secant_interval_robust(total_protein,fp_min,fp_max,tolerance=tolerance,p=p)


def differential_exp():
    """Can we actually solve this problem just by solving the
    associated system of differential equations? """
    energies = range(100)
    kds = [delta_g_to_kd(energy) for energy in energies]
    def s(i,p):
        return p/(kds[i] + p)
    def dpdt(p):
        return sum(kms[i]*s(i,p) - kp*p*(1-s(i,p)) for i in range(len(energies)))
    
def mu_to_fp(mu):
    return exp(mu/(R * temp))

def fp_to_mu(fp):
    return R * temp * log(fp)

def predict_mu(energies,k,z=None):
    """See Gerland & Hwa 2002, eq 14"""
    if z is None: #occasionally we will wish to compute this repeatedly for large z
        z = sum(exp(-beta*ep) for ep in energies)
    return (log(k) - log(z)) * R*temp

def predict_mu2(energies2,k):
    #What's going on here? Should beta be negative?  Mon Mar  4 11:18:56 EST 2013
    ln_z = sum([log(beta*ep) for ep in energies])
    print ln_z
    z = exp(ln_z)
    return (log(k) - log(z)) * R*temp

def free_energy(energies,fp):
    """This is only an approximation, since we don't consider the
    energy of the free protein'"""
    etas = etas_from_energies(energies,fp) 
    return sum(ep * eta for (ep,eta) in zip(energies,etas))

def free_energy_exact(energies,fp):
    etas = etas_from_energies(energies,fp)
    bound_protein_energy = sum(ep * eta for (ep,eta) in zip(energies,etas))
    free_protein_energy = fp * sum((1-eta) * ep for (ep,eta) in zip(energies,etas))
    return bound_protein_energy - free_protein_energy

def delta_free_energy(energies,k,delta_k = 0.01):
    """A reasonably good approximation to the chemical potential when
    system is dilute, hence more binding sites than copies  """
    fp = find_free_protein(energies,k)
    delta_fp = find_free_protein(energies,k + delta_k)
    return (free_energy(energies,delta_fp) - free_energy(energies,fp))/delta_k

def delta_free_energy_exact(energies,k,delta_k = 0.01):
    delta_k = 0.01
    fp = find_free_protein(energies,k)
    delta_fp = find_free_protein(energies,k + delta_k)
    return (free_energy_exact(energies,delta_fp) -
            free_energy_exact(energies,fp))/delta_k

def find_free_protein2(energies,copy_number):
    kds = [delta_g_to_kd(energy) for energy in energies]
    return copy_number/(sum(1/kd for kd in kds) + 1)
    
def p_eq(energies,p,k):
    kds = [delta_g_to_kd(energy) for energy in energies]
    return sum(p/(kd + p) for kd in kds) + p - k

def simplify_p_eq(energies):
    """Return K that summarizes p_eq"""
    kds = [delta_g_to_kd(energy) for energy in energies]
    n = len(energies)
    f = lambda p: sum(p/(kd + p) for kd in kds)
    ys = map(f,kds)
    g = lambda p,c1,c2:n*(p*c2)/((p*c2) + c1)
    h = lambda c1,c2:sum((y - g(kd,c1,c2))**2 for (y,kd) in zip(ys,kds))
    c = bisect_interval(h,-100,100)

def final_exp(energies):
    kds = [delta_g_to_kd(energy) for energy in energies]
    mean_log_kd = mean(map(log,kds))
    sd_log_kd = sqrt(variance(map(log,kds)))
    phi = lambda x:mlab.normpdf(x,mean_log_kd,sd_log_kd)
    X = lambda: random.gauss(mean_log_kd,sd_log_kd)
    f = lambda p: sum(p/(kd + p) for kd in kds)
    g = lambda x,p: p/(exp(x) + p)
    g_inv = lambda y,p: log(p*(1-y)/y)
    dg_inv_dy = lambda y,p:1/(y-1)
    fy = lambda y,p: phi(log(p*(1-y)/y))*abs(1/(y-1))
    
def expectation_of_sigmoid_of_normal(xs):
    """Solve k = sum p/(p + kd) by treating sum as expectation of
    sigmoid of normal"""
    mu = mean(xs)
    sigma = sqrt(variance(xs))
    def expect(mu,sigma):
        return (1/(exp(-mu) + 1) +
                (exp(-mu) + 1) * exp(-mu)/(2*(exp(-mu) + 1)**3)*sigma**2)

def naive_mu(k):
    """Return log of copy number scaled for biochemical constants.
    Approximation good for very high k.  If using for very low k, subtract """
    return log(k) * R * temp

# March experiments

def graph_no_2(energies=None,control_energies=None,filename=None,plotting=False):
    """The purpose of this experiment is to compare the behavior of
    the analytic and approximated chemical potentials on real and
    synthetic genomes.  The output of the experiment is a graph
    displaying chemical potential vs. copy number for analytic &
    approx. mu on e. coli and synthetic genomes."""
    lb = -5
    ub = 40
    relative_error = 10**-3
    pssm = sufficache.PSSM(Escherichia_coli.Crp)
    copy_numbers = [mantissa * 10**exponent for exponent in range(0,5+1)
                    for mantissa in [1] #for [1,5]
                    if mantissa * 10**exponent <= 10**5]
    mus = []
    approx_mus = []
    control_mus = []
    control_approx_mus = []
    if energies is None:
        genome = get_ecoli_genome()
        energies = pssm.slide_trap(genome)
    if control_energies is None:
        control_genome = random_site(len(genome))
        control_energies = pssm.slide_trap(control_genome)
    for copy_number in copy_numbers:
        mus.append(compute_mu_from_energies(energies,copy_number,relative_error,lb=lb,ub=ub))
        approx_mus.append(predict_mu(energies,copy_number))
        control_mus.append(compute_mu_from_energies(control_energies,
                                              copy_number,relative_error,lb=lb,ub=ub))
        control_approx_mus.append(predict_mu(control_energies,copy_number))
        print("Copy number:",copy_number,mus[-1],approx_mus[-1],
              control_mus[-1],control_approx_mus[-1])
    if plotting:
        plt.plot(copy_numbers,mus,label="E. coli Mu")
        plt.plot(copy_numbers,approx_mus,label="Approximated E. coli Mu")
        plt.plot(copy_numbers,control_mus,label="Control Mu")
        plt.plot(copy_numbers,control_approx_mus,label="Approximated Control Mu")
        plt.xlabel("Copy number")
        plt.ylabel("Chemical potential (kcal/mol + Const)")
        plt.legend(loc=0)
    if filename:
        plt.savefig(filename + "_analytic_vs_approx_genome_vs_control.png")
    return copy_numbers,mus,approx_mus,control_mus,control_approx_mus

def graph_no_2_from_copies(tf_name,mus,exact_copies,approx_copies,
                           exact_control_copies, approx_control_copies,
                           filename=None):
    """The purpose of this experiment is to compare the behavior of
    the analytic and approximated chemical potentials on real and
    synthetic genomes.  The output of the experiment is a graph
    displaying chemical potential vs. copy number for analytic &
    approx. mu on e. coli and synthetic genomes."""
    def clip(copies,mus):
        js = [j for j,copy in enumerate(copies) if 1 <= copy <= 10**6]
        return rslice(copies,js),rslice(mus,js)
    plt.plot(*clip(exact_copies,mus),label="E. coli $\mu$")
    plt.plot(*clip(approx_copies,mus),label="E. coli $\hat{\mu}$")
    plt.plot(*clip(exact_control_copies,mus),label="Control $\mu$")
    plt.plot(*clip(approx_control_copies,mus),label="Control $\hat{\mu}$")
    if tf_name in copy_numbers:
        copy_number = copy_numbers[tf_name]
        plt.plot([copy_number,copy_number],[0,50],linestyle="-.",
                 label="Estimated copy number")
    plt.xlabel("Copy number")
    plt.ylabel("Chemical potential (kcal/mol + Const)")
    plt.semilogx()
    plt.legend(loc=0)
    plt.title(tf_name)
    if filename:
        plt.savefig(filename + "_graph_no_2.png")
        plt.close()
    
def site_killing_exp_all_motifs(copy_dict):
    for tf in Escherichia_coli.tfs:
        motif = getattr(Escherichia_coli,tf)
        motif_length = len(motif)
        copy_number = 10 * motif_length
        exact_copies = copy_dict[tf][0]
        exact_i = min([i for i,copy in enumerate(exact_copies) if copy > copy_number])
        exact_copy = exact_copies[exact_i]
        mu = mus[exact_i]
        print copy_number,exact_copy,mu
        site_killing_exp(motif,tf,energies=None,copy_number=exact_copy,mu=mu)
        
def site_killing_exp(motif,motif_name,energies=None,copy_number=None,mu=None):
    """If we use the approximated vs. analytic mu, how many sites do we lose?"""
    genome = get_ecoli_genome()
    pssm = sufficache.PSSM(motif)
    if energies is None:
        fname = "/home/poneill/binding_sites/dats/%s_traps.dat" % motif_name
        energies = load_array(fname,'f')
    if copy_number is None:
        copy_number = len(motif)
    if mu is None:
        mu = compute_mu_from_energies(energies,copy_number,
                                      relative_error=10**-3,lb=-5,ub=30)
    approx_mu = predict_mu(energies,copy_number)
    motif_energies = [pssm.trap(site) for site in motif]
    probs = [fermi_dirac(energy,mu,beta) for energy in motif_energies]
    approx_probs = [fermi_dirac(energy,approx_mu,beta)
                        for energy in motif_energies]
    lost_sites = len(filter(lambda z:z > 0.1,
                                zipWith(lambda real,approx:real-approx,probs,
                                        approx_probs)))
    lost_percentage = lost_sites/float(len(motif))
    print "lost percentage:",lost_percentage
    if True:
    #plotting annotated
        plt.scatter(probs,approx_probs)
        plt.xlabel("Exact Probability")
        plt.ylabel("Approximated Probability")
        for i,site in enumerate(motif):
            plt.annotate(site.operon,(probs[i],approx_probs[i]))
            plt.title(motif_name + " Binding Sites")
            plt.savefig(motif_name + "_approximate_vs_exact.png")
        plt.close()
    #plotting un-annotated
        plt.scatter(probs,approx_probs)
        plt.xlabel("Probability given exact $\mu$")
        plt.ylabel("Probability given approximate $\mu$")
        plt.title(motif_name + " Binding Sites")
        plt.savefig(motif_name + "_approximate_vs_exact_unannotated.png")
        plt.close()
    return probs,approx_probs

def mcmc_mu(pssm,genome,iterations=1000):
    w = len(pssm)
    n = len(genome) - w + 1
    energies = []
    cur_energy = pssm.trap(genome[:w])
    for iteration in xrange(iterations):
        i = random.randrange(n)
        prop_energy = pssm.trap(genome[i:i+w])
        r = random.random()
        if r < cur_energy/prop_energy:
            cur_energy = prop_energy
            energies.append(cur_energy)
        print cur_energy
    return energies

def load_array(filename,typ):
    n = 5000000 # shitty hard-coded constants
    arr = array(typ)
    with open(filename) as f:
        try:
            arr.fromfile(f,n)
        except:
            pass
    return list(arr)
    
def graph_2_results_for_all_motifs():
    results = []
    for tf in Escherichia_coli.tfs:
        print tf
        fname = "/home/poneill/binding_sites/%s_traps.dat" % tf
        control_fname = "/home/poneill/binding_sites/%s_traps_control.dat" % tf
        try:
            energies = load_array(fname,'f')
            control_energies = load_array(control_fname,'f')
            results.append((tf,graph_no_2(energies,control_energies)))
        except:
            print "failed on:",tf
    return results

def graph_2_results_for_all_motifs_from_copies(copy_dict):
    for tf in Escherichia_coli.tfs:
        print tf
        graph_no_2_from_copies(tf,mus,copy_dict[tf][0],copy_dict[tf][1],
                   copy_dict[tf][2],copy_dict[tf][3], filename=tf+"_ns_4_27_2")

def graph_2_from_results(tf_name,
                         copy_numbers,mus,approx_mus,control_mus,control_approx_mus,plotting=True):
    plt.plot(copy_numbers,mus,label="E. coli $\mu$")
    plt.plot(copy_numbers,approx_mus,label="Approximated E. coli $\mu$")
    plt.plot(copy_numbers,control_mus,label="Control $\mu$")
    plt.plot(copy_numbers,control_approx_mus,label="Approximated Control $\mu$")
    plt.xlabel("Copy number")
    plt.ylabel("Chemical potential (kcal/mol + Const)")
    plt.semilogx()
    plt.title(tf_name)
    plt.legend(loc=0)
    if plotting:
        plt.savefig(tf_name + "_graph_no_2.png")
        plt.close()

ks = [random.gauss(0,1) + 10 for i in range(10)]
f = lambda p:sum(p/float(p+k) for k in ks)

def f_n(n):
    """Return nth derivative of f, evaluated at zero"""
    if n == 0:
        return 0
    else:
        return sum((-1)**(n+1) * factorial(n)/(k**n) for k in ks)

def taylor(x,n):
    """Carry out Taylor expansion of f to n terms"""
    return sum(f_n(i)*x**i/factorial(i) for i in range(n))

def how_many_copies_are_required_exp(motif,energies,mus,plotting=False,real_copies=None,approx_copies=None,title=None,filename=None,alphas=None):
    """Determine how many copies are required to regulate motif to
    desired alpha level using exact, approximate mu"""
    pssm = sufficache.PSSM(motif)
    motif_energies = [pssm.trap(site) for site in motif]
    if alphas is None:
        alphas = []
        for mu in mus:
            motif_probs = [fermi_dirac(energy,mu,beta)
                           for energy in motif_energies]
            alpha = sum(motif_probs)
            alphas.append(alpha)
            print mu,alpha
    if plotting:
        print len(real_copies),len(approx_copies),len(alphas)
        plt.plot(real_copies,alphas,label="$\mu$")
        plt.plot(approx_copies,alphas,label="$\hat{\mu}$")
        plt.xlabel("Copy number")
        plt.ylabel("Total Motif Occupancy")
        plt.legend(loc=0)
        plt.title(title)
        plt.semilogx()
        #plt.show()
        plt.savefig(filename,dpi=400)
        plt.close()
    return alphas

def compute_exact_copies(mus,energies):
    return [show(sum(fd_probs(energies,mu))) for mu in mus]

def compute_approx_copies(mus,energies):
    z = sum([exp(-beta*energy) for energy in energies])
    return map(lambda mu:exp(mu/(R*temp)+log(z)),mus)

def graph_all_how_many_copies_are_required_exp(copy_dict,plotting=True):
    mus = range(-5,50)
    for tf_name in Escherichia_coli.tfs:
        print tf_name
        fname = "/home/poneill/binding_sites/dats/%s_traps.dat" % tf_name
        control_fname = "/home/poneill/binding_sites/dats/%s_traps_control.dat" % tf_name
        try:
            energies = load_array(fname,'f')
            control_energies = load_array(control_fname,'f')
        except e:
            print "couldn't read array for:",tf_name
            raise e
        exact_copies = copy_dict[tf_name][0]
        approx_copies = copy_dict[tf_name][1]
        exact_control_copies = copy_dict[tf_name][2]
        approx_control_copies = copy_dict[tf_name][3]
        motif = getattr(Escherichia_coli,tf_name)
        if plotting:
            how_many_copies_are_required_exp(motif,energies,mus,
                                             True,exact_copies,approx_copies,
                                             tf_name,
                                             tf_name + "_how_many_copies.png")
    return copy_dict

def induction_ratio_exp(copy_dict,alpha_dict):
    """For each TF, compute the induction ratio for mu,mu_hat"""
    real_inductions = []
    approx_inductions = []
    for tf_name in Escherichia_coli.tfs:
        alphas = alpha_dict[tf_name]
        real_copies = copy_dict[tf_name][0]
        approx_copies = copy_dict[tf_name][1]
        lower_cutoff = len(getattr(Escherichia_coli,tf_name)) * 0.05
        upper_cutoff = len(getattr(Escherichia_coli,tf_name)) * 0.95
        lower_i = min([i for i,alpha in enumerate(alphas) if alpha > lower_cutoff])
        upper_i = max([i for i,alpha in enumerate(alphas) if alpha < upper_cutoff])
        lower_real = real_copies[lower_i]
        lower_approx = approx_copies[lower_i]
        upper_real = real_copies[upper_i]
        upper_approx = approx_copies[upper_i]
        real_induction_ratio = upper_real/lower_real
        approx_induction_ratio = upper_approx/lower_approx
        real_inductions.append(real_induction_ratio)
        approx_inductions.append(approx_induction_ratio)
        print tf_name,lower_real,upper_real,lower_approx,upper_approx,real_induction_ratio,approx_induction_ratio
    return real_inductions,approx_inductions

def compute_mode_dict():
    epsilon = 10**-10
    mode_dict = {}
    for tf in Escherichia_coli.tfs:
        regs = [site.regulation for site in getattr(Escherichia_coli,tf)]
        log_mode_ratio = (log(regs.count("activation") + epsilon) -
                          log(regs.count("repression") + epsilon))
 	print tf,log_mode_ratio
 	mode_dict[tf] = log_mode_ratio
    return mode_dict

def make_copy_dict(exact=True,approx=True,ns=False):
    """Make a dictionary of the form:
    {tf_name:(exact_copies,approx_copies,exact_control_copies,approx_control_copies)}
    where the elements of the tuple are lists of copy numbers
    corresponding to various values of mu (range(-5,50)) under the
    four conditions listed.  If ns (non-specific) binding is true,
    subtract 5 kcal/mol from the energy score of each site
    """
    copy_dict = {}
    for tf in Escherichia_coli.tfs:
        print tf
        fname = "/home/poneill/binding_sites/dats/%s_traps.dat" % tf
        control_fname = "/home/poneill/binding_sites/dats/%s_traps_control.dat" % tf
        try:
            energies = load_array(fname,'f')
            control_energies = load_array(control_fname,'f')
        except e:
            print "failed on:",tf
            raise e
        if ns:
            min_ep = min(energies)
            mean_ep = mean(energies)
            print "min_ep:",min_ep
            print "mean_ep:",mean_ep
            print "ns_gap:", mean_ep - min_ep
            motif = getattr(Escherichia_coli,tf)
            width = len(motif[0])
            print "width:",width
            absolute_ns_energy = -8 #kBT = -5 kca/mol
            ep_ns = 2*width + absolute_ns_energy #Assume binding energy is -2kbt/match
            print "ep_ns:",ep_ns
            pssm = sufficache.PSSM(motif)
            site_energies = [pssm.trap(site) for site in motif]
            problem_energies = filter(lambda ep:ep > ep_ns,site_energies)
            print "%s sites have energies greater than ep_ns" % len(problem_energies)
            offset = lambda ep:log(exp(-beta*ep) + exp(-beta*ep_ns))/-beta
            energies = map(offset,energies)
            control_energies = map(offset,control_energies)
        if exact:
            exact_copies = [sum(fd_probs(energies,mu,beta)) for mu in verbose_gen(mus)]
            exact_control_copies = [sum(fd_probs(control_energies,mu,beta))
                                    for mu in verbose_gen(mus)]
        if approx:
            z = sum(exp(-beta*ep) for ep in energies)
            control_z = sum(exp(-beta*ep) for ep in control_energies)
            approx_copies = [approximate_copy_number_from_mu(energies,mu,z)
                             for mu in verbose_gen(mus)]
            approx_control_copies = [approximate_copy_number_from_mu(control_energies,
                                                                     mu,control_z)
                                     for mu in verbose_gen(mus)]
        if not exact:
            copy_dict[tf] = (approx_copies,approx_control_copies)
        elif not approx:
            copy_dict[tf] = (exact_copies,exact_control_copies)
        else:
            copy_dict[tf] = (exact_copies,approx_copies,
                             exact_control_copies,approx_control_copies)
    return copy_dict

def lexa_ns_sanity_check():
    energies = get_energies("LexA")
    control_energies = get_energies("LexA",control=True)
    absolute_ns_energy = -8 #kBT = -5 kca/mol
    width = 16
    ep_ns = 2*width + absolute_ns_energy #Assume binding energy is -2kbt/match
    offset = lambda ep:log(exp(-beta*ep) + exp(-beta*ep_ns))/-beta
    ns_energies = map(offset,energies)
    ns_control_energies = map(offset,control_energies)
    z = sum(exp(-beta*ep) for ep in energies)
    control_z = sum(exp(-beta*ep) for ep in control_energies)
    ns_z = sum(exp(-beta*ep) for ep in ns_energies)
    ns_control_z = sum(exp(-beta*ep) for ep in ns_control_energies)
    exact_copies = [sum(fd_probs(energies,mu,beta)) for mu in verbose_gen(mus)]
    exact_control_copies = [sum(fd_probs(control_energies,mu,beta))
                            for mu in verbose_gen(mus)]
    approx_copies = [approximate_copy_number_from_mu(energies,mu,z)
                     for mu in verbose_gen(mus)]
    approx_control_copies = [approximate_copy_number_from_mu(control_energies,
                                                             mu,control_z)
                             for mu in verbose_gen(mus)]
    ns_exact_copies = [sum(fd_probs(ns_energies,mu,beta)) for mu in verbose_gen(mus)]
    ns_exact_control_copies = [sum(fd_probs(ns_control_energies,mu,beta))
                            for mu in verbose_gen(mus)]
    ns_approx_copies = [approximate_copy_number_from_mu(ns_energies,mu,ns_z)
                     for mu in verbose_gen(mus)]
    ns_approx_control_copies = [approximate_copy_number_from_mu(ns_control_energies,
                                                             mu,ns_control_z)
                             for mu in verbose_gen(mus)]
    

def approximate_copy_number_from_mu(energies,mu,z=None):
    """Solve eq 14 of Gerland & Hwa 2002 for k"""
    if z is None:
        z = sum(exp(-beta*ep) for ep in energies)
    return exp(mu/(R*temp) + log(z))

def get_energies(tf_name,control=False):
    fname = "/home/poneill/binding_sites/dats/%s_traps.dat" % tf_name
    control_fname = "/home/poneill/binding_sites/dats/%s_traps_control.dat" % tf_name
    try:
        if control:
            control_energies = load_array(control_fname,'f')
        else:
            energies = load_array(fname,'f')
    except:
        print "failed on:",tf_name
    return energies if not control else control_energies

def hill_coefficient_exp(tf_name,approx=False):
    """What is the effective hill coefficient of a binding site?"""
    motif = getattr(Escherichia_coli,tf_name)
    pssm = PSSM(motif)
    real_copies = copy_dict[tf_name][0]
    approx_copies = copy_dict[tf_name][1]
    copies = real_copies if not approx else approx_copies
    ns = []
    x_ks = []
    for site in motif:
        site_energy = pssm.trap(site)
        xs = copies
        ys = map(lambda mu:fermi_dirac(site_energy,mu),mus)
        plt.plot(xs,ys)
        x_k,n = fit_hill_function(xs,ys)
        print site,site.operon,n,x_k
        x_ks.append(x_k)
        ns.append(n)
    #plt.semilogx()
    #plt.show()
    return ns,x_ks

def fit_hill_function(xs,ys):
    """Assuming an (increasing) sigmoidal relationship between xs and
    ys, find the Hill coefficient n which best fits the model: y = ymax/(1 + (x_k/x)**n)"""
    ymax = max(ys)
    def error_fun(args):
        [x_k,n] = args
        y_hat = lambda x:ymax/(1 + (x_k/x)**n)
        return sum([(y - y_hat(x))**2 for (x,y) in zip(xs,ys)])
    init_guess = [1,1]
    x_k,n = fmin(error_fun,init_guess)
    return x_k,n

def ns_energy_sanity_check_exp():
    for tf in Escherichia_coli.tfs:
        print tf
        fname = "/home/poneill/binding_sites/dats/%s_traps.dat" % tf
        control_fname = "/home/poneill/binding_sites/dats/%s_traps_control.dat" % tf
        try:
            energies = load_array(fname,'f')
            control_energies = load_array(control_fname,'f')
        except e:
            print "failed on:",tf
            raise e
        min_ep = min(energies)
        mean_ep = mean(energies)
        print "min_ep:",min_ep
        print "mean_ep:",mean_ep
        print "ns_gap:", mean_ep - min_ep
        motif = getattr(Escherichia_coli,tf)
        width = len(motif[0])
        print "width:",width
        absolute_ns_energy = -8 #kBT = -5 kca/mol
        ep_ns = 2*width + absolute_ns_energy #Assume binding energy is -2kbt/match
        print "ep_ns:",ep_ns
        pssm = sufficache.PSSM(motif)
        site_energies = [pssm.trap(site) for site in motif]
        print "mean site energy:",mean(site_energies)
        print "mean diff:",mean_ep - mean(site_energies)
        problem_energies = filter(lambda ep:ep > ep_ns,site_energies)
        print "%s sites have energies greater than ep_ns" % len(problem_energies)

