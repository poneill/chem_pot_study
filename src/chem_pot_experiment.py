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
from sufficache import PSSM
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

def fast_sample_probs(probs):
    r = random.random()
    total = 0
    for i in range(len(probs)):
        total += probs[i]
        if total > r:
            return i

        
def mutate_site(site):
    pos = random.randrange(len(site))
    return site[:pos] + random.choice("ACGT") + site[pos+1:]
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
    
def compute_mu(b,beta):
    """compute mu from: b = exp(-beta*mu)"""
    return log(b)/-beta

def compute_b(mu,beta):
    """compute b from: b = exp(-beta*mu)"""
    return exp(-beta*mu)

def fermi_dirac(e,mu,beta=beta):
    return 1/(exp((e-mu)*beta) + 1)
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
    
def mb_probs_from_energies(energies,beta=beta):
    props = [exp(-beta*e) for e in energies]
    z = sum(props)
    return [p/z for p in props]

def fd_probs_from_energies(energies,num_tfs,beta=beta,
                           mu_tolerance = 1e-10,
                           approx=False):
    if approx:
        mu = predict_mu(energies,num_tfs)
    else:
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
def fd_probs(energies,mu,beta=beta):
    if type(energies) is list:
        return [fermi_dirac(e,mu,beta) for e in energies]
    elif isinstance(energies,dict):
        return {e:fermi_dirac(e,mu,beta) for e in energies}
    else:
        assert(False)
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
def entropy(probs):
    return -sum(p*log(p,2) for p in probs)
print("loaded chem_pot_experiment")

def countify(energies_or_counts):
    if type(energies_or_counts) is list:
        counts = Counter(energies_or_counts)
    else:
        counts = energies_or_counts
    return counts
def compute_mu_from_energies(energies_or_counts,num_tfs,relative_error,
                             method=fd_bisect,lb=-200,ub=200):
    counts = countify(energies_or_counts)
    def total_prob(mu):
        """sum probability over all sites; should sum to num_tfs;divide by n"""
        return sum(fermi_dirac(energy,mu,beta) * counts[energy]
                   for energy in counts)/float(num_tfs)
    mu = method(lb,ub,num_tfs,relative_error,total_prob)
    return mu
    
def log2(x):
    return log(x,2)
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
    # crp_scale_function = biochem_scale_output(crp)
    # crp_energies = sorted(map(crp_scale_function,crp.slide_score(get_ecoli_genome())))
    # crp_traps = crp.slide_trap(get_ecoli_genome())

def diffs(xs):
    return map(lambda (x,y): x-y,zip(xs+[None],[None]+xs)[1:len(xs)])
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
    
    
def config_energy(config,energies):
    #print config
    return sum([energies[i] for i in config]) #this (explicit list )beats a
                                                  #generator
                                                  #expression by a factor of two!
def prefix_sum(config,i):
    return sum(config[:i])
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

    
def expand_counts(counts):
    """Given a dictionary of counts, expand it into a list"""
    return [k for k in counts for i in range(counts[k])]

def zipWith(f,xs,ys):
    return map(lambda (x,y):f(x,y),zip(xs,ys))

def correct (x,y):
    return log(exp(-beta*x) + exp(-beta*y))/-beta
def how_many_geq_than(threshold,rev_sorted_xs):
    count = 0
    for x in rev_sorted_xs:
        if x >= threshold:
            count+=1
        else: return count
    return count
def tpr(sorted_motif_scores,threshold):
    tp = how_many_geq_than(threshold,sorted_motif_scores)
    fn = len(sorted_motif_scores) - tp
    return tp / float(tp + fn)

def fpr(sorted_motif_scores,sorted_genome_scores,threshold):
    tp = how_many_geq_than(threshold,sorted_motif_scores)
    fp = how_many_geq_than(threshold,sorted_genome_scores) - tp
    tn = len(sorted_motif_scores) - tp
    return fp / float(fp + tn)
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
def sliding_window(genome,w):
    for i in range(len(genome)-w+1):
        yield genome[i:i+w]
    
def l2_diff(n,k,m,b):
    energies = range(n)
    gibbs_result = gibbs_more_analytic_opt(energies,k)
    mu = compute_mu_from_energies(energies,k,1e-10)
    return l2(gibbs_result,fd_probs(energies,mu=mu*m,beta=beta*b))
    
def delta_g_to_kd(energy):
        return exp(energy/(R*temp)) #km,

def etas_from_energies(energies,p):
    kds = [delta_g_to_kd(energy) for energy in energies]
    return [p/(kd + p) for kd in kds]
def find_free_protein(energies,copy_number,fp_min=0,fp_max=None,tolerance=1e-10):
    if fp_max == None:
        fp_max = copy_number
    total_protein = lambda fp:sum(etas_from_energies(energies,fp)) + fp - copy_number
    return bisect_interval(total_protein,fp_min,fp_max,tolerance=tolerance)

def predict_mu(energies,k,z=None):
    """See Gerland & Hwa 2002, eq 14"""
    if z is None: #occasionally we will wish to compute this repeatedly for large z
        z = sum(exp(-beta*ep) for ep in energies)
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
        plt.ylabel("Chemical potential ($K_BT$ + Const)")
        plt.legend(loc=0)
    if filename:
        plt.savefig(filename + "_analytic_vs_approx_genome_vs_control_hi_res.png",dpi=300)
    return copy_numbers,mus,approx_mus,control_mus,control_approx_mus

def clip(copies,mus):
        js = [j for j,copy in enumerate(copies) if 1 <= copy <= 10**6]
        return rslice(copies,js),rslice(mus,js)

def graph_no_2_from_copies(tf_name,mus,exact_copies,approx_copies,
                           exact_control_copies, approx_control_copies,
                           filename=None):
    """The purpose of this experiment is to compare the behavior of
    the analytic and approximated chemical potentials on real and
    synthetic genomes.  The output of the experiment is a graph
    displaying chemical potential vs. copy number for analytic &
    approx. mu on e. coli and synthetic genomes."""
    plt.plot(*clip(exact_copies,mus),label="E. coli $\mu$")
    plt.plot(*clip(approx_copies,mus),label="E. coli $\hat{\mu}$")
    plt.plot(*clip(exact_control_copies,mus),label="Control $\mu$")
    plt.plot(*clip(approx_control_copies,mus),label="Control $\hat{\mu}$")
    if tf_name in copy_numbers:
        copy_number = copy_numbers[tf_name]
        plt.plot([copy_number,copy_number],[0,50],linestyle="-.",
                 label="Estimated copy number")
    plt.xlabel("Copy number")
    plt.ylabel("Chemical potential ($K_BT$ + Const)")
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
    lost_sites = len(filter(lambda z:z > 0.5,
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
            plt.savefig(motif_name + "_approximate_vs_exact.png",dpi=300)
        plt.close()
    #plotting un-annotated
        plt.scatter(probs,approx_probs)
        plt.xlabel("Probability given exact $\mu$")
        plt.ylabel("Probability given approximate $\mu$")
        plt.title(motif_name + " Binding Sites")
        plt.savefig(motif_name + "_approximate_vs_exact_unannotated.png",dpi=300)
        plt.close()
    return probs,approx_probs

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
                   copy_dict[tf][2],copy_dict[tf][3], filename=tf+"_ns_5_12_")

def graph_2_from_results(tf_name,
                         copy_numbers,mus,approx_mus,control_mus,control_approx_mus,plotting=True):
    plt.plot(copy_numbers,mus,label="E. coli $\mu$")
    plt.plot(copy_numbers,approx_mus,label="Approximated E. coli $\mu$")
    plt.plot(copy_numbers,control_mus,label="Control $\mu$")
    plt.plot(copy_numbers,control_approx_mus,label="Approximated Control $\mu$")
    plt.xlabel("Copy number")
    plt.ylabel("Chemical potential ($k_BT$ + Const)")
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

def how_many_copies_are_required_exp(motif,energies,mus,plotting=False,real_copies=None,approx_copies=None,title=None,filename=None,alphas=None,roc=False):
    """Determine how many copies are required to regulate motif to
    desired alpha level using exact, approximate mu"""
    pssm = sufficache.PSSM(motif)
    width = len(motif[0])
    print "width:",width
    absolute_ns_energy = -8 #kBT = -5 kca/mol
    ep_ns = 2*width + absolute_ns_energy #Assume binding energy is -2kbt/match
    print "ep_ns:",ep_ns
    offset = lambda ep:log(exp(-beta*ep) + exp(-beta*ep_ns))/-beta
    motif_energies = map(offset,[pssm.trap(site) for site in motif])
    if alphas is None:
        alphas = [sum([fermi_dirac(energy,mu,beta) for energy in motif_energies])
                  for mu in mus]
    
    if plotting and not roc:
        print len(real_copies),len(approx_copies),len(alphas)
        real_indices = [i for i,v in enumerate(real_copies) if 1 <= v <= 10**6]
        approx_indices = [i for i,v in enumerate(approx_copies) if 1 <= v <= 10**6]
        plt.plot(rslice(real_copies,real_indices),rslice(alphas,real_indices),
                 label="$\mu$")
        plt.plot(rslice(approx_copies,approx_indices),rslice(alphas,approx_indices),
                 label="$\hat{\mu}$")
        plt.xlabel("Copy number")
        plt.ylabel("Total Site Occupancy")
        plt.legend(loc=0)
        plt.title(title)
        plt.semilogx()
        #plt.show()
        plt.savefig(filename,dpi=400)
        plt.close()
    elif plotting and roc:
        n = float(len(motif))
        G = 4639675 # len E. coli K12 MG1655
        real_tprs = [alpha/float(k) for alpha,k in zip(alphas,real_copies)]
        real_fprs = [(n-alpha)/(G - float(k)) for alpha,k in zip(alphas,real_copies)]
        approx_tprs = [alpha/float(k) for alpha,k in zip(alphas,approx_copies)]
        approx_fprs = [(n-alpha)/(G - float(k)) for alpha,k in zip(alphas,approx_copies)]
        plt.plot(real_fprs,real_tprs,label="$\mu$")
        plt.plot(approx_fprs,approx_tprs,label="$\hat{\mu}$")
        plt.xlabel("FPR")
        plt.ylabel("TPR")
        plt.legend(loc=0)
        plt.title(title)
        plt.savefig("roc_" + filename,dpi=400)
        plt.close()
    return alphas



def compute_exact_copies(mus,energies):
    return [show(sum(fd_probs(energies,mu))) for mu in mus]

def compute_approx_copies(mus,energies):
    z = sum([exp(-beta*energy) for energy in energies])
    return map(lambda mu:exp(mu/(R*temp)+log(z)),mus)

def graph_all_how_many_copies_are_required_exp(copy_dict,plotting=True,roc=False):
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
                                             tf_name + "_how_many_copies.png",
                                             roc=roc)
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
    ys, find the Hill coefficient n which best fits the model: y = ymax/(1 + (x_k/x)**n).  Return (x_k, n)"""
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

def predict_last_gibbs_exp(epsilons,k):
    toy_beta = 1
    f = (lambda mu: sum(1/(exp(toy_beta*(ep-mu))+1)  for ep in epsilons))
    mu = bisect_interval(lambda x:f(x)-k,min(epsilons)-1,max(epsilons)+1)
    ps = [1/(exp(toy_beta*(ep-mu)) + 1) for ep in epsilons]
    return ps

def last_gibbs_exp(epsilons,k):
    """Consider a toy system containing two TFs on a graph with sites
    having energies 0, 1, 2.  At each time, the system is in some
    state X denoted as a list of its bound sites.  The TFs may also be
    unbound, so we consider the system to be the set of states along
    with a (cytoplasmic) reservoir which may contain an unlimited
    number of sites.  We perform Gibbs sampling in order to find the
    stationary occupancies for each site."""
    num_sites = len(epsilons)
    toy_beta = 1
    iterations = 100000
    bound_sites = range(k)
    history = [bound_sites]
    for iteration in xrange(iterations):
        unbound_sites = [i for i in range(num_sites) if i not in bound_sites]
        removed_copy = random.choice(bound_sites)
        bound_sites.remove(removed_copy)
        unbound_sites.append(removed_copy)
        probs = normalize([exp(-toy_beta * epsilons[i]) for i in unbound_sites])
        newly_bound_site = inverse_cdf_sample(unbound_sites,probs)
        bound_sites.append(newly_bound_site)
        history.append(bound_sites[:])
    return history

def site_probs_from_history(history,num_sites):
    return [sum(i in state for state in history)/float(len(history))
            for i in range(num_sites)]
        
        
