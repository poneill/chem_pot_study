"""
The purpose of this script is to check the analytic results for the
stochastic model to Gillespie simulation.
"""

# TODO: why does time-averaging differ from state-averaging?  Does it?

from random import expovariate
import sys
sys.path.append("../lib/utils")
from utils import *
from params import *
from chem_pot_experiment import fd_probs_from_energies
import itertools

def simulate(ks,p,t_final,time_modulus=1000,record_history=False):
    """Given a list of forward binding constants, ks, and a quantity
    of protein p, perform Gillespie simulation and return a trajectory
    up to time t.  If history is true, return the full trajectory.
    Otherwise, return a list of mean values for each species."""
    n = len(ks)
    # the space is structured as [complexes] + [free protein]
    X = [0] * n + [p]
    t = 0
    last_reported_time = 0
    if record_history:
        history = [(t,X)]
    else:
        mean_X = X[:]
    while t < t_final:
        p_t = X[-1] # current copy number of free protein track
        # propensities for the next reaction to involve the ith site,
        # whether through a binding or unbinding reaction.  If X[i] ==
        # 0, site is unbound and binding reaction occurs w/ propensity
        # p_t * k_i.  Otherwise, unbinding reaction occurs with propensity 1.
        
        propensities = [p_t*ks[i] if not X[i] else 1 for i in range(n)]
        #print "state:",X
        #print "propensities:",propensities
        t_next = expovariate(sum(propensities))
        t = t + t_next
        if t > last_reported_time + time_modulus:
            print "time:",t
            last_reported_time = t
        if record_history:
            history.append((t,X[:]))
        else:
            mean_X = zipWith(lambda x,y:x+y*t_next,mean_X,X)
        if t < t_final:
            r_i = inverse_cdf_sample(range(n),normalize(propensities)) # reaction index
            #print "chose reaction index:",r_i
            if X[r_i] == 0: #if the reaction is a binding event...
                X[r_i] = 1
                X[-1] -= 1
            else: # an unbinding event...
                X[r_i] = 0
                X[-1] += 1
    if record_history:
        return history
    else:
        return map(lambda x:x/float(t_final),mean_X)

def simulate_mean(ks,p,t_final,time_modulus=1000):
    """Given a list of forward binding constants, ks, and a quantity
    of protein p, perform Gillespie simulation and return a trajectory
    up to time t.  If history is true, return the full trajectory.
    Otherwise, return a list of mean values for each species."""
    n = len(ks)
    # the space is structured as [complexes] + [free protein]
    X = [0] * n + [p]
    t = 0
    last_reported_time = 0
    mean_X = X[:]
    last_X = None
    while t < t_final:
        p_t = X[-1] # current copy number of free protein track
        # propensities for the next reaction to involve the ith site,
        # whether through a binding or unbinding reaction.  If X[i] ==
        # 0, site is unbound and binding reaction occurs w/ propensity
        # p_t * k_i.  Otherwise, unbinding reaction occurs with propensity 1.
        
        propensities = [p_t*ks[i] if not X[i] else 1 for i in range(n)]
        #print "state:",X
        #print "propensities:",propensities
        t_next = expovariate(sum(propensities))
        t = t + t_next
        if t > last_reported_time + time_modulus:
            print "time:",t
            last_reported_time = t
        if t < t_final:
            r_i = inverse_cdf_sample(range(n),normalize(propensities)) # reaction index
            #print "chose reaction index:",r_i
            # the following commented line, when replaced by the one that replaces it,
            # results in a 3x speedup.  WTF.
            
            #mean_X = zipWith(lambda x,y:x+y*t_next,mean_X,X) #try to optimize...
            mean_X = [x+y*t_next for x,y in zip(mean_X,X)]
            if X[r_i] == 0: #if the reaction is a binding event...
                X[r_i] = 1
                X[-1] -= 1
            else: # an unbinding event...
                X[r_i] = 0
                X[-1] += 1
    return map(lambda x:x/float(t_final),mean_X)


def average(trajectory):
    """Compute the mean value of each component of the state."""
    t_f = trajectory[-1][0]
    durations = [t1 - t0 for ((t0,X0),(t1,X1)) in pairs(trajectory)]
    mean_vals = map(lambda xs:sum(xs)/t_f,transpose([map(lambda x:x*duration,X)
                          for (duration,(t,X)) in zip(durations,trajectory)]))
    return mean_vals

def average_imp(trajectory):
    """Alternate implementation of average.  Slower, but runs in
    constant space (in length of trajectory)."""
    t_f = float(trajectory[-1][0])
    cum_X = [0] * len(trajectory[0][1])
    t0,X0 = trajectory[0]
    for t1,X1 in trajectory[1:]:
        dt = t1 - t0
        cum_X = [x+y*t_next for x,y in zip(cum_X,X0)]
        #cum_X = zipWith(lambda x,y:x+y*dt,cum_X,X0)
        t0,X0 = t1,X1
    mean_X = map(lambda x:x/t_f,cum_X)
    return mean_X

def recover_states(trajectory):
    """Determine the probability of any given state"""
    state_dict = {}
    for ((t0,X0),(t1,X1)) in pairs(trajectory):
        if not tuple(X0) in state_dict:
            state_dict[tuple(X0)] = 0
        state_dict[tuple(X0)] += (t1 - t0)
    C = sum(state_dict.values())
    return {k:v/C for k,v in state_dict.items()}

def check_state_dict(state_dict,ks,p):
    denom = sum(T(p,q)*sum(product(rslice(ks,comb))
                                 for comb in itertools.combinations(k_range,q))
                                  for q in range(p+1))
    def prob(state):
        print "state:",state
        q = sum(state[:-1])
        constants = [k for (k,j) in zip(ks,state[:-1]) if j]
        print "slice:",constants
        numer = T(p,q)*product(constants)
        return numer/denom
    all_states = state_dict.keys()
    pred_probs = map(prob,all_states)
    obs_probs = [state_dict[state] for state in all_states]
    plt.scatter(pred_probs,obs_probs)
    plt.plot([0,1],[0,1])
    plt.xlabel("Predicted state probability")
    plt.ylabel("Observed state probability")
    
def marginals_from_state_dict(state_dict):
    n = len(state_dict.keys()[0])
    return [sum(states[state] for state in states if state[i] == 1)
            for i in range(n)]
    
def predict_deterministic(ks,p):
    """Given ks and total protein p, predict occupancies from
    deterministic model."""
    kds = [1.0/k for k in ks]
    probs = probs_from_kds(kds,p)
    return probs
    
def predict_stochastic(ks,p):
    """Predict occupancies from the CME using the methods developed in 08/13"""
    k_range = range(len(ks))
    def si_terms(i,q):
        """Return sum of all terms containing ki with exactly q copies bound"""
        return sum([T(p,q)*product(rslice(ks,comb))
                    for comb in itertools.combinations(k_range,q)
                    if i in comb])
    def non_si_terms(i,q):
        """Return sum of all terms containing ki with exactly q copies bound"""
        return sum([T(p,q)*product(rslice(ks,comb))
                    for comb in itertools.combinations(k_range,q)
                    if not i in comb])
    def ith_prob(i):
        """Return the occupancy of the ith site"""
        term = (sum(non_si_terms(i,q) for q in range(p + 1)) /
                sum(si_terms(i,q) for q in range(p + 1)))
        return 1/(1+term)
    return [ith_prob(i) for i in range(len(ks))]

def probs_from_kds(kds,total_protein,epsilon = 10**-10):
    def f(p):
        return sum(p/(kd+p) for kd in kds) + p
    lo, hi = 0, total_protein
    while (abs(hi-lo) > epsilon):
        guess = (hi + lo)/2.0
        total = f(guess)
        if total > total_protein:
            hi = guess
        else:
            lo = guess
    p = guess
    return [p/(kd+p) for kd in kds]

def T(n,k):
    """triangle number"""
    return fac(n)/fac(n-k)

def predict_stochastic_ref(ks,p):
    k_range = range(len(ks))
    print "denom"
    denom = sum(T(p,q)*sum(product(rslice(ks,comb))
                                 for comb in itertools.combinations(k_range,q))
                                 for q in range(p+1))
    print denom
    def numer(i):
        print "numer"
        return sum([T(p,q)*sum(product(rslice(ks,comb))
                                    for q in range(p+1)
             for comb in itertools.combinations(k_range,q)
                                   if i in comb)])
    def numer_test(i):
        #print "numer test"
        return sum(T(p,q)*sum(product(rslice(ks,comb)) if i in comb else 0
                              for comb in itertools.combinations(k_range,q))
                   for q in range(p+1))
    return [numer_test(i)/denom for i in range(len(ks))]

print "loaded"

def compare(xs,ys,color=None,show=True):
    n = min(len(xs),len(ys))
    print n
    if not color:
        plt.scatter(xs[:n],ys[:n])
    else:
        plt.scatter(xs[:n],ys[:n],color=color)
    plt.plot([0,1],[0,1])
    if show:
        plt.show()

def xi(ks,p,i):
    """return quotient of non-si terms over {si terms}/(pki)"""
    k_range = range(len(ks))
    non_si_terms = sum(T(p,q)*sum(product(rslice(ks,comb)) if not i in comb else 0
                              for comb in itertools.combinations(k_range,q))
                   for q in range(p+1))
    si_terms = sum(T(p,q)*sum(product(rslice(ks,comb)) if i in comb else 0
                              for comb in itertools.combinations(k_range,q))
                   for q in range(p+1))
    return non_si_terms/si_terms/(ks[i]*p)

def xi_ref(ks,p,i):
    """return quotient of non-si terms over {si terms}/(pki)"""
    k_range = range(len(ks))
    non_si_terms = sum(T(p,q)*sum(product(rslice(ks,comb)) if not i in comb else 0
                              for comb in itertools.combinations(k_range,q))
                   for q in range(p+1))
    si_terms = sum(T(p,q)*sum(product(rslice(ks,comb)) if i in comb else 0
                              for comb in itertools.combinations(k_range,q))
                   for q in range(p+1))
    return (non_si_terms,si_terms/(ks[i]*p))
    
def test_vieta_formula():
    ks = [2,3,4,5,6]
    n = len(ks)
    K = lambda x:sum(product(comb)*(x**(n-i))
                     for i in range(len(ks) + 1)
                     for comb in itertools.combinations(ks,i))
    print "K(1):",K(1)
    f = lambda x:product([x + k for k in ks])
    print "f(1):",f(1)
