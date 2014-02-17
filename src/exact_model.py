"""
This script implements gillespie sampling of the heterogeneous
langmuir model in order to compare occupancies to those predicted by
exact math modeling.
"""

from scipy.stats import expon
from utils import inverse_cdf_sample

def simulate(ks,p,tf=10):
    n = len(ks)
    X = [1] * n + [0] * n + [p] #state is a vector of n empty sites, n complexes, protein
    t = 0
    history = [(t,X[:])]
    while t < tf:
        p = X[-1]
        print "X:",X
        rates = [p*ks[i]*X[i] for i in range(n)] + [X[n + i] for i in range(n)]
        print "rates:",rates
        master_rate = float(sum(rates))
        dt = expon.rvs(1,1/master_rate)
        print "dt:",dt
        print "normalized rates:",normalize(rates)
        j = inverse_cdf_sample(range(len(X)),normalize(rates))
        print "chose reaction:",j
        if j < n:
            print "forming complex"
            #update state for complex formation
            X[j] = 0
            X[n+j] = 1
            X[-1] -= 1
        else:
            #update state for complex dissociation
            print "dissolving complex"
            X[j] = 0
            X[j-n] = 1
            X[-1] += 1
        t += dt
        history.append((t,X[:]))
    return history

def occupancies(history):
    """Compute average occupancies over history"""
    tf,Xf = history[-1]
    occupancies = map(sum,transpose([map(lambda x:(t1-t0)/tf*x,X0)
                   for ((t0,X0),(t1,X1)) in pairs(history)]))
    return occupancies
