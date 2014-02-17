import random

beta = 1.61

class System(object):
    """Unified system for simulating TF dynamics through heterogeneous
    Langmuir adsorption"""
    def __init__(self,eps,p):
        """Initialize system from list of forward binding rates,
        protein copy number"""
        self.eps = eps
        self.ks = [exp(-beta*ep) for ep in eps]
        self.xs = [0]*len(self.ks)
        self.p = p
        self.total_p = p
        self.time = 0
        self.history = [0]*len(self.ks)
        self.qZ = None # partition function for marginal probabilities

    def step(self):
        props = [k*self.p if not x else 1 for k,x in zip(self.ks,self.xs)]
        sum_rate = sum(props)
        # pick time of next reaction
        dt = random.expovariate(sum_rate)
        # record history
        self.time += dt
        for i,x in enumerate(self.xs):
            self.history[i] += x * dt
        # end record history
        # do quick and dirty inv_cdf sampling to avoid normalization
        r = random.random()*sum_rate
        #print "r:",r
        acc = 0
        for i,prop in enumerate(props):
            acc += prop
            #print "i,acc:",i,acc
            if acc > r:
                reaction_i = i
                break
        # update the state
        #print "chose reaction:",i
        if self.xs[reaction_i] == 0: # if site is empty...
            self.xs[reaction_i] = 1
            self.p -= 1
        else: # unbinding reaction
            self.xs[reaction_i] = 0
            self.p += 1
        #print self.xs
                
    def mean_occupation(self):
        return [t/self.time for t in self.history]

    def do_steps(self,n):
        for i in xrange(n):
            self.step()

    def configuration_weight(self,xs):
        n = len(xs)
        p_factor = product([p - j for j in range(n)])
        return p_factor * product(rslice(self.ks,xs))

    def set_qZ(self):
        acc = 0
        n = len(self.ks)
        for i in range(self.total_p + 1):
            print i
            for combo in itertools.combinations(range(n),i):
                acc += self.configuration_weight(combo)
        self.qZ = acc
        
    def q(self,i):
        """Compute marginal probabilities analytically"""
        if self.qZ is None:
            self.set_qZ()
        acc = 0
        n = len(self.ks)
        for num_p in range(self.total_p + 1):
            print i
            for combo in itertools.combinations(range(n),num_p):
                if i in combo:
                    acc += self.configuration_weight(combo)
        return acc/self.qZ
        
print "loaded System"
