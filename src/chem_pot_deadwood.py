# Dead code removed from chem_pot_experiment Sat May 18 20:28:50 EDT 2013


def delta_g_from_kd(kd,temp):
    return -R*temp*log(kd)


def index(xs, p):
    """Return index of first x satisfying p, or None if none do"""
    winners = filter(lambda (i, x): p(x), zip(range(len(xs)), xs))
    return winners[0][0] if winners else None


def normalize_dict(counts):
    return dict(zip(counts.keys(),normalize(counts.values())))



def motif_ic(motif):
    return sum(columnwise_ic(motif))



def motif_gini_coefficient(motif):
    ics = sorted(columnwise_ic(motif))
    return gini(ics)


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

def verbose_sum(xs):
    total = 0
    for x in xs:
        total += x
        print total
    return total


def fd(e_i,beta,mu):
    return 1/(exp(beta*(e_i - mu)) + 1)

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



def ghost_probs(energies,num_tfs,beta=beta):
    mb_probs = mb_probs_from_energies(energies,beta)
    return [1-(1-prob)**num_tfs for prob in mb_probs]


def get_copy_number(energies,mu,beta=beta):
    k = 0
    for energy in energies:
        k += fermi_dirac(energy,mu,beta)
    return k


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


def uncountify(counts):
    return sum([[k] * counts[k] for k in verbose_gen(counts,1e5)],[])
    
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


def log10(x):
    return log(x,10)


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
    
def genomic_nmers(genome,n):
    g = len(genome)
    for i in xrange(g - n + 1):
        yield genome[i:i+n]

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


def ith_remainder(cur_config,prop_config,i):
    return sum(cur_config) - prefix_sum(prop_config,i)

def index_by(xs,indices):
    return [xs[i] for i in indices]


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


def update_indices_test():
    n = 10
    k = 3
    init_config = tuple(range(k))
    configs = combinations(range(n),k)
    configs_by_update = iterate_list(lambda c:update_indices(c,n),
                                     init_config,choose(n,k)-1)


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


def test_recover_matrix():
    m = [[random.randrange(0,10) for i in range(10)] for j in range(10)]
    xs = [[random.randrange(0,10) for i in range(10)] for j in range(10)]
    ys = [concat(matrix_mult(m,transpose([x]))) for x in xs]
    m_test = recover_matrix(xs,ys)
    m_test2 = recover_matrix2(xs,ys)


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


def mask_motif_scores(genome_scores,motif_scores):
    for motif_score in verbose_gen(motif_scores):
        genome_scores.remove(motif_score)
    return genome_scores


def fpr_and_tpr(sorted_motif_scores,sorted_genome_scores):
    m = []
    for sms in sorted_motif_scores:
        m.append((fpr(sorted_motif_scores,sorted_genome_scores,sms),
                  tpr(sorted_motif_scores,sms)))
        print sms
    return transpose(m)

    
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


def find_free_protein_experimental(energies,copy_number):
    kds = [delta_g_to_kd(energy) for energy in energies]
    n = len(energies)
    k_mean = exp(mean(map(log,kds)))
    

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


def predict_mu2(energies2,k):
    #What's going on here? Should beta be negative?  Mon Mar  4 11:18:56 EST 2013
    ln_z = sum([log(beta*ep) for ep in energies])
    print ln_z
    z = exp(ln_z)
    return (log(k) - log(z)) * R*temp


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

