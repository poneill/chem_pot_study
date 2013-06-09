
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

def graph_2_results_for_all_motifs_from_copies(copy_dict):
    for tf in Escherichia_coli.tfs:
        print tf
        graph_no_2_from_copies(tf,mus,copy_dict[tf][0],copy_dict[tf][1],
                   copy_dict[tf][2],copy_dict[tf][3], filename=tf+"_ns_5_12_")

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
