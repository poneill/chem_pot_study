"""This file contains library code for the final draft of the chemical
potential paper."""

# The first step is to generate, for each TF, a list of binding
# energies to model the binding landscape of the chromosome.

# Next, we construct mu-k diagrams.  This requires summing occupancies
# for a given value of mu in order to recover the copy number.
# (Conceptually, we consider mu a function of copy number, not
# vice-versa, but it is computationally more efficient to compute the
# inverse function (viz: mu->k))

# To compute mu-mu_hat diagrams, we must compare binding probabilities
# for each known binding site at physiological copy number.

# The misclassification rate diagram is just a summary of the mu-mu_hat diagrams.

# the k-occupancy diagrams just generalize mu-mu_hat diagrams to
# arbitrary copy number.
