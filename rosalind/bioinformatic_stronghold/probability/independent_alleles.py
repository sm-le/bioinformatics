"""
# Probability 3

Given k (for generation) and N (number of organism with Aa Bb genotype), return a 
probability of N number of organisms with Aa Bb after kth generation. Assumption here is 
that self-body and its offspring always mate with individual with Aa Bb genotype and result 
two offsprings.  

Example,
IN: 2 1
OUT: 0.684
"""

from math import comb

def independent_allele_v1(generation:int, number:int) -> float:
    """Find probability of given "number" of Aa Bb after "generation"

    Args:
        generation: number of generation
        number: number of success after generation
    """
    
    """
    Probability,

    1st gen,
        AA - BB = (1/4) * (1/4) = 1 / 16
        AA - Bb = (1/4) * (1/2) = 2 / 16
        AA - bb = (1/4) * (1/4) = 1 / 16
        Aa - Bb = (1/2) * (1/2) = 4 / 16 * 
        Aa - BB = (1/2) * (1/4) = 2 / 16
        Aa - bb = (1/2) * (1/4) = 2 / 16
        aa - BB = (1/4) * (1/4) = 1 / 16
        aa - Bb = (1/4) * (1/2) = 2 / 16
        aa - bb = (1/4) * (1/4) = 1 / 16

    2nd gen,
        AA - BB,
            AA - BB = (1/2) * (1/2) = 4 / 16
            AA - Bb = (1/2) * (1/2) = 4 / 16
            AA - bb = (1/2) * 0 = 0
            Aa - Bb = (1/2) * (1/2) = 4 / 16 * 
            Aa - BB = (1/2) * (1/2) = 4 / 16
            Aa - bb = (1/2) * 0 = 0
            aa - BB = 0 * (1/2) = 0
            aa - Bb = 0 * (1/2) = 0
            aa - bb = 0 * 0 = 0

        AA - Bb,
            Aa - Bb = (1/2) * (1/4) = 4 / 16
        AA - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - Bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - BB,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - BB,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - Bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
    
    Summary,
        All probability leading to Aa-Bb genotype is 1/4 or 0.25. In 
        probability, the binomial distribution with param n and p is the discrete
        probability distribution of the number of successes in a sequence of n
        indepedent experiments,

        Corresponding formula is 
            P_{x} = (n x) * p^{x} * q^{n-x}
            (n x) = number of combinations
            x = number of times for a specific outcome within n trials
            p = probability of success on a single trial
            q = probability of failure on a single trial
            n = number of trials (two offspring and given as 2^generation)

        Here, we already stated that p is 0.25 and q is 1-0.25
    """

    # cum_prob = 0
    # for i in range(0, number):
    #     combination = comb(2**generation, i)
    #     raw_probability = ((0.25)**i) * ((0.75)**(2**generation-i))
    #     cum_prob += combination * raw_probability

    cum_prob = sum([comb(2**generation, i) * ((0.25)**i) * ((0.75)**(2**generation-i)) for i in range(number)])

    return 1 - cum_prob

def independent_allele_v2(generation:int, number:int) -> float:
    """Find probability of given "number" of Aa Bb after "generation"

    Args:
        generation: number of generation
        number: number of success after generation
    """
    
    """
    Probability,

    1st gen,
        AA - BB = (1/4) * (1/4) = 1 / 16
        AA - Bb = (1/4) * (1/2) = 2 / 16
        AA - bb = (1/4) * (1/4) = 1 / 16
        Aa - Bb = (1/2) * (1/2) = 4 / 16 * 
        Aa - BB = (1/2) * (1/4) = 2 / 16
        Aa - bb = (1/2) * (1/4) = 2 / 16
        aa - BB = (1/4) * (1/4) = 1 / 16
        aa - Bb = (1/4) * (1/2) = 2 / 16
        aa - bb = (1/4) * (1/4) = 1 / 16

    2nd gen,
        AA - BB,
            AA - BB = (1/2) * (1/2) = 4 / 16
            AA - Bb = (1/2) * (1/2) = 4 / 16
            AA - bb = (1/2) * 0 = 0
            Aa - Bb = (1/2) * (1/2) = 4 / 16 * 
            Aa - BB = (1/2) * (1/2) = 4 / 16
            Aa - bb = (1/2) * 0 = 0
            aa - BB = 0 * (1/2) = 0
            aa - Bb = 0 * (1/2) = 0
            aa - bb = 0 * 0 = 0

        AA - Bb,
            Aa - Bb = (1/2) * (1/4) = 4 / 16
        AA - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - Bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - BB,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        Aa - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - BB,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - Bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
        aa - bb,
            Aa - Bb = (1/2) * (1/2) = 4 / 16
    
    Summary,
        All probability leading to Aa-Bb genotype is 1/4 or 0.25. In 
        probability, the binomial distribution with param n and p is the discrete
        probability distribution of the number of successes in a sequence of n
        indepedent experiments,

        Corresponding formula is 
            P_{x} = (n x) * p^{x} * q^{n-x}
            (n x) = number of combinations
            x = number of times for a specific outcome within n trials
            p = probability of success on a single trial
            q = probability of failure on a single trial
            n = number of trials (two offspring and given as 2^generation)

        Here, we already stated that p is 0.25 and q is 1-0.25
    """

    # cum_prob = 0
    # for i in range(number, 2**generation):
    #     combination = comb(2**generation, i)
    #     raw_probability = ((0.25)**i) * ((0.75)**(2**generation-i))
    #     cum_prob += combination * raw_probability

    cum_prob = sum([comb(2**generation, i) * ((0.25)**i) * ((0.75)**(2**generation-i)) for i in range(number, (2**generation)+1)])

    return cum_prob