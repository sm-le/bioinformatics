"""
# Probability 2

Given six positive integers, return expected number of offspring with dominant phenotype

Genotype pairs are as following,
1. AA-AA
2. AA-Aa
3. AA-aa
4. Aa-Aa
5. Aa-aa
6. aa-aa

Example,
IN: 1 0 0 1 0 1
OUT: 3.5
"""

def expected_offspring(genotype_pairing:str) -> float:
    """Expected number of offspring given genotype pairing
    of mendalian unit of heredity.

    Args:
        genotype_pairing: a sequence of a string of genotype pairs

    Returns:
        Expected number of offspring with dominant phenotype
    """
    
    """
    Probability of possessing dominant phenotype,
    1. AA-AA, Pr(dominant) 1
    2. AA-Aa, Pr(dominant) 1
    3. AA-aa, Pr(dominant) 1
    4. Aa-Aa, Pr(dominant) 3/4
    5. Aa-aa, Pr(dominant) 1/2
    6. aa-aa, Pr(dominant) 0

    Expected number of two offspring per couple,
    E(dominant) = Pr(dominant_1) * 2 * no_1_genotype_pair +
                  Pr(dominant_2) * 2 * no_2_genotype_pair +
                  Pr(dominant_3) * 2 * no_3_genotype_pair +
                  Pr(dominant_4) * 2 * no_4_genotype_pair +
                  Pr(dominant_5) * 2 * no_5_genotype_pair +
                  Pr(dominant_6) * 2 * no_6_genotype_pair +
    """

    genotype_pairing = genotype_pairing.split(" ")
    
    # number of offspring
    offsprings = 2

    # Probability of dominant phenotype, 
    pr_dom_1 = 1
    pr_dom_2 = 1
    pr_dom_3 = 1
    pr_dom_4 = 3/4
    pr_dom_5 = 1/2
    pr_dom_6 = 0

    return (int(genotype_pairing[0]) * offsprings * pr_dom_1) + (int(genotype_pairing[1]) * offsprings * pr_dom_2) + \
            (int(genotype_pairing[2]) * offsprings * pr_dom_3) + (int(genotype_pairing[3]) * offsprings * pr_dom_4) + \
            (int(genotype_pairing[4]) * offsprings * pr_dom_5) + (int(genotype_pairing[5]) * offsprings * pr_dom_6)