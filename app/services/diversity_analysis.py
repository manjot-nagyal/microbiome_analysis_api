import numpy as np

def shannon_diversity(counts: np.ndarray) -> float:
    '''

    Calculate Shannon diversity index (H').
    Shannon diversity index takes into account both abundance and evenness of the species present.

    Formula: H' = -sum(p_i * ln(p_i)) for i = 1 to n
    where p_i is the proportion of individuals in the ith species.

    Args:
        data (np.ndarray): A 2D array where each row represents an individual and each column represents a species.

    Returns:
        Shannon diversity index value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    '''
    counts = np.asarray(counts, dtype=float)

    # account for negative abundance values
    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")
    

    counts = counts[counts > 0]
    # consider edge case where there are no counts
    if counts.size == 0:
        return np.nan

    
    proportions = counts / counts.sum()
    return -np.sum(proportions * np.log(proportions)) 
    
    

def simpson_diversity(counts: np.ndarray) -> float:
    """
    Calculate Simpson diversity index (1-D).

    The Simpson diversity index represents the probability that two randomly
    selected individuals in the habitat belong to different species. It ranges
    from 0 (no diversity) to 1 (infinite diversity).

    Formula: 1 - sum(p_i^2) where p_i is the proportion of species i.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Simpson diversity index value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    """
    counts = np.asarray(counts, dtype=float)

    # account for negative abundance values
    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed in diversity calculations")
    

    counts = counts[counts > 0]
    proportions = counts / counts.sum()

    # consider edge case where there are no counts
    if counts.size == 0:
        return np.nan
    
    proportions = counts / counts.sum()
    return float(1.0 - np.sum(proportions ** 2))
    
def pielou_evenness(counts: np.ndarray) -> float:
    '''
    Pielou's evenness measures how evenly the species are distributed in a sample.
    It ranges from 0 (no evenness i.e. no dominance by one species) to 1 (perfect evenness).

    Formula: J' = H' / ln(S) where H' is the Shannon diversity index and S is the number of species.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Pielou's evenness value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    '''
    counts = np.asarray(counts, dtype=float)

    # account for negative abundance values
    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")
    
    if counts.size == 0:
        return np.nan
    
    total = counts.sum()
    if total == 0:
        return np.nan
    
    observed_species = np.count_nonzero(counts)
    # account for when there is only one species
    if observed_species == 1:
        return 1.0
    #     return 0.0

    proportions = counts / total
    # ignore proportions that are too tiny to be considered
    tiny = np.finfo(float).tiny
    pos_proortions = proportions[proportions > tiny] 
    # calculate shannon diversity index
    h_prime = -np.sum(pos_proortions * np.log(pos_proortions)) 
    
    if np.isnan(h_prime) or np.isinf(h_prime):
        return np.nan
    
    # calculate maximum diversity
    max_diversity = np.log(observed_species) if observed_species > 1 else 0.0

    if h_prime == 0.0 or max_diversity == 0.0:
        return 0.0

    return h_prime / max_diversity


def chao1_estimator(counts: np.ndarray) -> float:
    '''
    Calculate Chao1 estimator for species richness.

    Chao1 score is a non-parametric estimation of the number of species in a sample.
    It is based on the number of observed species and the number of rare species.

    Formula: Chao1 = S + [(f1(f1 - 1)) / (2(f2 + 1))]
    where S is the number of observed species, f1 are singletons and f2 are doubletons.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Chao1 estimator value, or NaN for edge cases

    Raises: 
        ValueError: If any count value is negative
    '''
    counts = np.asarray(counts, dtype=float)

    # account for negative abundance values
    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")
    if counts.size == 0:
        return np.nan
    counts = counts[counts > 0]
    S = np.count_nonzero(counts)
    if S == 0:
        return 0.0
    
    f1 = np.count_nonzero(counts == 1)
    f2 = np.count_nonzero(counts == 2)

    if f1 == 0:
        return float(S)
    if f2 == 0:
        return float(S) + (f1 * (f1 - 1)) / 2.0
    
    return float(S) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))
