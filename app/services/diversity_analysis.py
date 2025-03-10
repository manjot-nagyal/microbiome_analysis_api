import numpy as np

def shannon_diversity(counts: np.ndarray) -> float:
    '''
    Calculate the Shannon diversity index for a given dataset. It takes into account
    both abundance and evenness of the species present.

    Formula: H = -sum(p_i * ln(p_i)) for i = 1 to n
    where p_i is the proportion of individuals in the ith species.

    Args:
        data (np.ndarray): A 2D array where each row represents an individual and each column represents a species.

    Returns:
        Shannon diversity index value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    '''
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]
    # consider edge case where there are no counts
    if counts.size == 0:
        return np.nan

    try: 
        proportions = counts / counts.sum()
        return -np.sum(proportions * np.log(proportions)) 
    except (ValueError, ZeroDivisionError, RuntimeWarning):
        return np.nan
    

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
    counts = counts[counts > 0]
    proportions = counts / counts.sum()

    # consider edge case where there are no counts
    if counts.size == 0:
        return np.nan
    try: 
        
        proportions = counts / counts.sum()
        return float(1.0 - np.sum(proportions ** 2))
    except (ValueError, ZeroDivisionError, RuntimeWarning):
        return np.nan
