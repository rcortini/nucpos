import numpy as np

def energy(probability,kT=1.0) :
    """
    Energy landscape from probability landscape
    """
    return -kT * np.log(probability)
