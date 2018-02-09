import numpy as np
import nucpos_vars
from nucleosome_positioning import NPBackend

def init_NP(order,mechanical_model,temperature) :
    """
    Define a function that inits the NP object which performs the calculations of the 
    probabilities, based on the "order" of the model and the "model". Based on
    the code by Marco Tompitak in his "MarkovModel" package
    """
    # constants
    datadir = '%s/mechanical_models'%(nucpos_vars.data_dir)
    extensions = {1:'nucdist',2:'dinucdist',3:'trinucdist'}
    # Load in the probability tensor for the 'long' oligonucleotides
    try :
        filelong = '%s/%s_%s.%s'%(datadir,mechanical_model,temperature,extensions[order])
        if (order > 1):
            fileshrt = '%s/%s_%s.%s'%(datadir,mechanical_model,temperature,extensions[order-1])
    except IOError :
        raise IOError('Data for mechanical model %s_%s not found'%(mechanical_model,temperature))
    rshptuplong = (148-order,) + (4,)*order
    Pl = np.genfromtxt(filelong).reshape(rshptuplong)
    if (order > 1):
        rshptupshrt = (149-order,) + (4,)*(order-1)
        Ps = np.genfromtxt(fileshrt).reshape(rshptupshrt)
        # Set up the backend
        return NPBackend(order, 147, Pl, Ps)
    else:
        return NPBackend(order, 147, Pl)


_NP = {}

def probability_landscape(seq, order, mechanical_model, temperature) :
    """
    Returns the probability landscape associated to the given order,
    mechanical_model, and temperature. Keeps track of all the NP models that
    were requested, so to avoid unnecessary re-invocations to "init_NP" and keep
    memory tidy
    """
    try :
        return _NP[(order,mechanical_model,temperature)].ProbLandscape(seq)
    except KeyError :
        NP =  init_NP(order, mechanical_model, temperature)
        _NP[(order,mechanical_model,temperature)] = NP
        return NP.ProbLandscape(seq)
