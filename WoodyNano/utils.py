import numpy as np
import os
import random
import edlib
from pathlib import Path


def change_suffix(fname, suf):
    if suf is None:
        return fname
    apath = Path(fname)
    ss = str(fname).split(apath.suffix)
    return ss[0] + suf

def location_convert_relative_2_absolute(loc, map_start):
    '''
    convert relative location to absolute location
    '''
    if type(loc) is tuple:
        if len(loc) == 2 and loc[0]!=None:
            return tuple([i+map_start for i in loc])
    return 0

def is_empty(c):
    '''
    check length of the input vector
    '''
    if len(c) == 0:
        return True
    else:
        return False

def aligner(primer, sequence, score, map_start=0, map_end=None):
    # do constrained alignment
    if score == -1:
        tmp = edlib.align(query=primer,
                          target=sequence[map_start:map_end],
                          mode='HW',
                          task='locations',
                          k=-1)
    else:
        tmp = edlib.align(query=primer,
                          target=sequence[map_start:map_end],
                          mode='HW',
                          task='locations',
                          k=score*len(primer))

    # convert relative locations to absolute locations
    if not is_empty(tmp['locations']):
        if map_start > 0:
            new_location = [location_convert_relative_2_absolute(i, map_start) for i in tmp['locations']]
            tmp['locations'] = new_location
    else:
        tmp['editDistance'] = 99
    return {'errorRate': tmp['editDistance']/len(primer),
            'locations': tmp['locations'],
            'editDistance': tmp['editDistance']}



def base_complement(k):
    """ Return complement of base.
    Performs the subsitutions: A<=>T, C<=>G, X=>X for both upper and lower
    case. The return value is identical to the argument for all other values.
    :param k: A base.
    :returns: Complement of base.
    :rtype: str
    """
    comp = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
        '-': '-'
    }
    try:

        return comp[k]

    except KeyError:

        print("WARNING: No reverse complement for {} found, returning argument.".format(k))

        return k


def reverse_complement(x):
    """Return reverse complement."""

    a = str()

    for i in x:

        a = base_complement(i)+a

    return a


