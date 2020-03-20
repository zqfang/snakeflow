from os.path import join, isfile
from itertools import combinations

def unique(seq):
    """Remove duplicates from a list in Python while preserving order.
    :param seq: a python list object.
    :return: a list without duplicates while preserving order.
    """
    seen = set()
    seen_add = seen.add

    return [x for x in seq if x not in seen and not seen_add(x)]

def parse_samples(tab=config['sample_meta']):
    """parse samples """
    SAMPLES=[]
    SAMPLES_ALIAS=[]
    GROUP=[]
    TIME=[]
    with open(tab, 'rU') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not len(line) or line.startswith('#'): continue #skip blank line or comment linne
        item = line.split(" ")
        SAMPLES.append(item[0])
        SAMPLES_ALIAS.append(item[1])
        GROUP.append(item[2])
        if len(item) >3: TIME.append(item[3])

    return SAMPLES, SAMPLES_ALIAS, GROUP, TIME