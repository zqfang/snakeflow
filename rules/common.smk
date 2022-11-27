import os
import re
import sys
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

def parse_samples(tab):
    """parse samples """
    SAMPLES=[]
    SAMPLES_ALIAS=[]
    GROUP=[]
    TIME=[]
    with open(tab, 'rU') as f:
        lines = f.readlines()
    for line in lines:
        if not len(line) or line.startswith('#'): continue #skip blank line or comment linne
        item = line.strip().split()
        SAMPLES.append(item[0])
        SAMPLES_ALIAS.append(item[1])
        GROUP.append(item[2])
        if len(item) >3: TIME.append(item[3])

    return SAMPLES, SAMPLES_ALIAS, GROUP, TIME


def get_sample_fastqs(sample_ids, fastq_dir):
    """This function walk through all possible subfolders to find fastq files given sample name.
    sample_ids: a list of sample ids
    fastq_dir: diectory

    return a tuple contained two dict. 
    - FASTQ_GE: {sample_id: {fastq: [], read_pattern: (),}, ...}
    
    read_pattern e.g. ('HTY-C-B-431A', 'R2', '.fq.gz')
    """
    #READ_PATTERN = "{sample}_{read}_{suffix}"
    pat = re.compile("(.*?)_(\w?\d)([.|_].*)$")

    
    FASTQ_GE = {s:{'fastq':[], 'pattern': []} for s in sample_ids}
    #READ_PATTERN  = {s:[] for s in sample_ids}
    # walk through all files
    for currentpath, folders, files in os.walk(fastq_dir):
        for f in files:
            #print(os.path.join(currentpath, f))
            for s in sample_ids:
                if f.find(s) >= 0: 
                    if f.endswith(("fq", "fq.gz", "fastq", "fastq.gz")):
                        fpath = os.path.join(currentpath, f)
                        FASTQ_GE[s]['fastq'].append(fpath)
                    pat_mat = pat.match(f)
                    if pat_mat:
                        FASTQ_GE[s]['pattern'].append(pat_mat.groups())

    # filter empty samples 
    for s in sample_ids:
        if FASTQ_GE[s]['fastq']:
            # sort files to align read 1 and read 2
            FASTQ_GE[s]['fastq'] = sorted(FASTQ_GE[s]['fastq'])
        else:
            print("Input Sample: %s has not FASTQ files ! Skip."%s, file=sys.stderr)
            del FASTQ_GE[s]

    return FASTQ_GE

