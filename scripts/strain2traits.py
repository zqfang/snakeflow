#!/usr/bin/env python
# coding: utf-8
import sys, os
import numpy as np
import pandas as pd

def strain2trait(trait, strain, ids, outdir):
    strains = pd.read_csv(strain)
    traits  = pd.read_csv(trait)
    with open(ids, 'r') as t:
        IDS = t.read().strip().split()

    mpd2strain = {v.loc['MPD'] : v.loc['Abbr'] for k,v in strains.iterrows()}
    traits['strain_abbr'] = traits.strainid.map(mpd2strain) 
    # drop strains could not map to SNP database
    traits.dropna(inplace=True)
    
    # save files
    traits = traits[traits.measnumSex.isin(IDS)]
    traits = traits[traits.strain_abbr.isin(strains['Abbr'])]
    trait_ids = traits.measnumSex.unique()
    if len(trait_ids) != len(np.unique(IDS)):
        sys.exit(f"Error! Some IDs Not Matched! Check file: {ids}")

    for tid in trait_ids:
        temp = traits[traits.measnumSex == tid]
        os.makedirs(os.path.join(outdir,f"MPD_{tid}"), exist_ok=True)
        temp.loc[:,['strain_abbr','strain']].to_csv(os.path.join(outdir,f"MPD_{tid}/strain.{tid}.txt"), index=False, header=None, sep="\t")
        temp.loc[:,['strain_abbr','mean']].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"), index=False, header=None, sep="\t", float_format="%.7f")


strain2trait(snakemake.input['trait'], snakemake.input['strain'], 
             snakemake.input['ids'],   snakemake.params['outdir'])