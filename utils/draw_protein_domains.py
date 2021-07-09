
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from bioservices.uniprot import UniProt
#pd.options.mode.chained_assignment = None  # default='warn'
from dna_features_viewer import GraphicFeature, GraphicRecord


class GeneVis:
    def __init__(self,):
        self.uniprot = UniProt()
        self.family_domains_columns = ['comment(DOMAIN)', 'comment(SIMILARITY)',
                                       'families', 'feature(COILED COIL)',
                                       'feature(COMPOSITIONAL BIAS)', 'feature(DOMAIN EXTENT)', 
                                       'feature(MOTIF)', 'feature(REGION)',
                                       'feature(REPEAT)','feature(ZINC FINGER)' ]
    def ensembl2Uniprot(self, ensembl_txid,):
        d = u.mapping(fr="ENSEMBL_TRS_ID", to="ACC", query=ensembl_txid)
        
  
    def search(self, uniprot_kw="Nid1_MOUSE"):
        #for col in family_domains_columns:
        results = self.uniprot.search(uniprot_kw, columns=','.join(self.family_domains_columns))
        self._results = results
        comm_pat = '(\d+)[.]+(\d+);[\t ]+/note=([0-9a-zA-Z"\- ]+);'
        domain_pat = re.compile('DOMAIN ' + comm_pat)
        motif_pat = re.compile('MOTIF ' + comm_pat)
        repeat_pat = re.compile('REPEAT ' + comm_pat)
        region_pat = re.compile('REGION ' + comm_pat)

        temp =  []
        for pat in [domain_pat, motif_pat, repeat_pat, region_pat]:
            ## convert to pandas df
            ext = pat.findall(results)
            if ext:
                temp += ext

        temp2 = pd.DataFrame.from_records(temp, columns=('start','end','group'))
        temp2.insert(0, column='chromosome', value='1')
        temp2.insert(3, column='strand', value=None)
        temp2['group'] = temp2['group'].str.strip('"')
        temp2['type'] = temp2['group'].str.replace(r" \d", "", regex=True)
        
        temp2.start = temp2.start.astype(int)
        temp2.end = temp2.end.astype(int)
        if temp2.empty:
            print(f"Warning: NO features found for {uniprot_kw}")
        self.features = temp2
        self._max_length = int(self.uniprot.search(uniprot_kw, columns='length').split("\n")[1])
        

    def add_mutation_feature(self, start, end, label, color="#FF1700"):
        self._features.append(GraphicFeature(start=start, end=end, strand=+1, color=color,
                               label=label))
        self._max_length=max(self._max_length, end)
    
    def add_domains(self, palette='tab10'):
        if self.features.empty:
            print("No found genes. do search again")
            return
        self._features = []
        ft = self.features['type'].unique()
        colors = sns.color_palette(palette=palette, n_colors=len(ft)).as_hex()
        self.features['color'] = self.features['type'].map({t:c for t,c in zip(ft, colors)})
        
        for i, row in self.features.iterrows():
            f = GraphicFeature(start=row.start, end=row.end, strand=+1, color=row.color,
                               label=row.group)
            self._features.append(f)
        
    def show_feature(self, ax=None, figure_width = 8, xlabel=""):
        if len(self._features) < 1: 
            print("No feautres to show")
            return
        record = GraphicRecord(sequence_length=self._max_length, features=self._features)
        ax2, _ = record.plot(ax=ax, figure_width=figure_width)
        ax2.set_xlabel(xlabel, fontweight="bold", fontsize=16)
        return ax2


if __name__ == "__main__":
    

    domains = GeneVis()
    # query Uniprot to get domains
    symbol = "Cubn"
    domains.search("Q9JLB4") # Q3URU2 Uniprot ID
    # draw protein domains
    domains.add_domains()
    
    # add mutation info: 
    domains.add_mutation_feature(start=437,end=437+1, 
                                 label=f"mutation: V/M, pos: 437")
    domains.add_mutation_feature(start=3562,end=3562+1, 
                                 label=f"mutation: G/D, pos: 3562")
    domains.add_mutation_feature(start=402,end=402+1, 
                                 label=f"mutation: H/Y, pos: 402")
    fig, ax =  plt.subplots(figsize=(20,4))
    domains.show_feature(ax=ax, figure_width=20, xlabel=symbol)
    fig.savefig("Q9JLB4.pdf")
    # ax = domains.show_feature(figure_width=20, xlabel=symbol)
