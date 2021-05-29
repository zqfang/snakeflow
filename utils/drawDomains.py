
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from geneviz.tracks import FeatureTrack, BiomartTrack, plot_tracks
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
        
    def show(self, region=None, show_label=False, figsize=(12, None)):
        # Plot track.
        feat_track = FeatureTrack(data=self.features, hue='type', label='group' if show_label else None)
        if self.features.empty:
            print("No found genes. do search again")
            return
        if region is None:
            region = ('1', 0, self.features.end.max()+1)
        fig, ax = plot_tracks([feat_track], region=region,
                    figsize=figsize, despine=True)
        #self.figure = fig
        return ax

    def add_mutation_feature(self, start, end, label, color="#FF1700"):
        self._features.append(GraphicFeature(start=start, end=end, strand=+1, color=color,
                               label=label))
        self._max_length=max(self._max_length, end)
    
    def add_feature(self, palette='tab10'):
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
        
    def show_feature(self, figure_width = 8, xlabel=""):
        if len(self._features) < 1: 
            print("No feautres to show")
            return
        record = GraphicRecord(sequence_length=self._max_length, features=self._features)
        ax, _ = record.plot(figure_width=figure_width)
        ax.set_xlabel(xlabel, fontweight="bold", fontsize=16)
        return ax
        
    


class VEPTrack:
    def __init__(self, metadata=None, gencode=None, columns=None):
        if isinstance(metadata, pd.DataFrame):
            self.metadata = metadata.copy()
        else:
            self.metadata = pd.read_table(metadata, header=None)
        if columns is not None:
            self.metadata.columns = columns 
        
        self.metadata['strand'] = self.metadata['STRAND'].astype(int)
        self.metadata['Consequence'] = self.metadata.Consequence.str.replace("_variant","")
        tmp = self.metadata.Location.str.split(":", expand=True)
        self.metadata.loc[:,'chromosome'] = tmp.iloc[:,0].astype(str)
        self.metadata.loc[:,'start'] = tmp.iloc[:,1].str.split("-").str[0].astype(int)
        self.metadata.loc[:,'end'] = tmp.iloc[:,1].str.split("-").str[-1].astype(int) + 1
        self._gene_track = BiomartTrack(
                            dataset='mmusculus_gene_ensembl',
                            hue='gene_name',
                            gene_id='gene_name',
                            transcript_id='transcript_id',
                            height=0.5,
                            spacing=0.1,
                            collapse='transcript')
        if isinstance(gencode, pd.DataFrame):
            self.gencode = gencode
        else:
            self.gencode = pd.read_table(gencode, header=None).iloc[:,:6]
            self.gencode.columns = ['chromosome','start','end','name','score','strand']
            self.gencode['chromosome'] = self.gencode['chromosome'].str.replace("chr","")
            self.gencode['name'] = self.gencode['name'].str.split(".").str[0]
            
    def _merge_interval(self, features):
        """
        features should be already sorted. or else will break
        covert start - end to pandas interval, then do the overlap
        """
        ## solution
        df = features.sort_values("start")
        ## This line compares if START of next row is greater than FINISH of current
        ## row ("shift" shifts down FINISH by one row). The value of expression before
        ## cumsum will be True if interval breaks (i.e. cannot be merged), so  
        ## cumsum will increment group value when interval breaks (cum sum treats True=1, False=0)
        df["groupx"]=(df["start"]>df["end"].shift()).cumsum()
        ## this returns min value of "START" column from a group and max value fro m "FINISH"
        result=df.groupby(['chromosome','Consequence','groupx'],as_index=False).agg({"start":"min", "end": "max"})
        if 'strand' not in result.columns:
            result['strand'] = np.nan
        return result
        
        
    def build_strain_tracks(self, transcript_id, strain=None, label=None, height=0.5):
        """
        only altnate genotype will be build
        """
        self._strains = []
        feat_tracks = []
        features = self.metadata.loc[self.metadata.Feature == transcript_id]
        #self._chrom = features.chromosome.values[0]
        if features.empty:
            print(f"transcript: ")
            return
        if strain is None:
            strain = features.IND.unique()
        
        # debug
        #self._features = []
        for s in strain:
            strain_f = features.loc[features.IND == s,] # only alternate
            if strain_f.empty: continue
            strain_f = strain_f.drop_duplicates()
            self._strains.append(s)
            if strain_f.shape[0] >1:
                strain_f = self._merge_interval(strain_f)
            #self._features.append(strain_f)
            if label is None:
                label = 'Consequence'
            ft = FeatureTrack(data=strain_f, height=0.5, hue='Consequence',label=label)# group='name')
            feat_tracks.append(ft)
        self.feature_tracks = feat_tracks
             
    def show_tracks(self, region, feature_tracks=[], figsize=(16,None), padding=[5000, 5000], reverse=False):
        """
        feature_tacks: a list of FeatureTacks
        """
        
        if feature_tracks:
            fts = feature_tracks
        else:
            fts = self.feature_tracks
        # test
        fig, axs = plot_tracks(fts + [self._gene_track],
                           region=region[:3], despine=True, 
                           padding=padding, reverse=reverse,
                           figsize=figsize)
        track_names = self._strains + ['Gene\nStrand '+region[-1]]
        for i, ax in enumerate(axs):
            ax.set_ylabel(track_names[i], rotation=0, fontweight='bold', labelpad=50)
        plt.tight_layout()
        plt.show()
        return fig
    
    def get_region(self, transcript_id):
        """
        find the best region to view 
        """
        temp = self.gencode[self.gencode['name'] == transcript_id]
        if temp.empty:
            print(f"{transcript_id} Not Found! It is correct?")
            return
        #self._gecode_track = FeatureTrack(temp) ## need exon data
        lower = temp.start.min()
        upper = temp.end.max()
        chrom = temp.chromosome.values[0]
        strand = temp.strand.values[0]
        return  (chrom, lower, upper, strand)
    def show(self, transcript_id, strains=None, label=None, figsize=(12, None)):
        self.build_strain_tracks(transcript_id=transcript_id, strain=strains, label=label, height=0.3)
        region = self.get_region(transcript_id)
        #padding_scale = (region[2] - region[1]) / region[1]
        fig= self.show_tracks(region=region, figsize=figsize)
        self.figure = fig
        plt.show()