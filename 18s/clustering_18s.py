"""
This script will concern the clustring/categorisation work that we will do
on both the 18S and the SNP data.

To start with, let's visualise the SNP data, through PCoA, then cluster with Kmeans and assess
using the shoulder method and shadow method.
"""

from base_18s import EighteenSBase
from skbio.stats.ordination import pcoa
import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

class Cluster18S(EighteenSBase):
    def __init__(self):
        # Let's pass in the distance class as parent so that we have access to the SNP df already
        super().__init__()
        self.poc_snp_df = self._generate_biallelic_snp_df_and_sample_list('Pocillopora')
        self.por_snp_df = self._generate_biallelic_snp_df_and_sample_list('Porites')
        self.poc_snp_pcoa_df = pcoa(self.poc_snp_df.to_numpy()).samples
        self.por_snp_pcoa_df = pcoa(self.por_snp_df.to_numpy()).samples
    
    def visalise_snp_df(self):
        # first let's just get the pcoA plotted up
        self.fig, self.axar = plt.subplots(1,2, figsize=(5,5))
        
        self.axar[0].scatter(self.poc_snp_pcoa_df['PC1'], self.poc_snp_pcoa_df['PC2'])
        self.axar[1].scatter(self.por_snp_pcoa_df['PC1'], self.por_snp_pcoa_df['PC2'])
        plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster.png'), dpi=600)

        # Now do a clustering and shoulder analysis
        poc_np = self.poc_snp_pcoa_df.to_numpy()
        poc_inertia = []
        for i in range(10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(poc_np)
            poc_inertia.append(kmeans.inertia_)
        
        por_np = self.por_snp_pcoa_df.to_numpy()
        por_inertia = []
        for i in range(10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(por_np)
            por_inertia.append(kmeans.inertia_)
            
        foo = 'bar'

    
    def _generate_biallelic_snp_df_and_sample_list(self, genus):
        if genus == 'Pocillopora':
            snp_dist_path = os.path.join(self.input_dir_18s, 'snp_dist_matrices', 'MP_PocG111_biallelic_gaps.tree.distances.txt')
        elif genus == 'Porites':
            snp_dist_path = os.path.join(self.input_dir_18s, 'snp_dist_matrices', 'MP_PorG111_biallelic_gaps.tree.distances.txt')
        else:
            raise NotImplementedError
        
        # Read in the distance file. There is no header so we will read in as list and then process
        with open(snp_dist_path, 'r') as f:
            dat = [line.rstrip().split('\t') for line in f]
        # Indices are currently in format (e.g. I01S01C011POR). Need to convert to sample-id.
        index = [_[0] for _ in dat]
        index = self._convert_index_to_sample_ids(index)
        dat = np.array([_[1:] for _ in dat], dtype=float).astype(int)
        # make the df
        return pd.DataFrame(data=dat, columns=index, index=index).astype(int)

    def _convert_index_to_sample_ids(self, index):
        # We want to convert these indices to samplie-id
        island_re = re.compile('I\d+')
        site_re = re.compile('S\d+')
        co_re = re.compile('C\d+')
        # A dict that converts from the current sample name (e.g. I01S01C011POR) to the proper sample-id
        # (e.g. TARA_CO-1016606)
        sample_name_dict = {}
        for ind in index:
            island = island_re.search(ind).group()
            site = site_re.search(ind).group()
            colony = co_re.search(ind).group()
            sample_id = self.sample_provenance_df[
                (self.sample_provenance_df['ISLAND#'] == island) & 
                (self.sample_provenance_df['SITE#'] == site) & 
                (self.sample_provenance_df['COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'] == colony) & 
                (self.sample_provenance_df['SAMPLE PROTOCOL LABEL, level 2'] == 'CS4L')].index.values.tolist()
            if len(sample_id) != 1:
                raise RunTimeError('More than one matching sample-id')
            else:
                sample_id = sample_id[0]
            sample_name_dict[ind] = sample_id
        
        # Convert the index to sample-id
        return [sample_name_dict[ind] for ind in index]

Cluster18S().visalise_snp_df()