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
import subprocess

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
        self.fig, self.axar = plt.subplots(2,2, figsize=(8,8))
        
        # self.axar[0,0].scatter(self.poc_snp_pcoa_df['PC1'], self.poc_snp_pcoa_df['PC2'])
        # self.axar[0,1].scatter(self.por_snp_pcoa_df['PC1'], self.por_snp_pcoa_df['PC2'])
        
        # self.axar[0,0].set_title('POC')
        # self.axar[0,1].set_title('POR')

        # Now do a clustering and shoulder analysis
        poc_np = self.poc_snp_pcoa_df.to_numpy()
        poc_inertia = []
        for i in range(10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(poc_np)
            poc_inertia.append(kmeans.inertia_)
        self.axar[0,0].plot([_ + 1 for _ in range(10)], poc_inertia, 'k-')
        self.axar[0,0].set_xticks([_ + 1 for _ in range(10)])
        self.axar[0,0].set_title('Pocillopora')
        self.axar[0,0].set_xlabel('k')
        self.axar[0,0].set_ylabel('inertia')

        por_np = self.por_snp_pcoa_df.to_numpy()
        por_inertia = []
        for i in range(10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(por_np)
            por_inertia.append(kmeans.inertia_)
        self.axar[0,1].plot([_ + 1 for _ in range(10)], por_inertia, 'k-')
        self.axar[0,1].set_xticks([_ + 1 for _ in range(10)])
        self.axar[0,1].set_title('Porites')
        self.axar[0,1].set_xlabel('k')
        self.axar[0,1].set_ylabel('inertia')


        # Elbows would suggest k==5 for POC and k==3 for POR.
        # Plot up a coloured version using the first two components.
        poc_kmeans = KMeans(n_clusters=5, n_init=100, algorithm='full').fit(self.poc_snp_pcoa_df)
        colours = []
        for i in range(poc_kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            scat = self.axar[1,0].scatter(self.poc_snp_pcoa_df.iloc[np.where(poc_kmeans.labels_==i)[0],0], self.poc_snp_pcoa_df.iloc[np.where(poc_kmeans.labels_==i)[0],1], s=16)
            colours.append(scat._original_facecolor[0])
            
        reset_xlim = self.axar[1,0].get_xlim()
        reset_ylim = self.axar[1,0].get_ylim()
        for i in range(poc_kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            # Plot the centroid as vline hline intersect
            #vline
            centroid_x = poc_kmeans.cluster_centers_[i][0]
            centorid_y = poc_kmeans.cluster_centers_[i][1]
            
            self.axar[1,0].plot([centroid_x, centroid_x],[self.axar[1,0].get_ylim()[0], self.axar[1,0].get_ylim()[1]], c=colours[i])
            self.axar[1,0].plot([self.axar[1,0].get_xlim()[0], self.axar[1,0].get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
        # lims need resetting as the plotting of the lines changes them.
        self.axar[1,0].set_xlim(reset_xlim)
        self.axar[1,0].set_ylim(reset_ylim)
        self.axar[1,0].set_title('k=5')
        
        por_kmeans = KMeans(n_clusters=3, n_init=100, algorithm='full').fit(self.por_snp_pcoa_df)
        colours = []
        for i in range(por_kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            scat = self.axar[1,1].scatter(self.por_snp_pcoa_df.iloc[np.where(por_kmeans.labels_==i)[0],0], self.por_snp_pcoa_df.iloc[np.where(por_kmeans.labels_==i)[0],1], s=16)
            colours.append(scat._original_facecolor[0])
        # need to loop again to plot the centroids. We have to do this in a seperate loop so that we can readjust the
        # x and y axis lims.
        reset_xlim = self.axar[1,1].get_xlim()
        reset_ylim = self.axar[1,1].get_ylim()
        for i in range(por_kmeans.cluster_centers_.shape[0]): # for each k
            # Plot the centroid as vline hline intersect
            # vline
            centroid_x = por_kmeans.cluster_centers_[i][0]
            centorid_y = por_kmeans.cluster_centers_[i][1]

            self.axar[1,1].plot([centroid_x, centroid_x],[self.axar[1,1].get_ylim()[0], self.axar[1,1].get_ylim()[1]], c=colours[i])
            self.axar[1,1].plot([self.axar[1,1].get_xlim()[0], self.axar[1,1].get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
            
        # lims need resetting as the plotting of the lines changes them.
        self.axar[1,1].set_xlim(reset_xlim)
        self.axar[1,1].set_ylim(reset_ylim)
        self.axar[1,1].set_title('k=3')

        plt.tight_layout()
        plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster.png'), dpi=600)
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

    def category_scores(self):
        """
        We now have a selection of parameter sets that look to be giving the highest correlation coefficients.
        We will work these initially to generate a set of categorisation agreement scores. I.e. we will
        work with the distance matrices generated from these scores to do clustering and then check the agreement
        of this clustering with the kmeans labels from the SNP biallelic data.

        hard code in a set of parameters to work with, then in series fetch their pcoas and perform kmeans clustering
        on them. We will not yet know what the Ks will be so the first step will be to plot these up.
        Perhaps we can start with one species method combinations (i.e. Porites braycurtis.). Produce a plot that
        is the first components of the pcoa, produce the elbow plot, and we can do this for each of the parameter 
        combinations that we want to try.

        NB it is here that we will want to be testing the with and without the SNP samples for the 18S clustering.
        """
        genus = 'Porites'
        dist_method = 'braycurtis'
        misco_range = np.arange(0.1, 0.6, 0.1)
        masco = 80
        figure, ax = plt.subplots(2, len(misco_range), figsize=(24,8))
        
        for misco_i, misco_v in enumerate(misco_range):
            # Read in the pcoa file of interest as a df
            pcoa_file_name = f'{genus}_True_True_True_False_biallelic_{dist_method}_dist_10000_pwr_False_{misco_v}_{masco}_3_pcoa.csv.gz'
            pcoa_path = os.path.join(self.output_dir_18s, pcoa_file_name)
            subprocess.run(['gzip', '-d', pcoa_path], check=True)
            pcoa_df = pd.read_csv(pcoa_path.replace('.gz', ''))
            pcoa_df = pcoa_df.iloc[:-1, 1:]
            subprocess.run(['gzip', pcoa_path.replace('.gz', '')], check=True)
            
            # Use the df to calculate k-means.
            inertia = []
            for i in range(10):
                kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(pcoa_df)
                inertia.append(kmeans.inertia_)

            # here we have the inertia ready to be plotted up
            ax[0,misco_i].plot([_ + 1 for _ in range(10)], inertia, 'k-')
            ax[0,misco_i].set_xticks([_ + 1 for _ in range(10)])
            ax[0,misco_i].set_title(f'misco: {misco_v}; masco {masco}')
            ax[0,misco_i].set_xlabel('k')
            ax[0,misco_i].set_ylabel('inertia')

            #For the time being just plot up the 2d pcoa on the bottom ax
            ax[1,misco_i].scatter(pcoa_df['PC1'], pcoa_df['PC2'])

            plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster_18s.png'), dpi=600)
            foo = 'bar'

        plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster_18s.png'), dpi=600)
        foo = 'bar'
        
        
        

if __name__ == "__main__":
    c = Cluster18S()
    # c.visalise_snp_df()
    c.category_scores()