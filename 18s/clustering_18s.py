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
from subprocess import CalledProcessError
from sklearn.metrics import silhouette_score
import compress_pickle
import itertools
import sys
import random

class Cluster18S(EighteenSBase):
    def __init__(self):
        # Let's pass in the distance class as parent so that we have access to the SNP df already
        super().__init__()
        self.poc_snp_df = self._generate_biallelic_snp_df_and_sample_list('Pocillopora')
        self.por_snp_df = self._generate_biallelic_snp_df_and_sample_list('Porites')
        self.poc_snp_pcoa_df = self._generate_biallelic_snp_pcoa_df('Pocillopora')
        self.por_snp_pcoa_df = self._generate_biallelic_snp_pcoa_df('Porites')
        self.poc_snp_kmeans_dict = {}
        self.poc_snp_kmeans_dict_pickle_path = os.path.join(self.cache_dir, 'poc_biallelic_snp_kmeans_dict.p.bz')
        self.por_snp_kmeans_dict = {}
        self.por_snp_kmeans_dict_pickle_path = os.path.join(self.cache_dir, 'por_biallelic_snp_kmeans_dict.p.bz')
        
        # Use the meta_info_table produced by output_tables.py
        self.meta_info_table_path = os.path.join(self.output_dir, 'coral_18S_meta_info_table_2020-04-15T12_53_51.189071UTC.csv')
        self.meta_info_df = self._get_meta_info_df()
    
    def _get_meta_info_df(self):
        return pd.read_csv(self.meta_info_table_path).set_index('readset')

    def _generate_biallelic_snp_pcoa_df(self, genus):
        if genus == 'Pocillopora':
            cache_path = os.path.join(self.cache_dir_18s, 'poc_biallelic_snp_pcoa_df.p.bz')
            try:
                return pd.read_pickle(cache_path)
            except FileNotFoundError:
                poc_snp_pcoa_df = pcoa(self.poc_snp_df.to_numpy()).samples
                poc_snp_pcoa_df.index = self.poc_snp_df.index
                poc_snp_pcoa_df.to_pickle(cache_path)
                return poc_snp_pcoa_df
        else:
            cache_path = os.path.join(self.cache_dir_18s, 'por_biallelic_snp_pcoa_df.p.bz')
            try:
                return pd.read_pickle(cache_path)
            except FileNotFoundError:
                por_snp_pcoa_df = pcoa(self.por_snp_df.to_numpy()).samples
                por_snp_pcoa_df.index = self.por_snp_df.index
                por_snp_pcoa_df.to_pickle(cache_path)
                return por_snp_pcoa_df

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
        poc_sil = {}
        # Hold the kmeans results in a dict and pickle out this dict so that we can use it
        # in comparing between the 18S and SNP categorisation agreements
        for i in range(1,10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(poc_np)
            poc_inertia.append(kmeans.inertia_)
            poc_sil[i + 1] = silhouette_score(self.poc_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
            self.poc_snp_kmeans_dict[i + 1] = kmeans
        
        compress_pickle.dump(self.poc_snp_kmeans_dict, self.poc_snp_kmeans_dict_pickle_path)

        self.axar[0,0].plot([_ + 1 for _ in range(1,10)], poc_inertia, 'k-')
        self.axar[0,0].set_xticks([_ + 1 for _ in range(1,10)])
        self.axar[0,0].set_title('Pocillopora')
        self.axar[0,0].set_xlabel('k')
        self.axar[0,0].set_ylabel('inertia')
        
        
        ax2 = self.axar[0,0].twinx()
        ax2.plot([_ + 1 for _ in range(1,10)], [poc_sil[_ + 1]  for _ in range(1,10)], 'r-')
        ax2.set_ylabel('silhouette score')
        self.axar[0,0].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')


        por_np = self.por_snp_pcoa_df.to_numpy()
        por_inertia = []
        por_sil = {}
        for i in range(1,10):
            kmeans = KMeans(n_clusters=i + 1, n_init=100, algorithm='full').fit(por_np)
            por_inertia.append(kmeans.inertia_)
            por_sil[i + 1] = silhouette_score(self.por_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
            self.por_snp_kmeans_dict[i + 1] = kmeans

        compress_pickle.dump(self.por_snp_kmeans_dict, self.por_snp_kmeans_dict_pickle_path)

        self.axar[0,1].plot([_ + 1 for _ in range(1, 10)], por_inertia, 'k-')
        self.axar[0,1].set_xticks([_ + 1 for _ in range(1, 10)])
        self.axar[0,1].set_title('Porites')
        self.axar[0,1].set_xlabel('k')
        self.axar[0,1].set_ylabel('inertia')

        ax2 = self.axar[0,1].twinx()
        ax2.plot([_ + 1 for _ in range(1, 10)], [por_sil[_ + 1] for _ in range(1, 10)], 'r-')
        ax2.set_ylabel('silhouette score')
        self.axar[0,1].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

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


    # TODO. This needs refactoring. But let's come back to this so that we can 
    # investigate the more important issues here.
    
    def _load_snp_kmeans_dicts(self):
        poc_snp_dict = compress_pickle.load(self.poc_snp_kmeans_dict_pickle_path)
        por_snp_dict = compress_pickle.load(self.por_snp_kmeans_dict_pickle_path)
        return poc_snp_dict, por_snp_dict

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
        
        self.poc_snp_kmeans_dict, self.por_snp_kmeans_dict = self._load_snp_kmeans_dicts()
        
        for genus in ['Pocillopora', 'Porites']:
            for dist_method in ['braycurtis', 'unifrac']:
                if genus == 'Porites' and dist_method == 'braycurtis':
                    misco_range = [f'{num:.1f}' for num in np.arange(0.1, 0.6, 0.1)]
                    masco_range = [80]
                elif genus == 'Porites' and dist_method == 'unifrac':
                    misco_range = ['0.02', '0.66']
                    masco_range = [50,100]
                elif genus == 'Pocillopora' and dist_method == 'braycurtis':
                    misco_range = ['0.08']
                    masco_range = range(100,250,50)
                elif genus == 'Pocillopora' and dist_method == 'unifrac':
                    #DONE
                    misco_range = ['0.08']
                    masco_range = range(100,300,50)
                
                
                num_param_combs = len(misco_range) * len(masco_range)
                figure, ax = plt.subplots(4, num_param_combs, figsize=(4*num_param_combs, 12))
        
                # Make a figure set for each of the misco values
                for misco_i, misco in enumerate(misco_range):
                    for masco_i, masco in enumerate(masco_range):
                        occ = OneClusterCol(genus=genus, dist_method=dist_method, misco=misco, masco=masco, ax=ax[:,(misco_i * len(masco_range)) + masco_i], parent=self)
                        occ.plot_col()

                # Once plotting is complete for all parameter combinatinos. Write out fig.
                plt.tight_layout()
                plt.savefig(os.path.join(self.eighteens_dir, f'temp_fig_cluster_{genus}_{dist_method}_18s.png'), dpi=600)
                foo = 'bar'

    def check_old_clustering_agreement(self):
        """
        We want to read in the old clustering assignments for Pocillopora
        that we were able to get high agreement with.
        We want to see how these correspond to the new SNP PCoA and to the old
        Plot up all points that don't have an annotation in the old set as grey 
        and the those that do with a colour corresponding to their categorisation.
        DONE
        """
        old_df = pd.read_csv("/home/humebc/projects/tara/tara_full_dataset_processing/18s/input/snp_dist_matrices/pocillopora_snp_clusters.csv")
        old_df.drop(columns='Colony', inplace=True)
        old_df.rename(columns={'code pangea': 'sample-id', 'snmf/PCA':'label'}, inplace=True)
        old_df.set_index('sample-id', inplace=True, drop=True)
        poc_snp_kmeans_dict = compress_pickle.load(os.path.join(self.cache_dir, 'poc_biallelic_snp_kmeans_dict.p.bz'))
        poc_kmeans = poc_snp_kmeans_dict[3]
        poc_label_series = pd.Series(poc_kmeans.labels_, index=self.poc_snp_df.index, name='label')
        in_common_sample_ids = [_ for _ in old_df.index if _ in poc_label_series.index]
        # I think the fastest and clearest way to display this is to plot up the 
        # poc SNP pcoa in 2D and colour (according to categories in old_df) the 27 in common samples
        # all others should be grey
        fig, ax = plt.subplots(3,1, figsize=(4,12))
        # first plot up in grey
        ax[0].scatter(self.poc_snp_pcoa_df.loc[[_ for _ in poc_label_series.index  if _ not in old_df.index],'PC1'], self.poc_snp_pcoa_df.loc[[_ for _ in poc_label_series.index  if _ not in old_df.index],'PC2'], c='lightgrey')
        # Then plot up each of the ks
        for i in range(3):
            samples_to_plot = old_df[old_df['label']==i+1].index
            ax[0].scatter(self.poc_snp_pcoa_df.loc[samples_to_plot,'PC1'], self.poc_snp_pcoa_df.loc[samples_to_plot,'PC2'])

        ax[0].set_title('Old POC SNP classifications')
        ax[0].set_xlabel('PC1')
        ax[0].set_xlabel('PC2')
        
        
        # Now do two plots below one same as above but within the context of the 18S pcoas
        # one for each dist method.
        # FOr the Poc braycurtis use misco 0.08 masco 200.
        # For the poc unifrac use misco 0.08 masco 200
        
        # First BC
        pcoa_file_name = f'Pocillopora_True_True_True_False_biallelic_braycurtis_dist_10000_pwr_False_0.08_200_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax[1].scatter(pcoa_df.loc[[_ for _ in pcoa_df.index  if _ not in old_df.index], 'PC1'], pcoa_df.loc[[_ for _ in pcoa_df.index  if _ not in old_df.index],'PC2'], c='lightgrey')
        for i in range(3):
            samples_to_plot = old_df[old_df['label']==i+1].index
            ax[1].scatter(pcoa_df.loc[samples_to_plot,'PC1'], pcoa_df.loc[samples_to_plot,'PC2'])

        # Then Unifrac
        pcoa_file_name = f'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_pwr_False_0.08_200_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax[2].scatter(pcoa_df.loc[[_ for _ in pcoa_df.index  if _ not in old_df.index], 'PC1'], pcoa_df.loc[[_ for _ in pcoa_df.index  if _ not in old_df.index],'PC2'], c='lightgrey')
        for i in range(3):
            samples_to_plot = old_df[old_df['label']==i+1].index
            ax[2].scatter(pcoa_df.loc[samples_to_plot,'PC1'], pcoa_df.loc[samples_to_plot,'PC2'])
        
        plt.savefig(os.path.join(self.eighteens_dir, 'temp_old_poc_snp_classification.png'), dpi=600)
        foo = 'bar'

    def _get_pcoa_df(self, pcoa_path=None):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        if pcoa_path is None:
            if self.dist_method == 'unifrac':
                pcoa_file_name = f'{self.genus}_True_True_True_False_biallelic_{self.dist_method}_dist_1000_pwr_False_{self.misco}_{self.masco}_3_pcoa.csv.gz'
            else:
                pcoa_file_name = f'{self.genus}_True_True_True_False_biallelic_{self.dist_method}_dist_10000_pwr_False_{self.misco}_{self.masco}_3_pcoa.csv.gz'

            pcoa_path = os.path.join(self.output_dir_18s, pcoa_file_name)
        
        pcoa_df = pd.read_csv(pcoa_path)
        pcoa_df.set_index('sample', drop=True, inplace=True)
        # Get rid of the proportion explained
        pcoa_df = pcoa_df.iloc[:-1,:]
        # Get rid of tech replicates and convert readset to sample-id
        drop_list = []
        for ind in pcoa_df.index:
            if not self.meta_info_df.at[ind, 'is_representative_for_sample']:
                drop_list.append(ind)
        
        # Now drop the rows and columns
        pcoa_df.drop(index=drop_list, inplace=True)

        # Now we need to convert these to sample-id format.
        sample_id_list = []
        for ind in pcoa_df.index:
            sample_id_list.append(self.fastq_info_df.at[ind, 'sample-id'])
        
        pcoa_df.index = sample_id_list
        
        
        return pcoa_df

    def clustering_overview(self):
        """
        Produce a 3 x 6 figure that gives an overview of the clustering for the SNP and 18S.
        Top row is SNP, first three cols = POC second three = POR. First col == no categories,
        second col = low k categoris thrid col = high k categories. Then repeat for POR
        Second row is Brayurtis 18S
        Thrid row is Unifrac 18S
        """
        fig, ax_arr = plt.subplots(3,6, figsize=(36,18))
        # SNP
        poc_snp_kmeans_dict = compress_pickle.load(os.path.join(self.cache_dir, 'poc_biallelic_snp_kmeans_dict.p.bz'))
        por_snp_kmeans_dict = compress_pickle.load(os.path.join(self.cache_dir, 'por_biallelic_snp_kmeans_dict.p.bz'))
        # POC do k at 5 and 10
        ax_arr[0,0].scatter(self.poc_snp_pcoa_df['PC1'], self.poc_snp_pcoa_df['PC2'], c='black')
        ax_arr[0,0].set_title(f'SNP_Pocillopora_biallelic')
        
        # k = 5
        kmeans = poc_snp_kmeans_dict[5]
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[0,1].scatter(self.poc_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], self.poc_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[0,1].set_title(f'SNP_Pocillopora_biallelic k=5')
        
        # k = 10
        kmeans = poc_snp_kmeans_dict[10]
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[0,2].scatter(self.poc_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], self.poc_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[0,2].set_title(f'SNP_Pocillopora_biallelic k=10')
        
        # POR do k at 3 and 7
        ax_arr[0,3].scatter(self.por_snp_pcoa_df['PC1'], self.por_snp_pcoa_df['PC2'], c='black')
        ax_arr[0,3].set_title(f'SNP_Porites_biallelic')
        
        # k = 3
        kmeans = por_snp_kmeans_dict[3]
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[0,4].scatter(self.por_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], self.por_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[0,4].set_title(f'SNP_Porites_biallelic k=3')

        # k = 7
        kmeans = por_snp_kmeans_dict[7]
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[0,5].scatter(self.por_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], self.por_snp_pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[0,5].set_title(f'SNP_Porites_biallelic k=7')

        # 18S
        # Bray
        # POC k at 4 and 7
        pcoa_file_name = f'Pocillopora_True_True_True_False_biallelic_braycurtis_dist_10000_pwr_False_0.08_200_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax_arr[1,0].scatter(pcoa_df['PC1'], pcoa_df['PC2'], c='black')
        ax_arr[1,0].set_title(f'18S_Pocillopora_braycurtis misco=0.08 masco=200')
        #k = 4
        kmeans = KMeans(n_clusters=4, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[1,1].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[1,1].set_title(f'18S_Pocillopora_braycurtis misco=0.08 masco=200 k=4')
        
        # k = 7
        kmeans = KMeans(n_clusters=7, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[1,2].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[1,2].set_title(f'18S_Pocillopora_braycurtis misco=0.08 masco=200 k=7')

        # POR k at 4 and 6 on misco 0.2 masco 80
        pcoa_file_name = f'Porites_True_True_True_False_biallelic_braycurtis_dist_10000_pwr_False_0.2_80_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax_arr[1,3].scatter(pcoa_df['PC1'], pcoa_df['PC2'], c='black')
        ax_arr[1,3].set_title(f'18S_Porites_braycurtis misco=0.20 masco=80')
        
        #k = 4
        kmeans = KMeans(n_clusters=4, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[1,4].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[1,4].set_title(f'18S_Porites_braycurtis misco=0.20 masco=80 k=4')
        
        # k = 6
        kmeans = KMeans(n_clusters=6, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[1,5].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[1,5].set_title(f'18S_Porites_braycurtis misco=0.20 masco=80 k=6')
        
        # Unifrac
        # POC misco at 0.08 masco at 200, k at 3 and 6
        pcoa_file_name = f'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_pwr_False_0.08_200_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax_arr[2,0].scatter(pcoa_df['PC1'], pcoa_df['PC2'], c='black')
        ax_arr[2,0].set_title(f'18S_Poocillopora_unifrac misco=0.08 masco=200')
        
        #k = 3
        kmeans = KMeans(n_clusters=3, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[2,1].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[2,1].set_title(f'18S_Poocillopora_unifrac misco=0.08 masco=200 k=3')
        
        # k = 6
        kmeans = KMeans(n_clusters=6, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[2,2].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[2,2].set_title(f'18S_Poocillopora_unifrac misco=0.08 masco=200 k=6')

        # POR misco at 0.66 masco at 100, k at 3 and 5
        pcoa_file_name = f'Porites_True_True_True_False_biallelic_unifrac_dist_1000_pwr_False_0.66_100_3_pcoa.csv.gz'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, pcoa_file_name))
        ax_arr[2,3].scatter(pcoa_df['PC1'], pcoa_df['PC2'], c='black')
        ax_arr[2,3].set_title(f'18S_Porites_unifrac misco=0.66 masco=100')

        #k = 3
        kmeans = KMeans(n_clusters=3, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[2,4].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[2,4].set_title(f'18S_Porites_unifrac misco=0.66 masco=100 k=3')
        
        # k = 5
        kmeans = KMeans(n_clusters=5, n_init=100, algorithm='full').fit(pcoa_df)
        for i in range(kmeans.cluster_centers_.shape[0]):
            ax_arr[2,5].scatter(pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 0], pcoa_df.iloc[np.where(kmeans.labels_ == i)[0], 1])
        ax_arr[2,5].set_title(f'18S_Porites_unifrac misco=0.66 masco=100 k=5')

        plt.savefig(os.path.join(self.eighteens_dir, 'temp_clustering_overview.png'), dpi=300)

    def produce_r_input(self):
        """
        We will run 6 PERMANOVAs. One for each species, marker, dist method combination.
        We will run these in R but we will output the input files here.
        For each PERMANOVA we will need a pair of files. One distance matrix and one meta info df.
        What details are included in the meta info files will differ between the SNP and 18S analyses.
        For both we will have Island and site. For the 18S we will also have, pre-QC seq depth, pre-QC richness
        post-QC richness. We will need to look back at the qc_output_dictionaries to get these values.
        """

        # SNP POC output
        # dist matrice
        self.poc_snp_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_snp_biallelic.dist.csv'), index=False, header=False)
        # meta df
        meta_df = self.sample_provenance_df.loc[self.poc_snp_df.index]
        meta_df = meta_df[meta_df['SAMPLE PROTOCOL LABEL, level 2'] == 'CS4L']
        meta_df = meta_df.loc[:,['ISLAND#', 'SITE#']]
        meta_df.rename(columns={'ISLAND#':'ISLAND', 'SITE#':'SITE'}, inplace=True)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_snp_biallelic.meta.csv'), index_label='sample_id', index=True, header=True)

        # SNP POR output
        # dist matrice
        self.por_snp_df.to_csv(os.path.join(self.output_dir_18s, 'porites_snp_biallelic.dist.csv'), index=False, header=False)
        # meta df
        meta_df = self.sample_provenance_df.loc[self.por_snp_df.index]
        meta_df = meta_df[meta_df['SAMPLE PROTOCOL LABEL, level 2'] == 'CS4L']
        meta_df = meta_df.loc[:,['ISLAND#', 'SITE#']]
        meta_df.rename(columns={'ISLAND#':'ISLAND', 'SITE#':'SITE'}, inplace=True)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'porites_snp_biallelic.meta.csv'), index_label='sample_id', index=True, header=True)

        # 18S
        # POC BrayCurtis
        dist_file_name = f'Pocillopora_True_True_True_False_biallelic_braycurtis_dist_10000_pwr_False_0.08_200_3.dist.gz'
        dist_df, meta_df = self._get_dist_meta_df(dist_path=os.path.join(self.output_dir_18s, dist_file_name), abundance_dict_file_name=dist_file_name.replace('.dist.gz', '_abundance_dict.p.bz'), genus='Pocillopora')
        dist_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_18S_braycurtis_0.08_200.dist.csv'), index=False, header=False)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_18S_braycurtis_0.08_200.meta.csv'), index=True, header=True, index_label='sample_id')

        # POC Unifrac
        dist_file_name = f'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_pwr_False_0.08_200_3.dist.gz'
        dist_df, meta_df = self._get_dist_meta_df(dist_path=os.path.join(self.output_dir_18s, dist_file_name), abundance_dict_file_name=dist_file_name.replace('.dist.gz', '_abundance_dict.p.bz'), genus='Pocillopora')
        dist_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_18S_unifrac_0.08_200.dist.csv'), index=False, header=False)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'pocillopora_18S_unifrac_0.08_200.meta.csv'), index=True, header=True, index_label='sample_id')

        # POR BrayCurtis
        dist_file_name = f'Porites_True_True_True_False_biallelic_braycurtis_dist_10000_pwr_False_0.2_80_3.dist.gz'
        dist_df, meta_df = self._get_dist_meta_df(dist_path=os.path.join(self.output_dir_18s, dist_file_name), abundance_dict_file_name=dist_file_name.replace('.dist.gz', '_abundance_dict.p.bz'), genus='Pocillopora')
        dist_df.to_csv(os.path.join(self.output_dir_18s, 'porites_18S_braycurtis_0.2_80.dist.csv'), index=False, header=False)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'porites_18S_braycurtis_0.2_80.meta.csv'), index=True, header=True, index_label='sample_id')

        # POR Unifrac
        dist_file_name = f'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_pwr_False_0.66_100_3.dist.gz'
        dist_df, meta_df = self._get_dist_meta_df(dist_path=os.path.join(self.output_dir_18s, dist_file_name), abundance_dict_file_name=dist_file_name.replace('.dist.gz', '_abundance_dict.p.bz'), genus='Pocillopora')
        dist_df.to_csv(os.path.join(self.output_dir_18s, 'porites_18S_unifrac_0.66_100.dist.csv'), index=False, header=False)
        meta_df.to_csv(os.path.join(self.output_dir_18s, 'porites_18S_unifrac_0.66_100.meta.csv'), index=True, header=True, index_label='sample_id')

        foo = 'bar'

    def _get_dist_meta_df(self, dist_path, abundance_dict_file_name, genus):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        
        dist_df = pd.read_csv(dist_path, index_col=0, header=None)

        # Get rid of tech replicates and convert readset to sample-id
        drop_list = []
        for ind in dist_df.index:
            if not self.meta_info_df.at[ind, 'is_representative_for_sample']:
                drop_list.append(ind)
        
        # Now drop the rows and columns
        dist_df.drop(index=drop_list, inplace=True)

        # Here we need to make the meta while we still have the readset index as 
        # qc info is stored with readset as key
        normalised_abundance_dict = compress_pickle.load(os.path.join(self.cache_dir, abundance_dict_file_name))
        meta_df = pd.DataFrame(index=normalised_abundance_dict.keys(), columns=['normalised_richness', 'raw_richness', 'raw_abund'])
        for k, v in normalised_abundance_dict.items():
            meta_df.at[k, 'normalised_richness'] = len(v.keys())
            # Here we can look up the consolidation dict
            sample_qc_dir = os.path.join(self.qc_dir, k)
            consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            absolute_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'abs_all_seq_abundance_dict.p.bz'))
            coral_keys = [k for k in absolute_abund_dict.keys() if k in coral_annotation_dict]
            genus_keys = [k for k in coral_keys if coral_annotation_dict[k] == genus]
            raw_abund = sum([absolute_abund_dict[k] for k in genus_keys])
            
            meta_df.at[k, 'raw_richness'] = len(consolidated_host_seqs_abund_dict.keys())
            meta_df.at[k, 'raw_abund'] = raw_abund
        
        foo = 'bar'

        meta_df = meta_df.reindex(dist_df.index)
        assert(list(meta_df.index) == list(dist_df.index))
        
        # Now we need to convert these to sample-id format.
        sample_id_list = []
        for ind in dist_df.index:
            sample_id_list.append(self.fastq_info_df.at[ind, 'sample-id'])
        dist_df.index = sample_id_list
        meta_df.index = sample_id_list

        # Now add the island and site info to the meta_df
        meta_df['ISLAND'] = self.sample_provenance_df.loc[meta_df.index, 'ISLAND#']
        meta_df['SITE'] = self.sample_provenance_df.loc[meta_df.index, 'SITE#']

        return dist_df, meta_df

class OneClusterCol:
    def __init__(self, genus, dist_method, misco, masco, ax, parent):
        self.parent = parent
        self.genus = genus
        self.dist_method = dist_method
        self.misco = misco
        self.masco = masco
        # Ax is an array of the axes (a single column worths)
        self.ax = ax
        self.figure_path = os.path.join(self.parent.eighteens_dir, 'temp_fig_cluster_{self.genus}_{self.dist_method}_18s.png')
        self.pcoa_df = self._get_pcoa_df()
        self.k_range = range(2,11,1)
        self.inertia = []
        self.sil = {}
        # The KMeans objects for the ks computed for the given genus, dist method and parameter set
        self.kmeans_dict = {}
        self._set_agreement_dicts()
        
        # These will be set dynamically after the kmeans have been calculated
        # and will be recalculated with every kmeans
        self.label_df_18s_data = None
        self.label_df_18s_random = None
        self.label_df_snp_data = None
        self.label_df_snp_random = None
        # The category mapping dicts will also be generated with every kmeans calculation
        self.cat_mapping_dicts = None
        self.agreements_data = None
        self.agreements_random = None
        self.max_k = None
        self.kmeans = None
        
    def _get_pcoa_df(self, pcoa_path=None):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        if pcoa_path is None:
            if self.dist_method == 'unifrac':
                pcoa_file_name = f'{self.genus}_True_True_True_False_biallelic_{self.dist_method}_dist_1000_pwr_False_{self.misco}_{self.masco}_3_pcoa.csv.gz'
            else:
                pcoa_file_name = f'{self.genus}_True_True_True_False_biallelic_{self.dist_method}_dist_10000_pwr_False_{self.misco}_{self.masco}_3_pcoa.csv.gz'

            pcoa_path = os.path.join(self.parent.output_dir_18s, pcoa_file_name)
        
        pcoa_df = pd.read_csv(pcoa_path)
        pcoa_df.set_index('sample', drop=True, inplace=True)
        # Get rid of the proportion explained
        pcoa_df = pcoa_df.iloc[:-1,:]
        # Get rid of tech replicates and convert readset to sample-id
        drop_list = []
        for ind in pcoa_df.index:
            if not self.parent.meta_info_df.at[ind, 'is_representative_for_sample']:
                drop_list.append(ind)
        
        # Now drop the rows and columns
        pcoa_df.drop(index=drop_list, inplace=True)

        # Now we need to convert these to sample-id format.
        sample_id_list = []
        for ind in pcoa_df.index:
            sample_id_list.append(self.parent.fastq_info_df.at[ind, 'sample-id'])
        
        pcoa_df.index = sample_id_list
        
        
        return pcoa_df
    
    def _set_agreement_dicts(self):
        """We have created a cache system for these dictionaries as they take considerable time to compute
        The dictionaries hold for a given k, 
        'data' - the agreement between the 18s and snp based on the data
        'random' - the same but using random labels
        'ratios' - the ratio of the best data agreement and the best random. We did this so that we can judge
        which k worked best taking into account the fact that agreement might be better soley due to the fact
        that k was smaller.
        """
        self.max_agreement_dict_data_pickle_path = os.path.join(self.parent.cache_dir_18s, f'{self.genus}_{self.dist_method}_biallelic_{self.misco}_{self.masco}_max_agreement_dict_data.p.bz')
        self.max_agreement_dict_random_pickle_path = os.path.join(self.parent.cache_dir_18s, f'{self.genus}_{self.dist_method}_biallelic_{self.misco}_{self.masco}_max_agreement_dict_random.p.bz')
        self.max_agreement_dict_ratios_pickle_path = os.path.join(self.parent.cache_dir_18s, f'{self.genus}_{self.dist_method}_biallelic_{self.misco}_{self.masco}_max_agreement_dict_ratios.p.bz')
        if os.path.exists(self.max_agreement_dict_data_pickle_path) and os.path.exists(self.max_agreement_dict_random_pickle_path) and os.path.exists(self.max_agreement_dict_ratios_pickle_path):
            self.max_agreement_dict_data = compress_pickle.load(self.max_agreement_dict_data_pickle_path)
            self.max_agreement_dict_random = compress_pickle.load(self.max_agreement_dict_random_pickle_path)
            self.max_agreement_dict_ratios = compress_pickle.load(self.max_agreement_dict_ratios_pickle_path)
        else:
            self.max_agreement_dict_data = {}
            self.max_agreement_dict_random = {}
            self.max_agreement_dict_ratios = {}
    
    def plot_col(self):
        # first calculate the inertial values and silhoutte scores to produce elbow plot
        # and silhoutte plot on same axes to inform fit of k values.
        self._calculate_inertia_silhoutte_scores()
        self._plot_inertia_silhoutte()
        
        # For calculating the agreements we will temporarily work with a reduced number of k values
        # as checking the agreements using the category mapping dictionaries is very expensive for larger
        # values of k
        self.k_range = range(2,7,1)
        self._compute_18s_snp_agreements()
        self._plot_agreement_values()
        
        # Plot up the first 2 coords coloured according to the clustering with the best agreement
        self._plot_scatter_cluster_agreement()

        # Then plot up coloured according to the snp clustering
        # For this we will only be able to plot up those points
        # that have an snp classficiation as coloured. For the others, plot as grey.
        self._plot_scatter_cluster_snp_cats()

    def _plot_scatter_cluster_snp_cats(self):
        """
        Plot up the 18s first two components but coloured according to the snp categories at the max_k k.
        Where categories are not available (most samples), plot as grey.
        """
        if self.genus == 'Pocillopora':
            label_df_snp = pd.Series(self.parent.poc_snp_kmeans_dict[self.max_k].labels_, index=self.parent.poc_snp_df.index, name='label')
        else:
            label_df_snp = pd.Series(self.parent.por_snp_kmeans_dict[self.max_k].labels_, index=self.parent.por_snp_df.index, name='label')

        colours = []
        # Firstly we want to scatter up the samples that aren't found in the snp
        # as light grey 
        scat_samples = [_ for _ in self.pcoa_df.index if _ not in label_df_snp.index]
        self.ax[3].scatter(self.pcoa_df.loc[scat_samples,'PC1'], self.pcoa_df.loc[scat_samples,'PC2'], s=16, c='lightgrey')
        for i in range(self.kmeans.cluster_centers_.shape[0]): # for each k
            # Get a list of samples that we want to plot up
            # This gives us the rows where the lables are i in the snp
            i_samples = label_df_snp[label_df_snp==i].index
            # Then we want to scatter those samples that are in the pcoa_df
            scat_samples = [_ for _ in self.pcoa_df.index if _ in i_samples]
            scat = self.ax[3].scatter(self.pcoa_df.loc[scat_samples,'PC1'], self.pcoa_df.loc[scat_samples,'PC2'], s=16)
            
            colours.append(scat._original_facecolor[0])
            
        reset_xlim = self.ax[3].get_xlim()
        reset_ylim = self.ax[3].get_ylim()
        
        for i in range(self.kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            # Plot the centroid as vline hline intersect
            #vline
            centroid_x = self.kmeans.cluster_centers_[i][0]
            centorid_y = self.kmeans.cluster_centers_[i][1]
            
            self.ax[3].plot([centroid_x, centroid_x],[self.ax[3].get_ylim()[0], self.ax[3].get_ylim()[1]], c=colours[i])
            self.ax[3].plot([self.ax[3].get_xlim()[0], self.ax[3].get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
        
        # lims need resetting as the plotting of the lines changes them.
        self.ax[3].set_xlim(reset_xlim)
        self.ax[3].set_ylim(reset_ylim)
        self.ax[3].set_title(f'k={self.max_k}')

    def _plot_scatter_cluster_agreement(self):
        """
        Plot up the first two comonents of the 18S pcoa and colour according to the kmeans categores
        given at the max_k.
        """
        self.max_k = sorted(self.max_agreement_dict_ratios, key=self.max_agreement_dict_ratios.get, reverse=True)[0]
        self.kmeans = self.kmeans_dict[self.max_k]
        colours = []
        for i in range(self.kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            scat = self.ax[2].scatter(self.pcoa_df.iloc[np.where(self.kmeans.labels_==i)[0],0], self.pcoa_df.iloc[np.where(self.kmeans.labels_==i)[0],1], s=16)
            colours.append(scat._original_facecolor[0])
            
        reset_xlim = self.ax[2].get_xlim()
        reset_ylim = self.ax[2].get_ylim()
        for i in range(self.kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            # Plot the centroid as vline hline intersect
            #vline
            centroid_x = self.kmeans.cluster_centers_[i][0]
            centorid_y = self.kmeans.cluster_centers_[i][1]
            
            self.ax[2].plot([centroid_x, centroid_x],[self.ax[2].get_ylim()[0], self.ax[2].get_ylim()[1]], c=colours[i])
            self.ax[2].plot([self.ax[2].get_xlim()[0], self.ax[2].get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
        # lims need resetting as the plotting of the lines changes them.
        self.ax[2].set_xlim(reset_xlim)
        self.ax[2].set_ylim(reset_ylim)
        self.ax[2].set_title(f'k={self.max_k}')

    def _plot_agreement_values(self):
        # Here we are ready to plot up the agreement values
        self.ax[1].plot([_ for _ in self.k_range], [self.max_agreement_dict_data[_] for _ in self.k_range], 'k-')
        self.ax[1].set_xticks([_ for _ in self.k_range])
        self.ax[1].set_xlabel('k')
        self.ax[1].set_ylabel('agreement score')
        ax2 = self.ax[1].twinx()
        ax2.plot([_ for _ in self.k_range], [self.max_agreement_dict_ratios[_] for _ in self.k_range], 'r-')
        ax2.set_ylabel('data/random ratio')
        self.ax[1].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

    def _generate_label_dfs(self, k):
        # create label dfs for both the 18s and the snp.
        
        # 18S
        self.label_df_18s_data = pd.Series(self.kmeans_dict[k].labels_, index=self.pcoa_df.index, name='label')
        # Also create a label df where the labels have been randomly shuffled
        # So that we can compare the agreement scores against agreement scores based
        # on chance alone.
        random_labels = random.sample(list(self.label_df_18s_data), len(self.label_df_18s_data))
        self.label_df_18s_random = pd.Series(random_labels, index=self.label_df_18s_data.index, name='label_random')

        # SNP
        if self.genus == 'Pocillopora':
            self.label_df_snp_data = pd.Series(self.parent.poc_snp_kmeans_dict[k].labels_, index=self.parent.poc_snp_df.index, name='label')
        else:
            self.label_df_snp_data = pd.Series(self.parent.por_snp_kmeans_dict[k].labels_, index=self.parent.por_snp_df.index, name='label')
        random_labels = random.sample(list(self.label_df_snp_data), len(self.label_df_snp_data))
        self.label_df_snp_random = pd.Series(random_labels, index=self.label_df_snp_data.index, name='label_random')
        
    def _agreement_data(self, cat_mapping_dict):
        count = 0
        agree = 0
        for sample, label in self.label_df_18s_data.items():
            if sample in self.label_df_snp_data:
                count += 1
                if cat_mapping_dict[self.label_df_18s_data[sample]] == self.label_df_snp_data[sample]:
                    agree += 1
        return agree/count

    def _agreement_random(self, cat_mapping_dict):
        count = 0
        agree = 0
        for sample, label in self.label_df_18s_random.items():
            if sample in self.label_df_snp_random:
                count += 1
                if cat_mapping_dict[self.label_df_18s_random[sample]] == self.label_df_snp_random[sample]:
                    agree += 1
        return agree/count

    def _calculate_agreement_for_given_k(self, k):
        self.agreements_data = []
        self.agreements_random = []
        print(f'Assessing mapped agreements for k = {k}')
        tot = len(self.cat_mapping_dicts)
        # We can actually perhaps calculate the best agreement between the 
        # the mappings by looking at the centroids of both kmeans.
        for cat_i, cat_mapping_dict in enumerate(self.cat_mapping_dicts):
            sys.stdout.write(f'\r{cat_i}/{tot}')
            # DATA
            # Calculate agreement for each of the category mappings
            agreement_data = self._agreement_data(cat_mapping_dict)
            self.agreements_data.append(agreement_data)
            
            agreement_random = self._agreement_random(cat_mapping_dict)
            self.agreements_random.append(agreement_random)
        print(f'\nComplete for k={k}')

    def _compute_18s_snp_agreements(self):
        for k in self.k_range:
            if k in self.max_agreement_dict_ratios:
                continue
            
            # pandas Series of sample to label
            self._generate_label_dfs(k)
            assert(set(self.label_df_18s_data.unique()) == set(self.label_df_snp_data.unique()))
            
            # For each 18S sample, look up the call and see if it agrees
            # We need to take into acount that the label names are likely different between the two clusterings.
            # As such we'll need to work out the best possible agreement.
            # The list of the differnt category mappings can be generated using permutations
            self.cat_mapping_dicts = [{k:v for k, v in zip(self.label_df_18s_data.unique(), _)} for _ in itertools.permutations(self.label_df_18s_data.unique(), len(self.label_df_18s_data.unique()))]
            
            self._calculate_agreement_for_given_k(k)
            
            self._populate_and_pickle_agreement_dicts(k)
            
            print(f'max agreement: {max(self.agreements_data)}')

    def _populate_and_pickle_agreement_dicts(self, k):
        self.max_agreement_dict_data[k] = max(self.agreements_data)
        self.max_agreement_dict_random[k] = max(self.agreements_random)
        self.max_agreement_dict_ratios[k] = max(self.agreements_data)/max(self.agreements_random)
        compress_pickle.dump(self.max_agreement_dict_data, self.max_agreement_dict_data_pickle_path)
        compress_pickle.dump(self.max_agreement_dict_random, self.max_agreement_dict_random_pickle_path)
        compress_pickle.dump(self.max_agreement_dict_ratios, self.max_agreement_dict_ratios_pickle_path)

    def _plot_inertia_silhoutte(self):
        # here we have the inertia ready to be plotted up
        self.ax[0].plot([_ for _ in self.k_range], self.inertia, 'k-')
        self.ax[0].set_xticks([_ for _ in self.k_range])
        self.ax[0].set_title(f'misco: {self.misco}; masco {self.masco}')
        self.ax[0].set_xlabel('k')
        self.ax[0].set_ylabel('inertia')
        ax2 = self.ax[0].twinx()
        ax2.plot([_ for _ in self.k_range], [self.sil[_] for _ in self.k_range], 'r-')
        ax2.set_ylabel('silhouette score')
        self.ax[0].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
    
    def _calculate_inertia_silhoutte_scores(self):
        # Save the kmeans calculations so that we don't have to re-calculate below
        for i in self.k_range:
            kmeans = KMeans(n_clusters=i, n_init=100, algorithm='full').fit(self.pcoa_df)
            self.kmeans_dict[i] = kmeans
            self.inertia.append(kmeans.inertia_)
            self.sil[i] = silhouette_score(self.pcoa_df, kmeans.labels_, metric = 'euclidean')
        

    



if __name__ == "__main__":
    c = Cluster18S()
    # c.visalise_snp_df()
    # c.category_scores()
    # c.check_old_clustering_agreement()
    # c.clustering_overview()
    c.produce_r_input()

# TODO
# Run permanovas of the dist matrices with meta data for site, island, 
# sequencing depth pre normalisation, richness post-normalisation.
# Do this for both the 18S and the SNP
# We will presumably need a distance matrix to pass in and then a meta info df.

# TODO check clustering of between old poc and new poc.
# DONE

# TODO, plot up ordination figures one for poc and one for por
# comparing between the 18S structure and the SNP structure.
# 2 ordinations for the SNP (at low and high k)
# and 4 for the 18S, one for each dist method also at high and low k.
# We can use the in detail cluster figs for the 18S to choose which can be representative.
# DONE

# TODO autocorrelation between different 18S distance methods. for k=3 
