#!/usr/bin/env python3

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
# https://stackoverflow.com/a/33616192/5516420 This gets rid of the deprecation warnings from sklearn
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
from sklearn.cluster import KMeans, AffinityPropagation, MiniBatchKMeans
import subprocess
from subprocess import CalledProcessError
from sklearn.metrics import silhouette_score
import compress_pickle
import itertools
import sys
import random
import json
import hashlib
from functools import partial
import ntpath
from sklearn import metrics
# I added fasta2net to the PYTHONPATH in the .bashrc
from fasta2net import Fasta2Net
from collections import defaultdict

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
    
    def visalise_snp_new_clusters(self):
        # TODO we are here.
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
        if os.path.exists(self.poc_snp_kmeans_dict_pickle_path):
            self.poc_snp_kmeans_dict = compress_pickle.load(self.poc_snp_kmeans_dict_pickle_path)
            k_range = sorted(self.poc_snp_kmeans_dict.keys())
            # print('Calculating kmeans')
            for i in k_range:
                print(f'k={i}')
                kmeans = self.poc_snp_kmeans_dict[i]
                poc_inertia.append(kmeans.inertia_)
                poc_sil[i] = silhouette_score(self.poc_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
        else:
            k_range = range(2, 101)
            for i in k_range:
                print(f'Computing k={i}')
                kmeans = KMeans(n_clusters=i, n_init=100, algorithm='full').fit(poc_np)
                poc_inertia.append(kmeans.inertia_)
                poc_sil[i] = silhouette_score(self.poc_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
                self.poc_snp_kmeans_dict[i] = kmeans
                try:
                    compress_pickle.dump(self.poc_snp_kmeans_dict, self.poc_snp_kmeans_dict_pickle_path)
                except:
                    # remove the last k that made it fail and write out again then break
                    k_to_del = sorted(self.poc_snp_kmeans_dict.keys(), reverse=True)[0]
                    del self.poc_snp_kmeans_dict[k_to_del]
                    poc_inertia = poc_inertia[:-1]
                    del poc_sil[k_to_del]
                    compress_pickle.dump(self.poc_snp_kmeans_dict, self.poc_snp_kmeans_dict_pickle_path)
                    break

        k_range = sorted(self.poc_snp_kmeans_dict.keys())
        self.axar[0,0].plot([_ for _ in k_range], poc_inertia, 'k-')
        self.axar[0,0].set_xticks([_ for _ in k_range])
        self.axar[0,0].set_title('Pocillopora')
        self.axar[0,0].set_xlabel('k')
        self.axar[0,0].set_ylabel('inertia')
        self.axar[0,0].grid(linewidth=0.5, color='lightgrey')
        
        ax2 = self.axar[0,0].twinx()
        ax2.plot([_ for _ in k_range], [poc_sil[_]  for _ in k_range], 'r-')
        ax2.set_ylabel('silhouette score')
        self.axar[0,0].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')


        por_np = self.por_snp_pcoa_df.to_numpy()
        por_inertia = []
        por_sil = {}
        if os.path.exists(self.por_snp_kmeans_dict_pickle_path +'boob'):
            self.por_snp_kmeans_dict = compress_pickle(self.por_snp_kmeans_dict_pickle_path)
            k_range = sorted(self.por_snp_kmeans_dict.keys())
            for i in k_range:
                kmeans = self.por_snp_kmeans_dict[i]
                por_inertia.append(kmeans.inertia_)
                por_sil[i] = silhouette_score(self.por_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
        else:
            k_range = range(2, 101)
            for i in k_range:
                print(f'Computing k={i}')
                kmeans = KMeans(n_clusters=i, n_init=100, algorithm='full').fit(por_np)
                por_inertia.append(kmeans.inertia_)
                por_sil[i] = silhouette_score(self.por_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
                self.por_snp_kmeans_dict[i] = kmeans

                try:
                    compress_pickle.dump(self.por_snp_kmeans_dict, self.por_snp_kmeans_dict_pickle_path)
                except:
                    # remove the last k that made it fail and write out again then break
                    k_to_del = sorted(self.por_snp_kmeans_dict.keys(), reverse=True)[0]
                    del self.por_snp_kmeans_dict[k_to_del]
                    por_inertia = por_inertia[:-1]
                    del por_sil[k_to_del]
                    compress_pickle.dump(self.por_snp_kmeans_dict, self.por_snp_kmeans_dict_pickle_path)
                    break

        k_range = sorted(self.por_snp_kmeans_dict.keys())
        self.axar[0,1].plot([_ for _ in k_range], por_inertia, 'k-')
        self.axar[0,1].set_xticks([_ for _ in k_range])
        self.axar[0,1].set_title('Porites')
        self.axar[0,1].set_xlabel('k')
        self.axar[0,1].set_ylabel('inertia')

        ax2 = self.axar[0,1].twinx()
        ax2.plot([_ for _ in k_range], [por_sil[_] for _ in k_range], 'r-')
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

    def _plot_up_kmeans_with_centroids(self, pcoa_df, kmeans, ax, labels=False):
        """
        A utility function for doing the kmeans plotting.
        It will take in a pcoa_df and a kmeans object
        and plot up the clusters on the ax including the centroids.
        """
        colours = []
        for i in range(kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            x = pcoa_df.iloc[np.where(kmeans.labels_==i)[0],0]
            y = pcoa_df.iloc[np.where(kmeans.labels_==i)[0],1]
            point_labels = y.index
            scat = ax.scatter(x, y, s=16)
            if labels==True:
                for j, lab in enumerate(point_labels):
                    if lab.replace('TARA_', '') in ['CO-0001668', 'CO-0002291']:
                        ax.annotate(lab.replace('TARA_', ''), (x[j], y[j]))
            colours.append(scat._original_facecolor[0])
            
        reset_xlim = ax.get_xlim()
        reset_ylim = ax.get_ylim()
        for i in range(kmeans.cluster_centers_.shape[0]): # for each k
            # plot up the centroid and the points in the same colour
            # Plot the centroid as vline hline intersect
            #vline
            centroid_x = kmeans.cluster_centers_[i][0]
            centorid_y = kmeans.cluster_centers_[i][1]
            
            ax.plot([centroid_x, centroid_x],[ax.get_ylim()[0], ax.get_ylim()[1]], c=colours[i])
            ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
        # lims need resetting as the plotting of the lines changes them.
        ax.set_xlim(reset_xlim)
        ax.set_ylim(reset_ylim)

    def check_original_three_agreement(self):
        """
        We have done a distance calculation with the string:
        'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_rai_False_0.08_200_3_original_three'
        Here we are only using islands I06, I10, I15. And we will compare these to the original
        snmf classifications.
        The idea of this figure was really to get our head around what could be going wrong and to do
        some sanity checking.
        1 - plot up the kmeans clustering of the 18S
        2 - plot up the kmeans clustering of the 'biallelic'
        3 - calculate and visualise the agreements between the 18S ordination/classification and the
        SNP-based classifications (both 'bialliclic' and 'old' snp.)
        4 - Check on the agreement between the 'biallelic' and 'old' (snmf) classifications
        5 - We also started to visualise the k means comparisons but it is not worth our time to finish this.
        We never fully finished this figure as we leant what we needed to from it
        and part way through making it we got a full set of classifications from Didier for both species
        using snmf and SVD methods.
        As part of making this figure we have come up with a relatively efficient algorithm for testing the agreement
        between cluster assignments and we will certainly be using this moving forwards.
        For now, we need to move onto optimising the 3 island and the 10+1 island towards the clustering.
        It could be that we should move away from the distance based optimisation.
        """
        fig, ax = plt.subplots(5,3, figsize=(12,24))
        title_font_size = 'x-small'
        pcoa_df = self._get_pcoa_df(pcoa_path=os.path.join(self.output_dir_18s, 'Pocillopora_True_True_True_False_biallelic_unifrac_dist_1000_rai_False_0.08_200_3_original_three_pcoa.csv.gz'))
        # Reordering
        # First plot up the 18S and show its classificaitons
        # Reorder - 1
        # Run kmeans and show what the clustering looks like for the 18S data
        # We should do this will all samples and with only the snp samples
        inertia = []
        silhouette_list = []
        k_range = range(2, len(pcoa_df.index)-1, 1)
        kmeans_dict_18s_cache_path = os.path.join(self.cache_dir_18s, '18s_kmeans_dict_original_three_islands.p.bz')
        if os.path.exists(kmeans_dict_18s_cache_path):
            kmeans_dict_18s = compress_pickle.load(kmeans_dict_18s_cache_path)
            for k in k_range:
                kmeans = kmeans_dict_18s[k]
                s_score = silhouette_score(pcoa_df, kmeans.labels_, metric = 'euclidean')
                silhouette_list.append(s_score)
                inertia.append(kmeans.inertia_)
        else:
            kmeans_dict_18s = {}
            print('Calculating Kmeans for 18s')
            for k in k_range:
                print(f'k={k}')
                kmeans = KMeans(n_clusters=k, n_init=100, algorithm='full').fit(pcoa_df)
                kmeans_dict_18s[k] = kmeans
                s_score = silhouette_score(pcoa_df, kmeans.labels_, metric = 'euclidean')
                silhouette_list.append(s_score)
                inertia.append(kmeans.inertia_)
            compress_pickle.dump(kmeans_dict_18s, kmeans_dict_18s_cache_path)


        # Here plot up the inertia and silloute score
        curr_ax = ax[0,0]
        curr_ax.plot(k_range, inertia)
        ax2 = curr_ax.twinx()
        ax2.plot(k_range, silhouette_list, 'r-')
        curr_ax.set_xlabel('k')
        curr_ax.set_ylabel('inertia')
        ax2.set_ylabel('silhouette score')
        curr_ax.spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        curr_ax.grid(linewidth=0.5, color='lightgrey')

        # Here plot up the inertia and silloute score
        # Now plot up for upto 10
        
        curr_ax = ax[0,1]
        curr_ax.plot(k_range, inertia)
        ax2 = curr_ax.twinx()
        ax2.plot(k_range, silhouette_list, 'r-')
        curr_ax.set_xlabel('k')
        curr_ax.set_ylabel('inertia')
        ax2.set_ylabel('silhouette score')
        curr_ax.spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        curr_ax.grid(linewidth=0.5, color='lightgrey')
        curr_ax.set_xlim(0,10)

        # plot k=4 as this looks like best candidate
        self._plot_up_kmeans_with_centroids(pcoa_df=pcoa_df, kmeans=kmeans_dict_18s[4], ax=ax[0,2], labels=True)
        
        
        # Reorder-2 Plot up the kmeans for the new snp
        poc_np = self.poc_snp_pcoa_df.to_numpy()
        poc_inertia = []
        poc_sil = {}
        # Hold the kmeans results in a dict and pickle out this dict so that we can use it
        # in comparing between the 18S and SNP categorisation agreements
        if os.path.exists(self.poc_snp_kmeans_dict_pickle_path):
            self.poc_snp_kmeans_dict = compress_pickle.load(self.poc_snp_kmeans_dict_pickle_path)
            k_range = sorted(self.poc_snp_kmeans_dict.keys())
            # print('Calculating kmeans')
            for i in k_range:
                # print(f'k={k}')
                kmeans = self.poc_snp_kmeans_dict[i]
                poc_inertia.append(kmeans.inertia_)
                poc_sil[i] = silhouette_score(self.poc_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
        else:
            k_range = range(2, 101)
            for i in k_range:
                print(f'Computing k={i}')
                kmeans = KMeans(n_clusters=i, n_init=100, algorithm='full').fit(poc_np)
                poc_inertia.append(kmeans.inertia_)
                poc_sil[i] = silhouette_score(self.poc_snp_pcoa_df, kmeans.labels_, metric = 'euclidean')
                self.poc_snp_kmeans_dict[i] = kmeans
                try:
                    compress_pickle.dump(self.poc_snp_kmeans_dict, self.poc_snp_kmeans_dict_pickle_path)
                except:
                    # remove the last k that made it fail and write out again then break
                    k_to_del = sorted(self.poc_snp_kmeans_dict.keys(), reverse=True)[0]
                    del self.poc_snp_kmeans_dict[k_to_del]
                    poc_inertia = poc_inertia[:-1]
                    del poc_sil[k_to_del]
                    compress_pickle.dump(self.poc_snp_kmeans_dict, self.poc_snp_kmeans_dict_pickle_path)
                    break

        k_range = sorted(self.poc_snp_kmeans_dict.keys())
        curr_ax = ax[1,0]
        curr_ax.plot([_ for _ in k_range], poc_inertia, 'k-')
        curr_ax.set_xticks([_ for _ in k_range])
        curr_ax.set_title('Pocillopora')
        curr_ax.set_xlabel('k')
        curr_ax.set_ylabel('inertia')
        curr_ax.grid(linewidth=0.5, color='lightgrey')
        
        ax2 = curr_ax.twinx()
        ax2.plot([_ for _ in k_range], [poc_sil[_]  for _ in k_range], 'r-')
        ax2.set_ylabel('silhouette score')
        curr_ax = ax[3,2].spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

        curr_ax = ax[1,1]
        curr_ax.plot([_ for _ in k_range], poc_inertia, 'k-')
        curr_ax.set_xticks([_ for _ in k_range])
        curr_ax.set_title("Kmeans clustering of the SNP 'biallelic' distance matrix")
        curr_ax.set_xlabel('k')
        curr_ax.set_ylabel('inertia')
        curr_ax.grid(linewidth=0.5, color='lightgrey')
        ax2 = curr_ax.twinx()
        ax2.plot([_ for _ in k_range], [poc_sil[_]  for _ in k_range], 'r-')
        ax2.set_ylabel('silhouette score')
        curr_ax.spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        curr_ax.set_xlim(0,25)

        # Finally, plot up the new snp dists with clusters coloured at a given k.
        self._plot_up_kmeans_with_centroids(pcoa_df=self.poc_snp_pcoa_df, kmeans=self.poc_snp_kmeans_dict[4], labels=True, ax=ax[1,2])
        
        # Now put the k=3 snmf classificaiton, the k=3 and the k=4 classification (new SNP) onto the 18S
        #### Reorder-3 - Put the snmf classifications onto this at k=3. 
        path_to_classifications = '/home/humebc/projects/tara/tara_full_dataset_processing/18s/input/snp_dist_matrices/pocillopora_snp_clusters.csv'
        c_df = pd.read_csv(path_to_classifications)
        c_df.drop(columns='Colony', inplace=True)
        c_df.rename(columns={'code pangea':'sample-id', 'snmf/PCA':'label'}, inplace=True)
        c_df.set_index('sample-id', drop=True, inplace=True)
        c_df['label'] = c_df['label'] -1
        # annotate the two that look like they should be a fourth.
        # We can find these out by looking at the coordinates
        c_df['label_4'] = list(c_df['label'])
        c_df.at['TARA_CO-0002135', 'label_4'] = 3
        c_df.at['TARA_CO-0002065', 'label_4'] = 3
        
        
        non_classified_samples = [_ for _ in pcoa_df.index if _ not in c_df.index]
        classified_samples = [_ for _ in pcoa_df.index if _ in c_df.index]
        ax[2,0].scatter(pcoa_df.loc[non_classified_samples, 'PC1'], pcoa_df.loc[non_classified_samples, 'PC2'], c='lightgrey')
        for i in c_df['label'].unique():
            x = pcoa_df.loc[c_df[c_df['label'] == i].index, 'PC1']
            y = pcoa_df.loc[c_df[c_df['label'] == i].index, 'PC2']
            lab = pcoa_df.loc[c_df[c_df['label'] == i].index, 'PC1'].index
            samples = ax[2,0].scatter(x, y)
            for j, lab in enumerate(lab):
                if lab.replace('TARA_', '') in ['CO-0001668', 'CO-0002291']:
                    ax[2,0].annotate(lab.replace('TARA_', ''), (x[j], y[j]))

        ax[2,0].set_title('location=18S colour=snmf categories (SNP) k=3\nAgreement with k=3 and k=3 = 22/27 = 0.81', fontsize=title_font_size)



            
        # new SNP onto 18s at k==3
        poc_biallelic_kmeans = compress_pickle.load(self.poc_snp_kmeans_dict_pickle_path)[3]
        poc_snp_labels_series = pd.Series(poc_biallelic_kmeans.labels_, index=self.poc_snp_pcoa_df.index, name='label')
        
        poc_snp_labels_series = poc_snp_labels_series.loc[[_ for _ in poc_snp_labels_series.index if _ in pcoa_df.index]]
        ax[2,1].scatter(pcoa_df.loc[[_ for _ in pcoa_df.index if _ not in poc_snp_labels_series], 'PC1'], pcoa_df.loc[[_ for _ in pcoa_df.index if _ not in poc_snp_labels_series], 'PC2'], c='lightgrey')
        for i in range(poc_biallelic_kmeans.cluster_centers_.shape[0]):
            x = pcoa_df.loc[poc_snp_labels_series[poc_snp_labels_series == i].index, 'PC1']
            y = pcoa_df.loc[poc_snp_labels_series[poc_snp_labels_series == i].index, 'PC2']
            labels = y.index
            samples = ax[2,1].scatter(x,y)
            for j, lab in enumerate(labels):
                if lab.replace('TARA_', '') in ['CO-0001668', 'CO-0002291']:
                    ax[2,1].annotate(lab.replace('TARA_', ''), (x[j], y[j]))
        ax[2,1].set_title(f'location=18S colour=biallelic (new SNP) k=3\nAgreement with k=3 and k=3 = 15/27 = 0.56', fontsize=title_font_size)

        # as a very simple proof of principle we can plot up 18S three islands with the
        # k=4 agreement (new SNP)
        poc_biallelic_kmeans = compress_pickle.load(self.poc_snp_kmeans_dict_pickle_path)[4]
        poc_snp_labels_series = pd.Series(poc_biallelic_kmeans.labels_, index=self.poc_snp_pcoa_df.index, name='label')
        
        poc_snp_labels_series = poc_snp_labels_series.loc[[_ for _ in poc_snp_labels_series.index if _ in pcoa_df.index]]
        ax[2,2].scatter(pcoa_df.loc[[_ for _ in pcoa_df.index if _ not in poc_snp_labels_series], 'PC1'], pcoa_df.loc[[_ for _ in pcoa_df.index if _ not in poc_snp_labels_series], 'PC2'], c='lightgrey')
        for i in range(poc_biallelic_kmeans.cluster_centers_.shape[0]):
            x = pcoa_df.loc[poc_snp_labels_series[poc_snp_labels_series == i].index, 'PC1']
            y = pcoa_df.loc[poc_snp_labels_series[poc_snp_labels_series == i].index, 'PC2']
            labels = y.index
            samples = ax[2,2].scatter(x,y)
            for j, lab in enumerate(labels):
                if lab.replace('TARA_', '') in ['CO-0001668', 'CO-0002291']:
                    ax[2,2].annotate(lab.replace('TARA_', ''), (x[j], y[j]))
        ax[2,2].set_title('location=18S colour=biallelic (new SNP) k=4\nAgreement with k=3 (18S) and k=4 (SNP) = 20/27 = 0.74', fontsize=title_font_size)

        
        # Work out which the two samples are that are looking like they should be a fourth by looking at the coordinates
        interest_df = self.poc_snp_pcoa_df[(self.poc_snp_pcoa_df['PC1'] < 0) & (self.poc_snp_pcoa_df['PC2'] < 200)]
        interest_df = interest_df.loc[[_ for _ in interest_df.index if _ in c_df.index]]
        print('The c_df points that are probably a fourth category are:')
        for _ in interest_df.index:
            print(f'\t{_}')
        
        #TODO we are here.
        # Reorder-4 plot up the new kmeans distances but with the old snmf classifications.
        # First plot up the k=3 for the biallelic. (k=4 is already plotted above)
        self._plot_up_kmeans_with_centroids(pcoa_df=self.poc_snp_pcoa_df, kmeans=self.poc_snp_kmeans_dict[3], labels=True, ax=ax[3,0])
        ax[3,0].set_title('Clustering of the biallelic SNP distance data at k=3\n(see above for k=4)', fontsize=title_font_size)
        
        # This shows us where the disagreement is happening and that there is a problem
        # projecting big onto small. But not vice versa.
        # The proof is to plot up the bad agreement on the second plot that shows
        # the bad agreement between old and new.
        # But then to show good agreement between a higher number of k on the new snp
        # but with 4 or 5 on the old.
        ax[3,1].scatter(self.poc_snp_pcoa_df.loc[[_ for _ in self.poc_snp_pcoa_df.index if _ not in c_df.index], 'PC1'], self.poc_snp_pcoa_df.loc[[_ for _ in self.poc_snp_pcoa_df.index if _ not in c_df.index], 'PC2'], c='lightgrey')
        for i in c_df['label'].unique():
            x = self.poc_snp_pcoa_df.loc[c_df[c_df['label'] == i].index, 'PC1']
            y = self.poc_snp_pcoa_df.loc[c_df[c_df['label'] == i].index, 'PC2']
            ax[3,1].scatter(x, y)
        ax[3,1].set_title('location=SNP biallelic colour=snmf k=3\nAgreement with k=3 (biallelic) and k=3 (snmf) = 18/27 = 0.67', fontsize=title_font_size)

        ax[3,2].scatter(self.poc_snp_pcoa_df.loc[[_ for _ in self.poc_snp_pcoa_df.index if _ not in c_df.index], 'PC1'], self.poc_snp_pcoa_df.loc[[_ for _ in self.poc_snp_pcoa_df.index if _ not in c_df.index], 'PC2'], c='lightgrey')
        for i in c_df['label_4'].unique():
            x = self.poc_snp_pcoa_df.loc[c_df[c_df['label_4'] == i].index, 'PC1']
            y = self.poc_snp_pcoa_df.loc[c_df[c_df['label_4'] == i].index, 'PC2']
            ax[3,2].scatter(x, y)
        ax[3,2].set_title('location=SNP biallelic colour=snmf k=4\nAgreement with k=4 (biallelic) and k=4 (snmf) = 27/27 = 1.00', fontsize=title_font_size)

        
        # Calculate the agreement here
        # We want to calculate the agreement for (for the 18s/snmf):
        # 3-->3
        # 3-->4
        # 4-->3
        # 4-->4
        # cat_mapping_dicts = [{k:v for k, v in zip(np.unique(kmeans.labels_), _)} for _ in itertools.permutations(np.unique(kmeans.labels_), len(np.unique(kmeans.labels_)))]
        # This is where we're going to have to thinkk about the logic. cat_mapping_dicts works when you're doing the comparisons
        # between the same number of categories (although we need to make this more efficient).
        # But when you're doing comparisons with different numbers of categories then things get pretty tricky.
        # We could iter our way through the smaller category list and see which of the cat it finds highest agreement with
        # and then remove this from the list and continue to the next etc. and find the highest agreement this way.
        # To be thorough about this I guess we could do the small list in the order of all permuations.
        # Perhaps we can see how time intensive this is. It is probably a good idea to work out how we can generalise this
        # so that we can pull it out into a method. I guess, we'll just want a label_df. i.e. a mapping of sample-id
        # to predicted category. So to map the three to the three then we ould pass in two 3 kmeans, and to 
        # map the 3 to the 4 we would pass in the 3 or 4 respective to which one we wanted the 4 four.
        
        # First do the agreement between the snmf and 18s original three islands.
        
        three_three = AgreementCalculator(lab_ser_one=pd.Series(kmeans_dict_18s[3].labels_, index=pcoa_df.index, name='label'), lab_ser_two=c_df['label'], parent=self).calculate_agreement()
        # self._test_kmeans_agreement(
        #     pd.Series(kmeans_dict_18s[3].labels_, index=pcoa_df.index, name='label'), 
        #     c_df['label'])
        three_four = AgreementCalculator(
            pd.Series(kmeans_dict_18s[3].labels_, index=pcoa_df.index, name='label'), 
            c_df['label_4'], parent=self).calculate_agreement()
        four_three = AgreementCalculator(
            pd.Series(kmeans_dict_18s[4].labels_, index=pcoa_df.index, name='label'), 
            c_df['label'], parent=self).calculate_agreement()
        four_four = AgreementCalculator(
            pd.Series(kmeans_dict_18s[4].labels_, index=pcoa_df.index, name='label'), 
            c_df['label_4'], parent=self).calculate_agreement()
        print('\n\nAgreement scores for snmf --> 18S original three')
        print(f'k=3 --> k=3: {three_three}')
        print(f'k=3 --> k=4: {three_four}')
        print(f'k=4 --> k=3: {four_three}')
        print(f'k=3 --> k=3: {four_four}')
        snmf_18s_agreement = np.empty((2,2))
        snmf_18s_agreement[0,0] = three_three
        snmf_18s_agreement[0,1] = three_four
        snmf_18s_agreement[1,0] = four_three
        snmf_18s_agreement[1,1] = four_four

        # Now we want to try doing the comparisons of the 18S and the snmf to the new snp
        print('\n\nAgreement scores for snmf SNP --> biallelic SNP ')
        # TODO here do two sets of nested for loops to do agreement comparisons
        # Agreement of snmf with biallelic
        snmf_agreement = np.empty((2,28))
        for s_i, snmf_index in enumerate(['label', 'label_4']):
            for k_i, k in enumerate(range(2,30)):
                snmf_agreement[s_i,k_i] = AgreementCalculator(c_df[snmf_index], pd.Series(self.poc_snp_kmeans_dict[k].labels_, index=self.poc_snp_df.index, name='label'), parent=self).calculate_agreement()
        
        # Agreement of 18S with biallelic
        print('\n\nAgreement scores for 18S SNP --> biallelic SNP ')
        e_agreement = np.empty((4,28))
        # TODO I think it is taking some time to make the series so create look up tables for these
        print('Creating look up dicts of label series for 18S. This may take some time...')
        e_series_lookup = {e_k: pd.Series(kmeans_dict_18s[e_k].labels_, index=pcoa_df.index, name='label') for e_k in range(2,30)}
        print('Creating look up dicts of label series for biallelic SNP. This may take some time...')
        b_series_lookup = {b_k: pd.Series(self.poc_snp_kmeans_dict[b_k].labels_, index=self.poc_snp_df.index, name='label') for b_k in range(2,30)}

        for e_i, e_k in enumerate(range(2,6)):
            for b_i, b_k in enumerate(range(2,30)):
                snmf_agreement[s_i,k_i] = AgreementCalculator(
                    e_series_lookup[e_k],
                    b_series_lookup[b_k], parent=self).calculate_agreement()
        
        #TODO we are here.
        # Plot up three heat maps to visualise the k agreements for the three blocks above.
        im = ax[4,0].imshow(snmf_18s_agreement)
        cbar = ax[4,0].figure.colorbar(im, ax=ax[4,0])
        cbar.ax.set_ylabel('agreement', rotation=-90, va="bottom")
        ax[4,0].set_xticks(np.arange(2))
        ax[4,0].set_yticks(np.arange(2))
        ax[4,0].set_xticklabels([3,4])
        ax[4,0].set_yticklabels([3,4])

        ax[4,1].imshow(snmf_agreement)
        ax[4,2].imshow(e_agreement)

        plt.tight_layout()
        plt.savefig(os.path.join(self.fig_output_dir_18s, 'snmf_agreement.png'), dpi=600)

class CheckPwrRai(EighteenSBase):
    """
    We are seeing some very strange disagreement between the PWR and the RAI normalisation in unifrac
    It seems like this was something going wrong in the ms_figure and agreement calculation scripting.
    I'm walking away from that now.
    We want to plot a network of the sequences. We will work with the misco and masco values used
    here to generate a network that we will attempt to colour according to the the snp classifications.
    This way we can hopefully see those sequnces that can be used to inform the classification
    and those that cannot.
    """
    def __init__(self, misco='0.01', masco='17'):
        super().__init__()
        self.misco = misco
        self.masco = masco
        self.meta_info_table_path = '/home/humebc/projects/tara/tara_full_dataset_processing/output/coral_18S_meta_info_table_2020-04-15T12_53_51.189071UTC.csv'
        self.meta_info_df = self._get_meta_info_df()
        
        self.genus = 'Pocillopora'
        self.distance_method = 'unifrac'
        self.pcoa_path_pwr = os.path.join(self.output_dir_18s, f'{self.genus}_True_True_True_False_biallelic_{self.distance_method}_dist_1000_pwr_False_{misco}_{masco}_3_original_three_pcoa.csv.gz')
        self.pcoa_path_rai = os.path.join(self.output_dir_18s, f'{self.genus}_True_True_True_False_biallelic_{self.distance_method}_dist_1000_rai_False_{misco}_{masco}_3_original_three_pcoa.csv.gz')
        self.pcoa_df_pwr = self._get_pcoa_df(self.pcoa_path_pwr)
        # self.pwr_var = self.pcoa_df_pwr.iloc[-1]
        self.pcoa_df_rai = self._get_pcoa_df(self.pcoa_path_rai)
        # self.pwr_var = self.pcoa_df_rai.iloc[-1]
        # self.pcoa_df_pwr = self.pcoa_df_pwr.iloc[:-1]
        # self.pcoa_df_rai = self.pcoa_df_rai.iloc[:-1]
        self.fig, self.ax_arr = plt.subplots(2,4, figsize=(16,8))
        self.fig_dir = '/home/humebc/projects/tara/tara_full_dataset_processing/18s/figures'
        self.n_clusters = 3
        self.fig_path = os.path.join(self.fig_dir, f'pwr_rai_{self.misco}_{self.masco}_{self.n_clusters}.png')
        self.pwr_kmeans_labels = pd.Series(KMeans(n_clusters=self.n_clusters, n_init=100).fit(self.pcoa_df_pwr).labels_, index=self.pcoa_df_pwr.index, name='label')
        self.rai_kmeans_labels = pd.Series(KMeans(n_clusters=self.n_clusters, n_init=100).fit(self.pcoa_df_rai).labels_, index=self.pcoa_df_rai.index, name='label')
        self.snp_classification = self.poc_snp_classifications
        self.agreement_list = []
        # To make the networks we will need to work with the abundance df
        # The abundance df can be got from the cache.
        # We will work with the pwr pcoa as this seems to be getting better agreement (slightly) than the rai
        self.abundance_dict_path = os.path.join(self.cache_dir, ntpath.basename(self.pcoa_path_pwr.replace('_pcoa.csv.gz', '_abundance_dict.p.bz')))
        # To keep things simple when checking the agreement results it will be helpful if the index
        # is the same for both the abundance_df and the pcoa_df so that the classification labels
        # are in the same order.
        self.abundance_df, self.abundance_dict = self._curate_abundance_df(self.abundance_dict_path)
        self.abundance_df = self.abundance_df.reindex(self.pcoa_df_pwr.index)

    def make_network(self):
        """
        Make a network from the sequences contained in the abundance_df. These are the sequences that the ordinations
        are being made with. We will colour them according to the snp classifications. For each sequence, we will
        get the cummulative abundance of the sqeuence and we will split this abundnace up by the different categories
        I am excited for this. It will be a new way to look at the data.
        Hopefully we will be able to see specific areas of the network that are specific to certain
        classifications. We will hopefully also be able to see areas that are not specific. We'll then in theory
        be able to imporove the ordination by selecting for the good sequences and throwing out the bad.
        We will hopefully be able to identify a pattern for doing this and possibly we will be able to use something
        like MED.
        """
        # First let's check that we can get a black and white network plotted up using the fasta2net.
        # This may not be so simple.
        # We have rewritten Fasta2Net so that it will take an abundance_dict as input
        # Currently our abundance dictionary is a dictionary of dictionaries where first key is the readset name
        # and second key is nucleotide sequence. We want to pass in a single abundance dictionary to Fasta2Net (tostart)
        # where it key is sequence and value is cummulative relative abundance
        net_abund_dict, c_dist_dict = self._make_fasta2net_abundance_and_color_dicts()
        net = Fasta2Net(abundance_dictionary=net_abund_dict, output_dir=os.path.join(self.output_dir_18s, 'network'), size_scaler=0.1, color_distribution_dict=c_dist_dict)
        net.make_network()
        foo = 'bar'
    
    def _make_fasta2net_abundance_and_color_dicts(self):
        """
        Create an abundance dictionary to be used in the network generation
        key = nuc seuqnece, value = value (cumrelabund * 1000)
        Also generate a colour disribution dictionary for use in the network
        key = nuc sequence, value is a list of the distribution of the sequence
        across as many categories as we want. Given that the categories we want
        to work with here are the SNP classifications the number of classifications
        will be the number of classification that cover the samples we are dealing with plus
        one for those samples that done have an snp classification
        """
        classification_index_dict = {classification:ind for ind, classification in enumerate(self.snp_classification['label'].unique())}
        classification_index_dict['none'] = len(classification_index_dict.keys())
        net_abund_dict = defaultdict(float)
        classification_distribution_dict = defaultdict(self.class_ddict_constructor)
        for k, v in self.abundance_dict.items(): # For each sample
            try:
                snp_classification = self.snp_classification.at[k, 'label']
            except KeyError:
                snp_classification = 'none'
            classification_index = classification_index_dict[snp_classification]
            tot = sum(v.values())
            for k_, v_ in v.items(): # sequence, absolute abundance pairing
                rel_abund = v_/tot
                net_abund_dict[k_] += rel_abund
                classification_distribution_dict[k_][classification_index] += rel_abund
        return {k: int(v*1000) for k, v in net_abund_dict.items()}, classification_distribution_dict

    def class_ddict_constructor(self):
        return [0 for i in range( len(self.snp_classification['label'].unique()) + 1 )]

    def plot(self):
        # Start by visualising the first 3 components of the PCOA's and colouring according to the kmeans clustering
        # at k==3
        # Plot up the pwr row
        self._plot_row(
            labels = self.pwr_kmeans_labels, pcoa_df = self.pcoa_df_pwr,
            pc2_ax_kmeans_cats = self.ax_arr[0,0], pc3_ax_kmeans_cats = self.ax_arr[0,1], 
            pc2_ax_snp_cats = self.ax_arr[0,2], pc3_ax_snp_cats = self.ax_arr[0,3], norm_m='pwr'
            )

        # Plot up the rai row
        self._plot_row(
            labels = self.rai_kmeans_labels, pcoa_df = self.pcoa_df_rai,
            pc2_ax_kmeans_cats = self.ax_arr[1,0], pc3_ax_kmeans_cats = self.ax_arr[1,1], 
            pc2_ax_snp_cats = self.ax_arr[1,2], pc3_ax_snp_cats = self.ax_arr[1,3], norm_m='rai'
            )

        plt.tight_layout()
        
        # Append the agreement results to the name
        self.fig_path = self.fig_path.replace('.png', f'_{self.agreement_list[0]:.2f}_{self.agreement_list[1]:.2f}.png')
        plt.savefig(self.fig_path, dpi=600)
        

    def _curate_abundance_df(self, abundance_dict_path):
        """
        Drop tech reps and convert readsets to sample-id
        """
        # Read in the abundance df
        abundance_dict = compress_pickle.load(abundance_dict_path)
        abundance_df = pd.DataFrame.from_dict(abundance_dict, orient='index')
        abundance_df[pd.isna(abundance_df)] = 0
        # Get rid of tech replicates and convert readset to sample-id
        drop_list = []
        for ind in abundance_df.index:
            if not self.meta_info_df.at[ind, 'is_representative_for_sample']:
                drop_list.append(ind)
        
        # Now drop the rows and columns
        abundance_df.drop(index=drop_list, inplace=True)

        # Drop in the dict too
        for k in drop_list:
            del abundance_dict['key']

        # Now we need to convert these to sample-id format.
        sample_id_list = []
        for ind in abundance_df.index:
            sample_id_list.append(self.fastq_info_df.at[ind, 'sample-id'])
        abundance_df.index = sample_id_list

        # Convert dict too
        abundance_dict = {self.fastq_info_df.at[k, 'sample-id'] : v for k, v in abundance_dict.items()}
        
        return abundance_df, abundance_dict

    def _get_meta_info_df(self):
        return pd.read_csv(self.meta_info_table_path).set_index('readset')

    def _plot_row(self, labels, pcoa_df, pc2_ax_kmeans_cats, pc3_ax_kmeans_cats, pc2_ax_snp_cats, pc3_ax_snp_cats, norm_m):
        for k in labels.unique():
            samples_pwr = labels[labels==k].index
            pc2_ax_kmeans_cats.scatter(pcoa_df.loc[samples_pwr, 'PC1'], pcoa_df.loc[samples_pwr, 'PC2'])
            pc3_ax_kmeans_cats.scatter(pcoa_df.loc[samples_pwr, 'PC1'], pcoa_df.loc[samples_pwr, 'PC3'])
        
        samples_in_common = list(set(self.snp_classification.index).intersection(set(labels.index)))
        labels = labels.reindex(index=samples_in_common)
        snp_classification_reindexed = self.snp_classification.reindex(index=samples_in_common)
        
        agreement = AgreementCalculator(lab_ser_one=labels, lab_ser_two=snp_classification_reindexed['label'], temp_dir_18s=self.temp_dir_18s, cache_dir_18s=self.cache_dir_18s).calculate_agreement()
        self.agreement_list.append(agreement)
        pc2_ax_kmeans_cats.set_title(f'{norm_m}, PC1 PC2\nmisco {self.misco} masco {self.masco}\nagreement={agreement:.2f}')
        pc3_ax_kmeans_cats.set_title(f'{norm_m}, PC1 PC3\nmisco {self.misco} masco {self.masco}\nagreement={agreement:.2f}')

        # Secondly, plot up in the next two axes the 18s ordiantion first 3 coordinates but with the snp classiciations
        samples_not_in_common = [_ for _ in pcoa_df.index if _ not in self.snp_classification.index]
        pc2_ax_snp_cats.scatter(pcoa_df.loc[samples_not_in_common, 'PC1'], pcoa_df.loc[samples_not_in_common, 'PC2'], c='lightgrey')
        pc3_ax_snp_cats.scatter(pcoa_df.loc[samples_not_in_common, 'PC1'], pcoa_df.loc[samples_not_in_common, 'PC3'], c='lightgrey')
        for svd in self.snp_classification['label'].unique():
            of_svd = self.snp_classification[self.snp_classification['label']==svd]
            samples_in_common = [_ for _ in pcoa_df.index if _ in of_svd.index]
            pc2_ax_snp_cats.scatter(pcoa_df.loc[samples_in_common, 'PC1'], pcoa_df.loc[samples_in_common, 'PC2'])
            pc3_ax_snp_cats.scatter(pcoa_df.loc[samples_in_common, 'PC1'], pcoa_df.loc[samples_in_common, 'PC3'])

        pc2_ax_snp_cats.set_title(f'{norm_m}, PC1 PC2\nmisco {self.misco} masco {self.masco}\nagreement={agreement:.2f}')
        pc3_ax_snp_cats.set_title(f'{norm_m}, PC1 PC3\nmisco {self.misco} masco {self.masco}\nagreement={agreement:.2f}')

    
    def _get_pcoa_df(self, pcoa_path):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        try:
            pcoa_df = pd.read_csv(pcoa_path, index_col=0)
        except Exception as e:
            foo = 'bar'
        # try:
        #     pcoa_df.set_index('sample', drop=True, inplace=True)
        # except KeyError:
        #     pass
        # Get rid of the proportion explained
        try:
            pcoa_df = pcoa_df.iloc[:-1,:]
        except Exception as e:
            foo = 'bar'
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


class AgreementCalculator:
    """
        Takes in a pairs of series where the index is sample name and the value is the classification.
        Performs a test to see the maximum possible agreement between the classficiations
        To do this it will work with the smallest set of labels, (either the '*_one' or '*_two')
        It will take a list of the classifications and work in order of those classifications
        For each classification, it will find the classification in the other dataset that
        corresponds to the largest agreement. It will then move on to the next classification
        untill it has found the best assignments for that list of classifications. Because order
        will matter, we will then do this for every permuation of the small classifications.

        Doing this for every permutation become inhibitive when we get to higher k values for the
        smaller series. For example k=8 is 40K permutations to test and the likelihood of finding a better
        result is very small. So let's implement a rule whereby we test 10 of the permutations randomly
        and then we test another 100 and if the best match doesn't change then we don't test any further.
    """
    def __init__(self, lab_ser_one, lab_ser_two, parent=None, cache_dir_18s=None, temp_dir_18s=None, silence=False):
        # We allow the input of either the parent object or the individual cache_dir_18s, temp_dir_18s
        if parent is None:
            assert((cache_dir_18s is not None) and (temp_dir_18s is not None))
            self.cache_dir_18s = cache_dir_18s
            self.temp_dir_18s = temp_dir_18s
        else:
            self.cache_dir_18s = parent.cache_dir_18s
            self.temp_dir_18s = parent.temp_dir_18s
        if len(lab_ser_one.unique()) < len(lab_ser_two.unique()):
            self.small_series = lab_ser_one
            self.large_series = lab_ser_two
        else:
            self.small_series = lab_ser_two
            self.large_series = lab_ser_one
        self.silence = silence
        try:
            # Need to attempt to convert to an int as int32 cannot be jsoned.
            self.small_unique_classifications = [int(_) for _ in self.small_series.unique()]
        except ValueError:
            self.small_unique_classifications = [_ for _ in self.small_series.unique()]
        self.large_unique_classifications = list(self.large_series.unique())
        # Calculation of the permutation lists takes considerable time.
        # We could work with the generator but ideally we would like to do the shuffle.
        # rather we will cache out and see if that speeds up.
        if len(self.small_unique_classifications) > 5:
            list_hash = self._md5sum_from_python_object(self.small_unique_classifications)
            perm_list_cache_path = os.path.join(self.cache_dir_18s, f'{list_hash}_perm_list.p.bz')
            if os.path.exists(perm_list_cache_path):
                self.ordered_clasification_lists = compress_pickle.load(perm_list_cache_path)
            else:
                if not self.silence:
                    print('Creating permutation list from scratch. This may take some time...')
                self.ordered_clasification_lists = list(itertools.permutations(self.small_unique_classifications, len(self.small_unique_classifications)))
                random.shuffle(self.ordered_clasification_lists)
                compress_pickle.dump(self.ordered_clasification_lists, perm_list_cache_path)
                if not self.silence:
                    print('Complete.')
        else:
            # Don't bother with the cache
            self.ordered_clasification_lists = list(itertools.permutations(self.small_unique_classifications, len(self.small_unique_classifications)))
            random.shuffle(self.ordered_clasification_lists)
        
        # We can make some look up dicts to speed this up. One for each series
        # Key will be the classification, value will be set of the sample_ids
        self.small_look_up = {}
        for classification in self.small_series.unique():
            self.small_look_up[classification] = set(self.small_series[self.small_series == classification].index)

        self.large_look_up = {}
        for classification in self.large_series.unique():
            self.large_look_up[classification] = set(self.large_series[self.large_series == classification].index)

        self.max_absolute_agreement = 0
        self.agreement_score_list = []
        self.samples_in_common = len(set(self.small_series.index).intersection(set(self.large_series.index)))

        self.done_slice_score = 0

    def calculate_agreement(self):
        if not self.silence:
            print('Running agreement calculations')
        self.tot_perm_lists = len(self.ordered_clasification_lists)
        if len(self.ordered_clasification_lists) > 200:
            # Do for 10, then 100. if the score doesn't change 
            self._agreement_slice(10)
            current_max = self.max_absolute_agreement
            # This means that we will continue to run in 100 slices until
            # we stop finding higher max agreements
            count = 0
            while True:
                self._agreement_slice(100)
                # Exit out if we didn't find a bigger score
                count += 1
                if current_max == self.max_absolute_agreement:
                    break
                else:
                    current_max = self.max_absolute_agreement
                # Exit out if we have done all of the permutations
                if self.done_slice_score > self.tot_perm_lists:
                    break
                if count == 4:
                    foo = 'bar'
        else:
            self._agreement_slice()
        if not self.silence:
            print('\rAll permutations complete')
            print(f'Max agreement k={len(self.small_unique_classifications)} and k={len(self.large_unique_classifications)}: {self.max_absolute_agreement}/{self.samples_in_common} = {max(self.agreement_score_list):.2f}')
        return max(self.agreement_score_list)

    def _agreement_slice(self, slice_size=None):
        if slice_size is not None:
            slice_list = self.ordered_clasification_lists[self.done_slice_score:self.done_slice_score + slice_size]
            self.done_slice_score += slice_size
        else:
            slice_list = self.ordered_clasification_lists
        for i, ordered_clasification_list in enumerate(slice_list):
            if not self.silence:
                sys.stdout.write(f'\r{i}/{self.tot_perm_lists} permutations complete')
            # For each calssification in order
            matched_classifications = []
            matched_score = 0
            for classification in ordered_clasification_list:
                # Search through the classifications of the other series to find the best agreement
                # For the time being if serveral matches have the same agreeemnt score,
                # we will go with the one found first.
                max_classification = None
                max_score = 0
                for other_classification in [_ for _ in self.large_unique_classifications if _ not in matched_classifications]:
                    intersection_score = len(self.small_look_up[classification].intersection(self.large_look_up[other_classification]))
                    if intersection_score > max_score:
                        max_classification = other_classification
                        max_score = intersection_score
                if max_classification is not None:
                    # Then add the union_score to the matched_score
                    matched_score += max_score
                    matched_classifications.append(max_classification)
                else:
                    pass
                # If there is no match we just don't put anything in the matched_clasifications
            # The agreement score should be given as a proportion of the maximum possible score.
            # The maximum possible score is of course the number of samples in common between the
            # two datasets.
            if matched_score > self.max_absolute_agreement:
                self.max_absolute_agreement = matched_score
            self.agreement_score_list.append(matched_score/self.samples_in_common)

    def _md5sum_from_python_object(self, p_obj):
        """ 
        Get a hash code for a given python object
        https://stackoverflow.com/questions/5417949/computing-an-md5-hash-of-a-data-structure
        """
        return hashlib.md5(json.dumps(sorted(p_obj), sort_keys=True).encode('utf-8')).hexdigest()

class ComputeClassificationAgreement:
    """
    Class for computing the classificaiton agreement between 18S and the SNP classifications
    Alongside all the other objects, the method takes in either a distance matrix or a pcoa path.
    If the distance method is 'braycurtis' or 'jaccard' then the method expects a distance matrix
    to be provided. In this case, the samples will be filtered according to the island list and a
    pcoa will be computed and the classification agreement will be conducted from this.
    If distance_method is 'unifrac' or 'PCA' then a pcoa_path is expected as the island_list-specific
    pcoa/pca has already been computed. In this case we will simply readin the pcoa/pca from this 
    path.
    """
    def __init__(self, distance_method, island_list, genus, output_dir, cache_dir, cache_dir_18s, input_dir_18s, output_dir_18s, fastq_info_df_path, temp_dir_18s, distance_matrix_path=None, pcoa_path=None):
        self.distance_method = distance_method
        self.cache_dir = cache_dir
        self.island_list = island_list
        self.output_dir = output_dir
        self.cache_dir_18s = cache_dir_18s
        self.input_dir_18s = input_dir_18s
        self.output_dir_18s = output_dir_18s
        self.fastq_info_df = self._make_fastq_info_df(fastq_info_df_path=fastq_info_df_path)
        self.temp_dir_18s = temp_dir_18s
        self.genus = genus
        self.snp_class_df = self._load_snp_classification_df()
        if self.distance_method in ['braycurtis', 'jaccard']:
            assert(pcoa_path is None and distance_matrix_path is not None)
            self.distance_matrix_path = distance_matrix_path
            self.pcoa_path = self.distance_matrix_path.replace('.dist.gz', f'_{self.island_list}_pcoa.csv.gz')
            self.abundance_dict_path = os.path.join(self.cache_dir, ntpath.basename(self.distance_matrix_path.replace('.dist.gz', '_abundance_dict.p.bz'))) 
        elif self.distance_method in ['unifrac', 'PCA']:
            assert(pcoa_path is not None and distance_matrix_path is None)
            self.pcoa_path = pcoa_path
            if self.distance_method == 'unifrac':
                # A distance matrix is available for the 'unifrac' but not for PCA obviously.
                # We can use this for calculating some hierarchical clustering.
                # TODO actually turn this into a distance matrix df below.
                self.distance_matrix_path = pcoa_path.replace('_pcoa.csv.gz', '.dist.gz')
                self.abundance_dict_path = os.path.join(self.cache_dir, ntpath.basename(self.pcoa_path.replace('_pcoa.csv.gz', '_abundance_dict.p.bz'))) 
                self.results_path = self.pcoa_path.replace('_pcoa.csv.gz', '_classification_result.txt')
                self.kmeans_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pcoa.csv.gz', '_kmeans_dict.p.bz'))
                self.kmeans_raw_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pcoa.csv.gz', '_kmeans_raw_dict.p.bz'))
                self.agreement_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pcoa.csv.gz', '_agreement_dict.p.bz'))
            else:
                self.abundance_dict_path = os.path.join(self.cache_dir, ntpath.basename(self.pcoa_path.replace('_pca.csv.gz', '_abundance_dict.p.bz'))) 
                self.results_path = self.pcoa_path.replace('_pca.csv.gz', '_classification_result.txt')
                self.kmeans_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pca.csv.gz', '_kmeans_dict.p.bz'))
                self.kmeans_raw_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pca.csv.gz', '_kmeans_raw_dict.p.bz'))
                self.agreement_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.pcoa_path).replace('_pca.csv.gz', '_agreement_dict.p.bz'))
        else:
            raise NotImplementedError
        
        # Use the meta_info_table produced by output_tables.py
        self.meta_info_table_path = os.path.join(self.output_dir, 'coral_18S_meta_info_table_2020-04-15T12_53_51.189071UTC.csv')
        self.meta_info_df = self._get_meta_info_df()
        # If braycurtis, then we need to remove the additional samples from the matrix
        # and then compute a pcoa.
        # If unifrac then the pcoa should already have been calculated and we can read it in
        if self.distance_method in ['braycurtis', 'jaccard']:
            # Then we need to strip out the samples that do not belong to the islands
            # specified in the island_list
            try:
                dist_df = pd.read_csv(self.distance_matrix_path, index_col=0, header=None).astype(float)
                dist_df.columns = dist_df.index
            except ValueError:
                dist_df = pd.read_csv(self.distance_matrix_path, index_col=0).astype(float)
            # Create the abundance df from the abundance dictionary so that we can compute 
            # clusters based on the raw count data.
            # Same as the distance matrix we will need to drop the samples that aren't found
            self.abundance_df = pd.DataFrame.from_dict(compress_pickle.load(self.abundance_dict_path), orient='index')
            self.abundance_df[pd.isna(self.abundance_df)] = 0
            
            drop_list = []
            if self.island_list == 'original_three':
                island_list = ['I06', 'I10', 'I15']
            elif self.island_list == 'ten_plus_one':
                island_list = ['I01','I02','I03','I04','I05', 'I06','I07','I08','I09','I10', 'I15']
            else:
                raise NotImplementedError
            samples_from_islands = self.meta_info_df[self.meta_info_df['ISLAND#'].isin(island_list)].index
            to_drop = [_ for _ in dist_df.index if _ not in samples_from_islands]
            dist_df.drop(index=to_drop, columns=to_drop, inplace=True)
            self.abundance_df.drop(index=to_drop, inplace=True)
            # Now we have the dist_df in shape and we should calculate the pcoa
            print('executing pcoa')
            pcoa_result = pcoa(dist_df)
            print('completed pcoa')
            self.pcoa_df = pcoa_result.samples
            self.pcoa_df.index = dist_df.index
            self.pcoa_df.index.name = 'sample'
            # now add the variance explained as a final row to the renamed_dataframe
            self.pcoa_df = self.pcoa_df.append(pcoa_result.proportion_explained.rename('proportion_explained'))
            # At this point we should write out the pcoa_df and then use the 
            # get_pcoa method to read it in with sample-id format names
            # by writing this out we can use it again when we come to do the other clustering methods
            
            self.pcoa_df.to_csv(self.pcoa_path)
            self.pcoa_df = self._get_pcoa_df(pcoa_path=self.pcoa_path)
            
            # To keep things simple when checking the agreement results it will be helpful if the index
            # is the same for both the abundance_df and the pcoa_df so that the classification labels
            # are in the same order.
            self.abundance_df = self._curate_abundance_df(self.abundance_df)
            self.abundance_df = self.abundance_df.reindex(self.pcoa_df.index)
            
            # Finally the results path should
            self.results_path = self.distance_matrix_path.replace('.dist.gz', f'_{self.island_list}_classification_result.txt')
            self.kmeans_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.distance_matrix_path).replace('.dist.gz', f'_{self.island_list}_kmeans_dict.p.bz'))
            self.kmeans_raw_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.distance_matrix_path).replace('.dist.gz', f'_{self.island_list}_kmeans_raw_dict.p.bz'))
            self.agreement_dict_pickle_path = os.path.join(self.cache_dir_18s, ntpath.basename(self.distance_matrix_path).replace('.dist.gz', f'_{self.island_list}_agreement_dict.p.bz'))

        elif self.distance_method in ['unifrac', 'PCA']:
            self.pcoa_df = self._get_pcoa_df(pcoa_path=self.pcoa_path)
            
            # Also get the abundance dataframe
            # In theory the samples should already have been removed that aren't in the island list
            self.abundance_df = pd.DataFrame.from_dict(compress_pickle.load(self.abundance_dict_path), orient='index')
            self.abundance_df[pd.isna(self.abundance_df)] = 0
            # To keep things simple when checking the agreement results it will be helpful if the index
            # is the same for both the abundance_df and the pcoa_df so that the classification labels
            # are in the same order.
            self.abundance_df = self._curate_abundance_df(self.abundance_df)
            self.abundance_df = self.abundance_df.reindex(self.pcoa_df.index)
        else:
            raise NotImplementedError
        # Dict for holding the kmeans conducted on the pcoa
        self.kmeans_dict = {}
        # Dict for holding the kmeans conducted on the raw sequence counts
        self.kmeans_raw_dict = {}
        self.agreement_dict = {}
        
        # The snp classifications are now associated to the base_18s.py class.
    
    def _load_snp_classification_df(self):
        if self.genus == 'Pocillopora':
            # First check to see if the cached version exists
            snp_cache_dir = os.path.join(self.input_dir_18s, 'snp_classifications', f'poc_snp_class_df.p.bz')
        elif self.genus == 'Porites':
            snp_cache_dir = os.path.join(self.input_dir_18s, 'snp_classifications', f'por_snp_class_df.p.bz')
        if os.path.exists(snp_cache_dir):
            return compress_pickle.load(snp_cache_dir)
        else:
            raise RuntimeError('SNP classification df not found')

    @staticmethod
    def _make_fastq_info_df(fastq_info_df_path):
        df = pd.read_csv(fastq_info_df_path)
        return df.set_index(keys='readset', drop=False)

    def _get_meta_info_df(self):
        return pd.read_csv(self.meta_info_table_path).set_index('readset')

    def _curate_abundance_df(self, abundance_df):
        """
        Drop tech reps and convert readsets to sample-id
        """

        # Get rid of tech replicates and convert readset to sample-id
        drop_list = []
        for ind in abundance_df.index:
            if not self.meta_info_df.at[ind, 'is_representative_for_sample']:
                drop_list.append(ind)
        
        # Now drop the rows and columns
        abundance_df.drop(index=drop_list, inplace=True)

        # Now we need to convert these to sample-id format.
        sample_id_list = []
        for ind in abundance_df.index:
            sample_id_list.append(self.fastq_info_df.at[ind, 'sample-id'])
        
        abundance_df.index = sample_id_list
        return abundance_df

    def _get_pcoa_df(self, pcoa_path):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        try:
            pcoa_df = pd.read_csv(pcoa_path, index_col=0)
        except Exception as e:
            foo = 'bar'
        # try:
        #     pcoa_df.set_index('sample', drop=True, inplace=True)
        # except KeyError:
        #     pass
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
    

    def compute_classficiation_agreement(self, k_range=range(2,20)):
        """
        For each of the .dist matrixes we produced we want to compute several clustering mechansims
        on them and work out the best agreement with the classifications that Didier gave.
        It will likely be good to work with a number of different clustering algorithms
        and we will likely want to multiprocess this given how many comparisons we will have to do.

        For the time being let's just work with a single clustering mechansims and then we can easily implement
        multiple at a later stage.

        We want to calculate agreements for the 'original_three' and the 'ten_plus_one'.
        For the Braycurtis this is as simple as working with a starting distance matrix and removing
        the corresponding samples then clustering.
        For the Unifrac we have to use the precise starting distance matrix that has the island_list
        in the name of the distance matrix.

        We need to work out the framework within which to do these comparisons.
        Perhaps we can do them based on each of the tests we want to run
        """
        
        # We should work with a k that is not less than the clusters assigned according to the SNP
        # but we should allow for the k to be higher for the 18S as there may be additional lineages
        # sampled that were not sampled for the SNP. Given that the complexity of calculating
        # the agreement between two classification levels is goverened by the smaller of the levels
        # and that the smaller level will generally be the k of the SNP data, then we can work with
        # quite a large k. let's start with k=20 to start with
        
        # A dataframe that will hold the results of the clustering from the vairous methods
        classification_df = pd.DataFrame(index=self.abundance_df.index)
        # KMEANS Clustering on PCoA
        # Let's implement a cache for this so that this effort is not wasted
        if os.path.exists(self.kmeans_dict_pickle_path + 'poo'):
            # Get cached kmeans dict
            self.kmeans_dict = compress_pickle.load(self.kmeans_dict_pickle_path)
            for k in k_range:
                kmeans = self.kmeans_dict[k]
                classification_df[f'kmeans_pcoa_{k}'] = kmeans.labels_
        else:
            # Compute kmeans from scratch and pickle out
            for k in k_range:
                kmeans = KMeans(n_clusters=k, n_init=10, algorithm='full').fit(self.pcoa_df)
                self.kmeans_dict[k] = kmeans
                classification_df[f'kmeans_pcoa_{k}'] = kmeans.labels_
            compress_pickle.dump(self.kmeans_dict, self.kmeans_dict_pickle_path)
        
        # # KMEANS Clustering on abundance df
        # # Only do this if we are working with 'unifrac' as it will just be doubling work to do
        # # it with braycurtis too a we will be doing exactly the same calculations
        # # if self.distance_method == 'unifrac':
        # if os.path.exists(self.kmeans_raw_dict_pickle_path):
        #     # Get cached kmeans dict
        #     try:
        #         self.kmeans_raw_dict = compress_pickle.load(self.kmeans_raw_dict_pickle_path)
        #     except EOFError:
        #         for k in k_range:
        #             kmeans = KMeans(n_clusters=k, n_init=10, algorithm='full', n_jobs=1).fit(self.abundance_df)
        #             self.kmeans_raw_dict[k] = kmeans
        #             classification_df[f'kmeans_abund_{k}'] = kmeans.labels_
        #         try:
        #             compress_pickle.dump(self.kmeans_raw_dict, self.kmeans_raw_dict_pickle_path)
        #         except TypeError as e:
        #             pass
        #     for k in k_range:
        #         kmeans = self.kmeans_raw_dict[k]
        #         classification_df[f'kmeans_abund_{k}'] = kmeans.labels_
        # else:
        #     # Compute kmeans from scratch and pickle out
        #     for k in k_range:
        #         kmeans = KMeans(n_clusters=k, n_init=10, algorithm='full').fit(self.abundance_df)
        #         self.kmeans_raw_dict[k] = kmeans
        #         classification_df[f'kmeans_abund_{k}'] = kmeans.labels_
        #     try:
        #         compress_pickle.dump(self.kmeans_raw_dict, self.kmeans_raw_dict_pickle_path)
        #     except TypeError:
        #         # Not sure what is causing this so we just won't dump it out if we have an issue
        #         pass

        # # AffinityPropagation
        # ap_pcoa = AffinityPropagation().fit(self.pcoa_df)
        # classification_df[f'affinitiyp_pcoa'] = ap_pcoa.labels_
        # # if self.distance_method == 'unifrac':
        # ap_abund = AffinityPropagation().fit(self.abundance_df)
        # classification_df[f'affinitiyp_abund'] = ap_pcoa.labels_
        
        # We only want to work with the samples that are found in both snp and 18S samples
        samples_in_common = list(set(classification_df.index).intersection(set(self.snp_class_df.index)))
        self.snp_class_df = self.snp_class_df.reindex(index=samples_in_common)
        classification_df = classification_df.reindex(index=samples_in_common)
        
        max_agreement = 0
        # max_agreement_old = 0
        max_method = None
        for classification_set in list(classification_df):
            # The ARI is great but its also useful to know exactly what percent of agreement you're getting
            # agreement = metrics.adjusted_rand_score(self.snp_class_df['label'], classification_df[classification_set])
            agreement = AgreementCalculator(lab_ser_one=pd.Series(self.kmeans_dict[k].labels_, index=self.pcoa_df.index, name='label'), lab_ser_two=self.snp_class_df['label'], cache_dir_18s=self.cache_dir_18s, temp_dir_18s=self.temp_dir_18s, silence=True).calculate_agreement()
            if agreement > max_agreement:
                max_agreement = agreement
                max_method = classification_set
                # max_agreement_old = agreement_old
            # self.agreement_dict[k] = agreement
        # self.agreement_dict['max'] = max(self.agreement_dict.values())
        # compress_pickle.dump(self.agreement_dict, self.agreement_dict_pickle_path)

        # For the time being we will write out the max agreement the cluster method used to obtain it and the 
        # number of classifications
        with open(self.results_path, 'w') as f:
            f.write(f'{max_agreement:.2f},{max_method}')

        # Here we are done.
        return max_agreement, max_method

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
    # c = Cluster18S()
    # c.visalise_snp_new_clusters()
    # c.category_scores()
    # c.check_old_clustering_agreement()
    # c.clustering_overview()
    # c.produce_r_input()
    # This produces the only partially complete snmf_agreement.png plot
    # that we used as a sanity check to see that we could still get good
    # classification agreement at a smaller number of islands.
    # c.check_original_three_agreement()
    c = CheckPwrRai()
    c.make_network()

    
