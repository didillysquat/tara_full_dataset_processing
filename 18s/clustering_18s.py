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
    
    def _get_pcoa_df(self, genus, dist_method, misco_v, masco):
        # Read in the pcoa file of interest as a df
        # get rid of tech reps and convert readset names to sample-id
        pcoa_file_name = f'{genus}_True_True_True_False_biallelic_{dist_method}_dist_10000_pwr_False_{misco_v}_{masco}_3_pcoa.csv.gz'
        pcoa_path = os.path.join(self.output_dir_18s, pcoa_file_name)
        try:
            subprocess.run(['gzip', '-d', pcoa_path], check=True)
        except CalledProcessError:
            # The file may already be unzipped
            pass
        pcoa_df = pd.read_csv(pcoa_path.replace('.gz', ''))
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
        
        subprocess.run(['gzip', pcoa_path.replace('.gz', '')], check=True)
        return pcoa_df

    def _calculate_inertia_silhoutte_scores(self, pcoa_df, k_range):
        inertia = []
        sil = {}
        # Save the kmeans calculations so that we don't have to re-calculate below
        kmeans_dict = {}
        for i in k_range:
            kmeans = KMeans(n_clusters=i, n_init=100, algorithm='full').fit(pcoa_df)
            kmeans_dict[i] = kmeans
            inertia.append(kmeans.inertia_)
            sil[i] = silhouette_score(pcoa_df, kmeans.labels_, metric = 'euclidean')
        return inertia, sil, kmeans_dict

    def _plot_inertia_silhoutte(self, ax, k_range, inertia, sil, misco_v, masco):
        # here we have the inertia ready to be plotted up
        ax.plot([_ for _ in k_range], inertia, 'k-')
        ax.set_xticks([_ for _ in k_range])
        ax.set_title(f'misco: {misco_v}; masco {masco}')
        ax.set_xlabel('k')
        ax.set_ylabel('inertia')
        ax2 = ax.twinx()
        ax2.plot([_ for _ in k_range], [sil[_] for _ in k_range], 'r-')
        ax2.set_ylabel('silhouette score')
        ax.spines['right'].set_color('red')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

    def _create_agreement_dict(self, genus, dist_method, misco_v, masco):
        """We have created a cache system for these dictionaries as they take considerable time to compute
        The dictionaries hold for a given k, 
        'data' - the agreement between the 18s and snp based on the data
        'random' - the same but using random labels
        'ratios' - the ratio of the best data agreement and the best random. We did this so that we can judge
        which k worked best taking into account the fact that agreement might be better soley due to the fact
        that k was smaller.
        """
        max_agreement_dict_data_pickle_path = os.path.join(self.cache_dir_18s, f'{genus}_{dist_method}_biallelic_{misco_v}_{masco}_max_agreement_dict_data.p.bz')
        max_agreement_dict_random_pickle_path = os.path.join(self.cache_dir_18s, f'{genus}_{dist_method}_biallelic_{misco_v}_{masco}_max_agreement_dict_random.p.bz')
        max_agreement_dict_ratios_pickle_path = os.path.join(self.cache_dir_18s, f'{genus}_{dist_method}_biallelic_{misco_v}_{masco}_max_agreement_dict_ratios.p.bz')
        if os.path.exists(max_agreement_dict_data_pickle_path) and os.path.exists(max_agreement_dict_random_pickle_path) and os.path.exists(max_agreement_dict_ratios_pickle_path):
            max_agreement_dict_data = compress_pickle.load(max_agreement_dict_data_pickle_path)
            max_agreement_dict_random = compress_pickle.load(max_agreement_dict_random_pickle_path)
            max_agreement_dict_ratios = compress_pickle.load(max_agreement_dict_ratios_pickle_path)
        else:
            max_agreement_dict_data = {}
            max_agreement_dict_random = {}
            max_agreement_dict_ratios = {}
        return max_agreement_dict_data, max_agreement_dict_random, max_agreement_dict_ratios

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
        
        genus = 'Porites'
        dist_method = 'braycurtis'
        misco_range = [f'{num:.1f}' for num in np.arange(0.1, 0.6, 0.1)]
        masco = 80
        figure, ax = plt.subplots(4, len(misco_range), figsize=(24,12))
        
        # Make a figure set for each of the misco values
        for misco_i, misco_v in enumerate(misco_range):
            pcoa_df = self._get_pcoa_df(genus, dist_method, misco_v, masco)

            # Use the df to calculate k-means.
            k_range = range(2,11,1)
            inertia, sil, kmeans_dict = self._calculate_inertia_silhoutte_scores(pcoa_df=pcoa_df, k_range=range(2,11,1))
            
            self._plot_inertia_silhoutte(ax=ax[0,misco_i], k_range=k_range, inertia=inertia, sil=sil, misco_v=misco_v, masco=masco)
            
            # Now plot up a line plot underneath the elbow and silhoutte plot that gives
            # the percentage agreement between the 18S and SNP categorsations.
            max_agreement_dict_data, max_agreement_dict_random, max_agreement_dict_ratios = self._create_agreement_dict(genus, dist_method, misco_v, masco)
            
            # for k in range(2,11,1):
            k_range = range(2,7,1)
            for k in k_range:
                if k in max_agreement_dict_ratios:
                    continue
                # We need to clean up the list of 18S samples we are working with first to get rid 
                # of non-representative sequences
                # and then to convert the readset index to a sample name index. Then finally to get rid of
                # and samples that we don't have snp information for.
                
                # First get rid of any samples that are non-representative tech reps
                # To do this, work through the present readsets and identify those that are not representative
                label_df_18s = pd.Series(kmeans_dict[k].labels_, index=pcoa_df.index, name='label')
                
                # create a dict of sample_name to label
                if genus == 'Pocillopora':
                    label_df_snp = pd.Series(self.poc_snp_kmeans_dict[k].labels_, index=self.poc_snp_df.index, name='label')
                else:
                    label_df_snp = pd.Series(self.por_snp_kmeans_dict[k].labels_, index=self.por_snp_df.index, name='label')
                
                # For each 18S sample, look up the call and see if it agrees
                # We need to take into acount that the label names are likely different between the two clusterings.
                # As such we'll need to work out the best possible agreement.
                # The list of the differnt category mappings can be generated using permutations
                assert(set(label_df_18s.unique()) == set(label_df_snp.unique()))
                cat_mapping_dicts = [{k:v for k, v in zip(label_df_18s.unique(), _)} for _ in itertools.permutations(label_df_18s.unique(), len(label_df_18s.unique()))]
                agreements_data = []
                # As a sanity check let's also do this for a random number of points to see what
                # random agreement looks like
                agreements_random = []
                print(f'Assessing mapped agreements for k = {k}')
                tot = len(cat_mapping_dicts)
                # We can actually perhaps calculate the best agreement between the 
                # the mappings by looking at the centroids of both kmeans.
                for cat_i, cat_mapping_dict in enumerate(cat_mapping_dicts):
                    sys.stdout.write(f'\r{cat_i}/{tot}')
                    # DATA
                    # Calculate agreement for each of the category mappings
                    # And work with the largest
                    count = 0
                    agree = 0
                    for sample, label in label_df_18s.items():
                        if sample in label_df_snp:
                            count += 1
                            if cat_mapping_dict[label_df_18s[sample]] == label_df_snp[sample]:
                                agree += 1
                    
                    agreements_data.append(agree/count)
                    
                    label_df_18s_random = pd.Series(
                        random.sample(
                            list(label_df_18s), len(label_df_18s)
                            ), index=label_df_18s.index, name='label_random')
                    label_df_snp_random = pd.Series(
                        random.sample(
                            list(label_df_snp), len(label_df_snp)
                            ), index=label_df_snp.index, name='label_random')
                    
                    count = 0
                    agree = 0
                    for sample, label in label_df_18s_random.items():
                        if sample in label_df_snp_random:
                            count += 1
                            if cat_mapping_dict[label_df_18s_random[sample]] == label_df_snp_random[sample]:
                                agree += 1
                    agreements_random.append(agree/count)
                print(f'\nComplete for k={k}')
                foo = 'bar'
                max_agreement_dict_data[k] = max(agreements_data)
                max_agreement_dict_random[k] = max(agreements_random)
                max_agreement_dict_ratios[k] = max(agreements_data)/max(agreements_random)
                compress_pickle.dump(max_agreement_dict_data, max_agreement_dict_data_pickle_path)
                compress_pickle.dump(max_agreement_dict_random, max_agreement_dict_random_pickle_path)
                compress_pickle.dump(max_agreement_dict_ratios, max_agreement_dict_ratios_pickle_path)
                print(f'max agreement: {max(agreements_data)}')
                #TODO we're here. Let's colour up according to best agreement and according to best silhoutte
                # because it looks like the agreement is pretty bad in general.
                # Alternative is to colour up according to the labels of the snp assignment.
                # I like that! Also remember to plot using only the SNP in the cluster calculation.
            
            # Here we are ready to plot up the agreement values
            ax[1,misco_i].plot([_ for _ in k_range], [max_agreement_dict_data[_] for _ in k_range], 'k-')
            ax[1,misco_i].set_xticks([_ for _ in k_range])
            # ax[1,misco_i].set_title(f'misco: {misco_v}; masco {masco}')
            ax[1,misco_i].set_xlabel('k')
            ax[1,misco_i].set_ylabel('agreement score')
            ax2 = ax[1,misco_i].twinx()
            ax2.plot([_ for _ in k_range], [max_agreement_dict_ratios[_] for _ in k_range], 'r-')
            ax2.set_ylabel('data/random ratio')
            ax[1,misco_i].spines['right'].set_color('red')
            ax2.spines['right'].set_color('red')
            ax2.yaxis.label.set_color('red')
            ax2.tick_params(axis='y', colors='red')

            # Finally plot up the first 2 coords coloured according to the clustering with the best agreement
            max_k = sorted(max_agreement_dict_ratios, key=max_agreement_dict_ratios.get, reverse=True)[0]
            
            # For each of the two k values with the highest sillouette scores
            # plot up the clustering coloured according to the 18s clustering
            
            kmeans = kmeans_dict[max_k]
            colours = []
            for i in range(kmeans.cluster_centers_.shape[0]): # for each k
                # plot up the centroid and the points in the same colour
                scat = ax[2, misco_i].scatter(pcoa_df.iloc[np.where(kmeans.labels_==i)[0],0], pcoa_df.iloc[np.where(kmeans.labels_==i)[0],1], s=16)
                colours.append(scat._original_facecolor[0])
                
            reset_xlim = ax[2,misco_i].get_xlim()
            reset_ylim = ax[2,misco_i].get_ylim()
            for i in range(kmeans.cluster_centers_.shape[0]): # for each k
                # plot up the centroid and the points in the same colour
                # Plot the centroid as vline hline intersect
                #vline
                centroid_x = kmeans.cluster_centers_[i][0]
                centorid_y = kmeans.cluster_centers_[i][1]
                
                ax[2,misco_i].plot([centroid_x, centroid_x],[ax[2,misco_i].get_ylim()[0], ax[2,misco_i].get_ylim()[1]], c=colours[i])
                ax[2,misco_i].plot([ax[2,misco_i].get_xlim()[0], ax[2,misco_i].get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
            # lims need resetting as the plotting of the lines changes them.
            ax[2,misco_i].set_xlim(reset_xlim)
            ax[2,misco_i].set_ylim(reset_ylim)
            ax[2,misco_i].set_title(f'k={max_k}')
            # ax[2,misco_i].scatter(pcoa_df['PC1'], pcoa_df['PC2'])

            # Then plot up coloured according to the snp clustering
            # For this we will only be able to plot up those points
            # that have an snp classficiation
            if genus == 'Pocillopora':
                label_df_snp = pd.Series(self.poc_snp_kmeans_dict[max_k].labels_, index=self.poc_snp_df.index, name='label')
            else:
                label_df_snp = pd.Series(self.por_snp_kmeans_dict[max_k].labels_, index=self.por_snp_df.index, name='label')

            colours = []
            # Firstly we want to scatter up the samples that aren't found in the snp
            # as light grey 
            scat_samples = [_ for _ in pcoa_df.index if _ not in label_df_snp.index]
            ax[3, misco_i].scatter(pcoa_df.loc[scat_samples,'PC1'], pcoa_df.loc[scat_samples,'PC2'], s=16, c='lightgrey')
            for i in range(kmeans.cluster_centers_.shape[0]): # for each k
                # Get a list of samples that we want to plot up
                # This gives us the rows where the lables are i in the snp
                i_samples = label_df_snp[label_df_snp==i].index
                # Then we want to scatter those samples that are in the pcoa_df
                scat_samples = [_ for _ in pcoa_df.index if _ in i_samples]
                scat = ax[3, misco_i].scatter(pcoa_df.loc[scat_samples,'PC1'], pcoa_df.loc[scat_samples,'PC2'], s=16)
                
                colours.append(scat._original_facecolor[0])
                
            reset_xlim = ax[2,misco_i].get_xlim()
            reset_ylim = ax[2,misco_i].get_ylim()
            current_ax = ax[3,misco_i]
            for i in range(kmeans.cluster_centers_.shape[0]): # for each k
                # plot up the centroid and the points in the same colour
                # Plot the centroid as vline hline intersect
                #vline
                centroid_x = kmeans.cluster_centers_[i][0]
                centorid_y = kmeans.cluster_centers_[i][1]
                
                current_ax.plot([centroid_x, centroid_x],[current_ax.get_ylim()[0], current_ax.get_ylim()[1]], c=colours[i])
                current_ax.plot([current_ax.get_xlim()[0], current_ax.get_xlim()[1]],[centorid_y, centorid_y], c=colours[i])
            
            # lims need resetting as the plotting of the lines changes them.
            current_ax.set_xlim(reset_xlim)
            current_ax.set_ylim(reset_ylim)
            current_ax.set_title(f'k={max_k}')

            # plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster_18s.png'), dpi=600)
            # foo = 'bar'

        plt.tight_layout()
        plt.savefig(os.path.join(self.eighteens_dir, 'temp_fig_cluster_18s.png'), dpi=600)
        foo = 'bar'
        
        
        

if __name__ == "__main__":
    c = Cluster18S()
    # c.visalise_snp_df()
    c.category_scores()