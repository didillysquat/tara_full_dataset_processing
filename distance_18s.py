from base_18s import EighteenSBase
import os
import compress_pickle

class EighteenSDistance(EighteenSBase):
    """
    Class concerned with the calculation and presentation of the between sample
    distances calculated based on the low level sequence diversity found within each of the samples.
    We will work with UniFrac and BrayCurtis versions of all of the distances we calculate.
    Eventually we will want to be comparing these distances to the clusters that Didier is going
    to produce. However, for the time being we will just look to see if we can find any natural clustering patterns.
    I think we will likely have to work in a somewhat nested manner of fashion.
    Fistly we should split by host species. From this we should verify that the samples that have
    the secondary sequences (i.e. different most abundant sequences to the 'normal' most abundant sequence for the
    genus) can be differentiated in the PCoA. It should be quite obvious.
    If this is the case then the second set of distances that we should compute should be seperated according to 
    the most abundant sequences (i.e. seperate out those samples that have the secondary sequences as their most abunant sequence).
    Once we have done this we should look at the plots and see how the structuring is looking.
    """
    def __init__(self):
        super().__init__()
        # We will need to which coral species the consolidated coral dict belongs
        # to for each of the samples. This differs from the TARA logsheet annotation
        # sometimes. We will add this information to the info_df and then pickle out
        # a new info_df so that we only have to do this once. We will also need to know
        # what the most abundant sequence is in the consolidated host dictionary for 
        # each sample so that we can later differentiate according to this and also
        # check for separation accorrding to this on the first PCoA that we will produce
        self._add_additional_info_to_info_df()
        
    def _add_additional_info_to_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz')):
            self.info_df = compress_pickle.load(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))
        else:
            most_abund_coral_genus_df_list = []
            most_abund_seq_of_coral_genus_df_list = []
            for sample_name in self.info_df.index:
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
                coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
                most_abund_coral_genus_df_list.append(self._identify_most_abund_coral_genus(rel_all_seq_abundance_dict, coral_annotation_dict))
                consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                most_abund_seq_of_coral_genus_df_list.append(sorted(_ for _ in consolidated_host_seqs_abund_dict.items(), key=lambda x: x[0], reverse=True)[0][0])
            self.info_df['most_abund_coral_genus'] = most_abund_coral_genus_df_list
            self.info_df['most_abund_seq_of_coral_genus'] = most_abund_seq_of_coral_genus_df_list
            compress_pickle.dump(self.info_df, os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))

    def _identify_most_abund_coral_genus(self, rel_all_seq_abundance_dict, coral_annotation_dict):
        for sorted_tup in sorted(
            [(seq_name, rel_abund) for seq_name, rel_abund in rel_all_seq_abundance_dict.items()], 
            key=lambda x: x[1], 
            reverse=True
            ):
                try:
                    genus = coral_annotation_dict[sorted_tup[0]]
                    if genus == 'Porites':
                        return 'Porites'
                    elif genus == 'Pocillopora':
                        return 'Pocillopora'
                    elif genus == 'Millepora':
                        return 'Millepora'
                except KeyError:
                    continue
    
    def make_and_plot_dist_and_pcoa(self, resolution_type):
        # Call the class that will be responsible for calculating the distances
        # Calculating the PCoA and then plotting
        da = DistanceAnlyses(resolution_type=resolution_type, info_df=self.info_df)
        da.compute_distances()

class DistanceAnlyses:
    def __init__(self, resolution_type, info_df):
        self.resolution_type = resolution_type
        self.info_df = info_df
        # get the category names and list of samples that will be used in calculating
        # the distances for each of those categories. For example when
        # resolution_type == 'host_only', categories will be Porites, Millepora, Pocillopora
        # When host_secondary_seq == e.g. Porites_first, Porites_second etc.
        self.categories, self.sample_lists = _generate_distance_categories()
        self.dist_methods = ['unifrac', 'braycurtis']
    
    def compute_distances(self):
        for distance_cat, sample_list in zip(self.categories, self.sample_lists):
            for dist_method in self.dist_methods:
                indi_dist_creator = IndiDistanceAnalysis()
            
    def _generate_distance_categories(self):
        if self.resolution_type == 'host_only':
            sample_lists = []
            categories = list(info_df['most_abund_coral_genus'].values())
            for category in categories:
                sample_lists.append(list(info_df[info_df['most_abund_coral_genus'] == category].values()))
            return categories, sample_lists

class IndiDistanceAnalysis:
    def __init__(self, distance_cat, sample_list, dist_method, qc_dir):
        """This is the class that will take care of creating a single between sample
        distance method"""
        self.qc_dir = qc_dir
        self.category = distance_cat
        self.samples = sample_list
        self.dist_method = dist_method
        # In the initial 18s we were using a normalisation of 10000 seqs for the braycurtis
        # but 1000 for the unifrac due to the amount of time it was taking to create
        # the trees. We will start here with 1000 seqs for the unifrac and see
        # how long it takes to make the tree. If its not too bad then we can up the number
        if self.dist_method == 'unifrac':
            self.num_seqs_to_normalise_to = 1000
        else:
            self.num_seqs_to_normalise_to = 10000
        self.abundance_df = _create_abundance_df()
        # The pseudo code for this is
        # Create a abundance dataframe (we need to choose whether to normalise this or not)
        # We probably should normalise.
        # Do this by creating a dict of dicts and using the dicts methods that exist
        # in the initial 18S code
        # master fasta file of the sequences from the samples,
        # If we are doing unifrac then we will need to compute a tree
        # then we can do either the braycurtis calculations or
        # the unifrac calculations.
        # then from this produce pcoa
        # then plot this
        # The tree creation will likely be one of the most expensive parts
        # so we should aim to pickle out the trees if possible, same with alignments.
        # We can possibly pickle these items out as their hashes.

    def _create_abundance_df(self):
        """For each sample in self.samples, load up the consolidated_host_seqs_abund_dict
        and create a dict that is sequence to normalised abundance. Then add this dictionary
        to a master dictionary where sample_name is key. Then finally create a df from this
        dict of dicts"""
        dict_to_create_df_from = {}
        for sample_name in samples:
            temp_sample_dict = {}
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            
            for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                normalised_abund = int(rel_abund*self.num_seqs_to_normalise_to)
                if normalised_abund:
                    temp_sample_dict[sequence] = normalised_abund
            dict_to_create_df_from[sample_name] = temp_sample_dict
        
        df = pd.DataFrame.from_dict(dict_to_create_df_from, orient='index')
        df[pd.isna(self.abundance_df)] = 0
        return df
            

        