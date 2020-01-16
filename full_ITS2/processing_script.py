"""
Script to cover all of the processing for the TARA ITS2 samples from the corals and from the surface waters
It will use a large part of the code that was used when we had only the data from the first three islands
available. This previous script was call tara_processing_script.py

To start with we will want to create an info df for the samples. This in turn can be used to generate a
datasheet that we will then be able to use to run the samples through SymPortal.

This is probaly a good stage to get to for a first attempt
"""
import os
import pickle
import sys
import pandas as pd

class ITS2Processing:
    def __init__(self):
        self.base_directory_of_sequencing_data = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200116_tara_pacific_its2/"
        self.cache_dir = os.path.abspath(os.path.join('.', 'cache'))
        self.info_df = self._make_info_df()

    def _make_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df.pickle')):
            info_df = pickle.load(open('{}/info_df.pickle'.format(os.getcwd()), 'rb'))
        else:
            # TODO add lat and long to this. We will need to make a site and island and OA to lat long dict at some point.
            gicdf = GenerateInfoCollectionDF(parent=self)
            info_df = gicdf.generate_df()

            pickle.dump(info_df, open('{}/info_df.pickle'.format(os.getcwd()), 'wb'))

            # here we have a dataframe with all of the samples in it
            # we can now move on to do the analysis of them.
            # this hsould be done in a separte method
        
        # Now add a column that gives us the coral number
        individual_list = []
        for i, ind in enumerate(info_df.index.values.tolist()):
            if 'CORAL' in info_df.at[ind, 'fastq_fwd_file_path']:
                # add extract and add the individual info
                individual_list.append(int(info_df.at[ind, 'fastq_fwd_file_path'].split('/')[11][-3:]))
            else:
            # add a 'NONE'
                individual_list.append(0)
        info_df['individual'] = individual_list
        return info_df

class GenerateInfoCollectionDF:
    """Class concerned with creating the information dataframe"""
    def __init__(self, parent):
        self.parent = parent
        self.base_dir = parent.base_directory_of_sequencing_data
        self.rows = []

    def generate_df(self):
        self._create_data_rows()
        return pd.DataFrame(data=self.rows, columns=['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site'])
    
    def _create_data_rows(self):
        """ Parse through the seq data file structure gathering 
        inferring the info for each sequecning file pair as we descend the structure/
        Collect a single row of data for each sample"""
        for location in os.listdir(self.parent.base_directory_of_sequencing_data):
            if 'ISLAND' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location)
                for site in os.listdir(parsing_dir_loc):
                    parsing_dir_site = os.path.join(parsing_dir_loc, site)
                    for sample_type in os.listdir(parsing_dir_site):
                        parsing_dir_sample_type = os.path.join(parsing_dir_site, sample_type)
                        if sample_type == 'CORAL':
                            for species in os.listdir(parsing_dir_sample_type):
                                parsing_dir_species = os.path.join(parsing_dir_sample_type, species)
                                for individual in os.listdir(parsing_dir_species):
                                    parsing_dir_indi = os.path.join(parsing_dir_species, individual, 'CS4L')
                                    # now we are in the directory that contains the actual paired fastq.gz files for a
                                    # given coral individual
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=species)
                                    self.rows.append(srg.create_sample_row())

                        elif sample_type == 'PLANKTON':
                            for water_type in os.listdir(parsing_dir_sample_type):
                                if water_type == 'CSW':
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type)
                                    for individual in os.listdir(parsing_dir_water_type):
                                        parsing_dir_indi = os.path.join(parsing_dir_water_type, individual, 'S320')
                                        # now we are in the directory that contains the actual paired fastq.gz files for a
                                        # given water sample
                                        # collect the information we need
                                        srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=water_type)
                                        self.rows.append(srg.create_sample_row())

                                elif water_type == 'SURFACE':
                                    # then this is a SURFACE sample and there are no individuals
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type, 'S320')
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_water_type, sample_type=sample_type, site=site, water_species=water_type)
                                    self.rows.append(srg.create_sample_row())

            elif 'OA' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location, 'PLANKTON', 'SURFACE', 'S320')
                # NB the arguments are a little messed up here on purpose as we don't have seperate site and location
                srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_water_type, sample_type='OA', site=location, water_species='PLANKTON')
                self.rows.append(srg.create_sample_row())

class SampleRowGenerator:
    def __init__(self, location, parsing_dir, sample_type, site, water_species):
        self.parsing_dir = parsing_dir
        self.location = location
        self.sample_type = sample_type
        self.site = site
        self.water_type = water_species
        self.sample_name = os.listdir(self.parsing_dir)[0].split('_')[0]
        self.files = os.listdir(self.parsing_dir)
        if len(self.files) != 2:
            raise NotImplementedError('Implement a function that concatenates multiple fwd and rev')
        self.fwd_found = False
        self.rev_found = False
        self.fwd_path = None
        self.rev_path = None
        self._get_seq_paths()

    def _get_seq_paths(self):
        for file_name in self.files:
                if 'R1' in file_name:
                    self.fwd_path = os.path.join(self.parsing_dir, file_name)
                    self.fwd_found = True
                elif 'R2' in file_name:
                    self.rev_path = os.path.join(self.parsing_dir, file_name)
                    self.rev_found = True

    def create_sample_row(self):
        if not self.fwd_found or not self.rev_found:
            print('fwd or rev read not found')
            sys.exit(1)
        else:
            return [self.sample_name, self.fwd_path, self.rev_path, self.sample_type, self.water_type, self.location, self.site]