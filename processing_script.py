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
import subprocess

class ITS2Processing:
    def __init__(self):
        self.base_directory_of_sequencing_data = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200116_tara_pacific_its2/"
        self.input_dir = os.path.abspath(os.path.join('.', 'input'))
        self.sample_provenance_path = os.path.join(self.input_dir, 'tara_samples_provenance.csv')
        self.sample_provenance_df = self._make_sample_provenance_df()
        self.cache_dir = os.path.abspath(os.path.join('.', 'cache'))
        self.info_df = self._make_info_df()
        # At this point we have an information dataframe that can be used to create datasheets.
        # This information dataframe contains the information for all of the samples
        # It does not contain information for the negative controls
        # Ideally we want to run the negative samples that correspond to the coral samples
        # in the same dataloading and analysis
        # The negative control sequencing files are all located in a single directory
        # and there is a mapping file that is space delimited. In this file the sequencing file name
        # i.e. TARA_AW-0000064_METAB-ITS2_CP8RR-12BA001-1
        # and then in the second and third column there is an abbreviated version
        # of the two sets of negative control sequencing files that map to the samples
        # i.e. CEB_AWN CEB_AWO. It seems that a very large number of samples
        # map to two negative controls. As such I'm guessing that the negative controls are
        # for the master mixes or something.
        # Because we are on a bit of a tight deadline until the next tara meeting
        # and becuase I don't think these negative are going to be of much use,
        # I'm going to discount them for the time being.

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path)
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        df.rename(columns={'EVENT latitude start (North)':'lat', 'EVENT longitude start (East)': 'lon'}, inplace=True)
        return df

    def _make_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df.pickle')):
            info_df = pickle.load(open(os.path.join(self.cache_dir, 'info_df.pickle'), 'rb'))
        else:
            # TODO add lat and long to this. We will need to make a site and island and OA to lat long dict at some point.
            gicdf = GenerateInfoCollectionDF(base_dir=self.base_directory_of_sequencing_data, provenance_df=self.sample_provenance_df)
            info_df = gicdf.generate_df()

            pickle.dump(info_df, open(os.path.join(self.cache_dir, 'info_df.pickle'), 'wb'))

        return info_df

    def make_sp_data_sheet(self):
        """Run this method to create a .csv file that can then be copy and pasted into a SymPortal
        datasheet template. We will then use this to start the symportal data. We are going to implement
        a new functionality to SymPortal that will allow for full paths to be provided in the datasheet in which
        case the the sequencing file pair will be looked for at that directory else the file name will be looked
        for in the directory that was provided as the argument to --load."""

class GenerateInfoCollectionDF:
    """Class concerned with creating the information dataframe"""
    def __init__(self, base_dir, provenance_df):
        self.base_dir = base_dir
        self.provenance_df = provenance_df
        self.rows = []

    def generate_df(self):
        self._create_data_rows()
        return pd.DataFrame(data=self.rows, columns=['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site', 'lat', 'lon'])
    
    def _create_data_rows(self):
        """ Parse through the seq data file structure gathering 
        inferring the info for each sequecning file pair as we descend the structure/
        Collect a single row of data for each sample"""
        for location in os.listdir(self.base_dir):
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
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=species, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_indi}')

                        elif sample_type == 'PLANKTON':
                            for water_type in os.listdir(parsing_dir_sample_type):
                                if water_type == 'CSW':
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type)
                                    for individual in os.listdir(parsing_dir_water_type):
                                        parsing_dir_indi = os.path.join(parsing_dir_water_type, individual, 'S320')
                                        # now we are in the directory that contains the actual paired fastq.gz files for a
                                        # given water sample
                                        # collect the information we need
                                        srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                        self.rows.append(srg.create_sample_row())
                                        print(f'Processed: {parsing_dir_indi}')

                                elif water_type == 'SURFACE':
                                    # then this is a SURFACE sample and there are no individuals
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type, 'S320')
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_water_type, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_water_type}')

            elif 'OA' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location, 'PLANKTON', 'SURFACE', 'S320')
                # NB the arguments are a little messed up here on purpose as we don't have seperate site and location
                srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_loc, sample_type='OA', site=location, water_species='PLANKTON', provenance_df=self.provenance_df)
                self.rows.append(srg.create_sample_row())
                print(f'Processed: {parsing_dir_loc}')


class SampleRowGenerator:
    def __init__(self, location, parsing_dir, sample_type, site, water_species, provenance_df):
        self.parsing_dir = parsing_dir
        self.location = location
        self.sample_type = sample_type
        self.site = site
        self.water_type = water_species
        self.provenance_df=provenance_df
        self.sample_name = '_'.join(os.listdir(self.parsing_dir)[0].split('_')[:2])
        self.lat = self.provenance_df.at[self.sample_name, 'lat']
        self.lon = self.provenance_df.at[self.sample_name, 'lon']
        self.files = os.listdir(self.parsing_dir)
        if len(self.files) != 2:
            self._concatenate_seq_files()
            self.files = os.listdir(self.parsing_dir)
            if len(self.files) != 2:
                raise RuntimeError(f'Error in concatenation of seq files in {self.parsing_dir}')
        self.fwd_found = False
        self.rev_found = False
        self.fwd_path = None
        self.rev_path = None
        self._get_seq_paths()

    def _concatenate_seq_files(self):
        """ Sometime there is more than 1 set of sequencing files. In this case
        we want to concatenate all of the R1 reads together and all of the R2 reads
        together. Importantly we want to merge the respective memebers of the sequencing
        pairs in the same order.
        We will get a list of R1 files in the director (arbitrary order)
        We will then get a corresponding list of R2 files in the same realtive order
        by doing a replace function on the R1 list.
        We will then concatenate each list of files. Importantly, we will do this
        directly on the .gz files as this is apparently OK.
        https://www.biostars.org/p/81924/
        """
        # get files
        r1_files_list = [file_name for file_name in self.files if 'R1' in file_name]
        r2_files_list = [file_name.replace("R1", "R2") for file_name in r1_files_list]
        
        # Create the list that holds the cat command and the arguments
        exe_1 = [os.path.join(self.parsing_dir, _) for _ in r1_files_list]
        exe_1.insert(0, 'cat')
        exe_2 = [os.path.join(self.parsing_dir, _) for _ in r2_files_list]
        exe_2.insert(0, 'cat')
        
        # outpaths that will be used to create a stream that will be passed as stdout
        out_path_1 = os.path.join(self.parsing_dir, r1_files_list[0].replace("R1.fastq.gz", "merged.R1.fastq.gz"))
        out_path_2 = os.path.join(self.parsing_dir, r2_files_list[0].replace("R2.fastq.gz", "merged.R2.fastq.gz"))

        # do cat
        with open(out_path_1, 'wb') as f:
            subprocess.run(exe_1, stdout=f)
        with open(out_path_2, 'wb') as f:
            subprocess.run(exe_2, stdout=f)
        
        # rm the old files that have now been concatenated
        for file_path in exe_1[1:]:
            os.remove(file_path)
        for file_path in exe_2[1:]:
            os.remove(file_path)

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
            return [self.sample_name, self.fwd_path, self.rev_path, self.sample_type, self.water_type, self.location, self.site, self.lat, self.lon]

its2processing = ITS2Processing()