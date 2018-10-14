import pandas as pd
import os
import sys
import pickle

def write_out_stats_and_reorder_files():
    ''' The purpose of this will be to produce a couple of very basic stats of how many samples
    we have for the different types and then to write out the files into two directories.
    We will write out the coral samples into one directory and the water samples into another.
    We do this as the coral samples will go all the way through the SP analysis and have predicted.
    The water samples can just be submitted to the db and we can get the sequences through QC. They
    will still need further processing but this will be more tertiary analysis.
    '''
    return
def generate_info_df_for_samples:
    if os.path.isfile('{}/info_df.pickle'.format(os.getcwd())):
        info_df = pickle.load(open('{}/info_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        tara_data_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific'

        info_df = pd.DataFrame()
        # info_df.colums = ['sample_name', 'fastq_fwd_file_name', 'fastq_rev_file_name', 'sample_type', 'host_phylum',
        #                   'host_class', 'host_order', 'host_family', 'host_genus', 'host_species', 'collection_latitude',
        #                   'collection_longitude', 'collection_data', 'collection_depth']
        columns_for_df = ['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site']
        info_df.colums = columns_for_df

        # lets create a dict where the key will be the sample_name and the value will be dict with each of the above values
        info_collection_dict = {}

        # now lets parse through the directories using them to get some of the information facts above
        generate_info_collection_dict(info_collection_dict, tara_data_dir)

        # here we should have the info_collection_dict populated. We can now turn each of these into series and then
        # add them to the info_df
        info_df = create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict, info_df)

        pickle.dump(info_df, open('{}/info_df.pickle'.format(os.getcwd()), 'wb'))

        # here we have a dataframe with all of the samples in it
        # we can now move on to do the analysis of them.
        # this hsould be done in a separte method
    return info_df


def create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict, info_df):
    series_list = []
    for sample_key, sample_info_dict in info_collection_dict.items():
        data_for_series = [sample_info_dict[ind] for ind in columns_for_df]
        temp_series = pd.Series(data_for_series, index=columns_for_df, name=sample_key)
    # now we can populate the info df using the series list
    info_df = pd.DataFrame.from_items([(s.name, s) for s in series_list]).T
    return info_df


def generate_info_collection_dict(info_collection_dict, tara_data_dir):
    for location in os.listdir(tara_data_dir):
        parsing_dir = '{}/{}'.format(tara_data_dir, location)
        for site in os.listdir(parsing_dir):
            parsing_dir = '{}/{}'.format(parsing_dir, site)
            for sample_type in parsing_dir:
                parsing_dir = '{}/{}'.format(parsing_dir, sample_type)
                if sample_type == 'CORAL':
                    for species in parsing_dir:
                        parsing_dir = '{}/{}'.format(parsing_dir, species)
                        for individual in parsing_dir:
                            parsing_dir = '{}/{}/CS4L'.format(parsing_dir, individual)
                            # now we are in the directory that contains the actual paired fastq.gz files for a
                            # given coral individual
                            # collect the information we need
                            sample_dict, sample_name = create_sample_dict(location, parsing_dir, sample_type, site,
                                                                          water_type)
                            info_collection_dict[sample_name] = sample_dict
                elif sample_type == 'PLANKTON':
                    for water_type in parsing_dir:
                        parsing_dir = '{}/{}'.format(parsing_dir, water_type)
                        for individual in parsing_dir:
                            parsing_dir = '{}/{}/S320'.format(parsing_dir, individual)

                            # now we are in the directory that contains the actual paired fastq.gz files for a
                            # given water sample

                            # collect the information we need
                            sample_dict, sample_name = create_sample_dict(location, parsing_dir, sample_type, site,
                                                                          water_type)

                            info_collection_dict[sample_name] = sample_dict


def create_sample_dict(location, parsing_dir, sample_type, site, water_type):
    # SAMPLE NAME
    sample_name = os.listdir(parsing_dir)[0].split('_')[0]
    # FWD and REV PATHS
    files = os.listdir(parsing_dir)
    if len(files) != 2:
        print('more than 2 files in individual\'s directory')
        sys.exit(1)
    fwd_found = False
    rev_found = False
    for file_name in files:
        if 'R1' in file_name:
            fwd_path = file_name
            fwd_found == True
        elif 'R2' in file_name:
            rev_path = file_name
            rev_found = True
    # make sure that both the fwd and rev paths have been identified
    if not fwd_found or not rev_found:
        print('fwd or rev read not found')
        sys.exit(1)
    sample_dict = {'sample_name': sample_name,
                   'fastq_fwd_file_path': fwd_path,
                   'fastq_rev_file_path': rev_path,
                   'coral_plankton': sample_type,
                   'spp_water': water_type,
                   'location': location,
                   'site': site}
    return sample_dict, sample_name

generate_info_df_for_samples()