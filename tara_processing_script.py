import pandas as pd
import os
import sys
import pickle


def figure_making_corals():
    '''The aim of this will be to make a figure that splits up the samples into species, site and island
    so that we can get a better idea about how the symbiodiniaceae changes across the locations and species.'''

    # I think we should be able to hack some of the code used in the symportal repo for this.
    # we esentially just want to re-jig the symportal plotting outputs so that we have the same plotting order of samples
    # which are based on the types and then we can just plot subsamples of these according to the order of
    # islands, then site then species. I think we can plot 3 levels of lines below each plot that can indicate
    # which island, site and species they are from.

    
    return

def write_out_stats_and_reorder_files():
    ''' The purpose of this will be to produce a couple of very basic stats of how many samples
    we have for the different types and then to write out the files into two directories.
    We will write out the coral samples into one directory and the water samples into another.
    We do this as the coral samples will go all the way through the SP analysis and have predicted.
    The water samples can just be submitted to the db and we can get the sequences through QC. They
    will still need further processing but this will be more tertiary analysis.
    '''

    info_df = generate_info_df_for_samples()

    # one off adding of two samples to the df
    # sample one

    # sample_name_one = 'CO0002044'
    # fastq_fwd_file_path_one = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific/ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL026/CS4L/CO0002044_BG8KK_12BA152_1_R1.fastq.gz'
    # fastq_rev_file_path_one = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific/ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL026/CS4L/CO0002044_BG8KK_12BA152_1_R2.fastq.gz'
    # coral_plankton_one = 'CORAL'
    # spp_water_one = 'MILLEPORA'
    # location_one = 'ISLAND06'
    # site_one = 'SITE01'
    # s_one = pd.Series(data=[sample_name_one, fastq_fwd_file_path_one, fastq_rev_file_path_one, coral_plankton_one, spp_water_one, location_one, site_one], name=sample_name_one, index=list(info_df))
    #
    # sample_name_two = 'CO0002041'
    # fastq_fwd_file_path_two = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific/ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL029/CS4L/CO0002041_H3YMJBCX2_12BA122_2_R1.fastq.gz'
    # fastq_rev_file_path_two = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific/ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL029/CS4L/CO0002041_H3YMJBCX2_12BA122_2_R2.fastq.gz'
    # coral_plankton_two = 'CORAL'
    # spp_water_two = 'MILLEPORA'
    # location_two = 'ISLAND06'
    # site_two = 'SITE01'
    #
    # s_two = pd.Series(
    #     data=[sample_name_two, fastq_fwd_file_path_two, fastq_rev_file_path_two, coral_plankton_two, spp_water_two,
    #           location_two, site_two], name=sample_name_two, index=list(info_df))
    #
    # info_df = info_df.append(s_one)
    # info_df = info_df.append(s_two)
    # apples = 'asdf'
    # pickle.dump(info_df, open('{}/info_df.pickle'.format(os.getcwd()), 'wb'))

    coral_output_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20181014_tara_initial_corals'

    non_coral_output_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20181014_tara_initial_non_corals'

    # get the coral rows
    coral = info_df.loc[info_df['coral_plankton'] == 'CORAL']
    num_coral = len(coral.index.values.tolist())

    # get the OA
    OA = info_df.loc[info_df['coral_plankton'] == 'OA']
    num_oa = len(OA.index.values.tolist())

    # CSW
    csw = info_df.loc[info_df['spp_water'] == 'CSW']
    num_csw = len(csw.index.values.tolist())

    # surface
    surface = info_df.loc[info_df['spp_water'] == 'SURFACE']
    num_surface = len(surface.index.values.tolist())

    # # num islands
    # num_islands = len(set([loc for loc in info_df.loc[:,'location'].values.toslist() if 'ISLAND' in loc]))

    # # let's now move the samples to their new directories os.rename is the command for this
    #
    # #coral
    # move_files_to_new_dir(coral_output_dir, coral)
    #
    # #oa
    # move_files_to_new_dir(non_coral_output_dir, OA)
    #
    # #CSW
    # move_files_to_new_dir(non_coral_output_dir, csw)
    #
    # #surface
    # move_files_to_new_dir(non_coral_output_dir, surface)



    apples = 'asdf'

    return

def move_files_to_new_dir(new_dir_base, df ):
    for ind in df.index.values.tolist():
        indi_series = df.loc[ind]
        fwd_file_name = indi_series['fastq_fwd_file_path'].split('/')[-1]
        fwd_current_path = indi_series['fastq_fwd_file_path']
        fwd_destination_path = '{}/{}'.format(new_dir_base, fwd_file_name)

        os.rename(fwd_current_path, fwd_destination_path)

        rev_file_name = indi_series['fastq_rev_file_path'].split('/')[-1]
        rev_current_path = indi_series['fastq_rev_file_path']
        rev_destination_path = '{}/{}'.format(new_dir_base, rev_file_name)

        os.rename(rev_current_path, rev_destination_path)

def generate_info_df_for_samples():

    if os.path.isfile('{}/info_df.pickle'.format(os.getcwd())):
        info_df = pickle.load(open('{}/info_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        tara_data_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180907_tara_pacific'


        # info_df.colums = ['sample_name', 'fastq_fwd_file_name', 'fastq_rev_file_name', 'sample_type', 'host_phylum',
        #                   'host_class', 'host_order', 'host_family', 'host_genus', 'host_species', 'collection_latitude',
        #                   'collection_longitude', 'collection_data', 'collection_depth']

        # lets create a dict where the key will be the sample_name and the value will be dict with each of the above values
        info_collection_dict = {}

        # now lets parse through the directories using them to get some of the information facts above
        generate_info_collection_dict(info_collection_dict, tara_data_dir)

        # here we should have the info_collection_dict populated. We can now turn each of these into series and then
        # add them to the info_df
        columns_for_df = ['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site']

        info_df = create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict)

        pickle.dump(info_df, open('{}/info_df.pickle'.format(os.getcwd()), 'wb'))

        # here we have a dataframe with all of the samples in it
        # we can now move on to do the analysis of them.
        # this hsould be done in a separte method
    return info_df


def create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict):
    series_list = []
    for sample_key, sample_info_dict in info_collection_dict.items():
        data_for_series = [sample_info_dict[ind] for ind in columns_for_df]
        temp_series = pd.Series(data_for_series, index=columns_for_df, name=sample_key)
        series_list.append(temp_series)
    # now we can populate the info df using the series list
    info_df = pd.DataFrame.from_items([(s.name, s) for s in series_list]).T
    return info_df


def generate_info_collection_dict(info_collection_dict, tara_data_dir):
    for location in os.listdir(tara_data_dir):
        if 'ISLAND' in location:
            parsing_dir_loc = '{}/{}'.format(tara_data_dir, location)
            for site in os.listdir(parsing_dir_loc):
                parsing_dir_site = '{}/{}'.format(parsing_dir_loc, site)
                for sample_type in os.listdir(parsing_dir_site):
                    parsing_dir_sample_type = '{}/{}'.format(parsing_dir_site, sample_type)
                    if sample_type == 'CORAL':
                        for species in os.listdir(parsing_dir_sample_type):
                            parsing_dir_species = '{}/{}'.format(parsing_dir_sample_type, species)
                            for individual in os.listdir(parsing_dir_species):
                                parsing_dir_indi = '{}/{}/CS4L'.format(parsing_dir_species, individual)
                                # now we are in the directory that contains the actual paired fastq.gz files for a
                                # given coral individual
                                # collect the information we need
                                create_sample_dict(location, parsing_dir_indi, sample_type, site,
                                                                              species, info_collection_dict)

                    elif sample_type == 'PLANKTON':
                        for water_type in os.listdir(parsing_dir_sample_type):
                            if water_type == 'CSW':
                                parsing_dir_water_type = '{}/{}'.format(parsing_dir_sample_type, water_type)
                                for individual in os.listdir(parsing_dir_water_type):
                                    parsing_dir_indi = '{}/{}/S320'.format(parsing_dir_water_type, individual)

                                    # now we are in the directory that contains the actual paired fastq.gz files for a
                                    # given water sample

                                    # collect the information we need
                                    create_sample_dict(location, parsing_dir_indi, sample_type, site,
                                                                                  water_type, info_collection_dict)


                            elif water_type == 'SURFACE':
                                # then this is a SURFACE sample and there are no individuals
                                parsing_dir_water_type = '{}/{}/S320'.format(parsing_dir_sample_type, water_type)

                                # collect the information we need
                                create_sample_dict(location, parsing_dir_water_type, sample_type, site, water_type, info_collection_dict)


        elif 'OA' in location:
            parsing_dir_loc = '{}/{}/PLANKTON/SURFACE/S320/'.format(tara_data_dir, location)
            # here we are already in a directory that contains the actual paried fastq.gz files
            # SAMPLE NAME
            sample_name = os.listdir(parsing_dir_loc)[0].split('_')[0]
            # FWD and REV PATHS
            files = os.listdir(parsing_dir_loc)
            if len(files) != 2:
                print('more than 2 files in individual\'s directory')
                sys.exit(1)
            fwd_found = False
            rev_found = False
            for file_name in files:
                if 'R1' in file_name:
                    fwd_path = '{}/{}'.format(parsing_dir_loc, file_name)
                    fwd_found = True
                elif 'R2' in file_name:
                    rev_path = '{}/{}'.format(parsing_dir_loc, file_name)
                    rev_found = True
            # make sure that both the fwd and rev paths have been identified
            if not fwd_found or not rev_found:
                print('fwd or rev read not found')
                sys.exit(1)
            sample_dict = {'sample_name': sample_name,
                           'fastq_fwd_file_path': fwd_path,
                           'fastq_rev_file_path': rev_path,
                           'coral_plankton': 'OA',
                           'spp_water': 'PLANKTON',
                           'location': location,
                           'site': location}
            info_collection_dict[sample_name] = sample_dict
    return

def create_sample_dict(location, parsing_dir, sample_type, site, water_type, info_collection_dict):
    # SAMPLE NAME
    sample_name = os.listdir(parsing_dir)[0].split('_')[0]
    # FWD and REV PATHS
    files = os.listdir(parsing_dir)
    if len(files) != 2:
        # we don't know why this is happening and it only happens in two samples
        # ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL029/CS4L/
        # ISLAND06/SITE01/CORAL/MILLEPORA/INDIVIDUAL026/CS4L
        # so for the time being we will just ignore this sample
        return

        # # sometimes there appear to be two sets of files instead of just the one
        # # for the time being I will put these down as two seperate samples and give them names with _0 and _1
        # if len(files) == 4:
        #     sample_base_0 = '_'.join(files[0].split('_')[:-1])
        #     print('more than 2 files in individual\'s directory')
        #
        #     fwd_found_one = False
        #     rev_found_one = False
        #     fwd_found_two = False
        #     rev_found_two = False
        #     for file_name in files:
        #         if sample_base_0 in file_name:
        #             if 'R1' in file_name:
        #                 fwd_path_one = '{}/{}'.format(parsing_dir, file_name)
        #                 fwd_found_one = True
        #             elif 'R2' in file_name:
        #                 rev_path_one = '{}/{}'.format(parsing_dir, file_name)
        #                 rev_found_one = True
        #         else:
        #             if 'R1' in file_name:
        #                 fwd_path_two = '{}/{}'.format(parsing_dir, file_name)
        #                 fwd_found_two = True
        #             elif 'R2' in file_name:
        #                 rev_path_two = '{}/{}'.format(parsing_dir, file_name)
        #                 rev_found_two = True
        #     if not fwd_found_one or not rev_found_one or not fwd_found_two or not rev_found_two:
        #         print('fwd or rev read not found')
        #         sys.exit(1)
        #     for i, tup in enumerate([(fwd_path_one, rev_path_one), (fwd_path_two, rev_path_two)]):
        #
        #         sample_dict = {'sample_name': '{}_{}'.format(sample_name, i),
        #                        'fastq_fwd_file_path': tup[0],
        #                        'fastq_rev_file_path': tup[1],
        #                        'coral_plankton': sample_type,
        #                        'spp_water': water_type,
        #                        'location': location,
        #                        'site': site}
        #         info_collection_dict['{}_{}'.format(sample_name, i)] = sample_dict
        # else:
        #     print('weird number of files found in directory {}'.format(parsing_dir))
        #     sys.exit(1)
    else:
        fwd_found = False
        rev_found = False
        for file_name in files:
            if 'R1' in file_name:
                fwd_path = '{}/{}'.format(parsing_dir, file_name)
                fwd_found = True
            elif 'R2' in file_name:
                rev_path = '{}/{}'.format(parsing_dir, file_name)
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
        info_collection_dict[sample_name] = sample_dict
    return

write_out_stats_and_reorder_files()