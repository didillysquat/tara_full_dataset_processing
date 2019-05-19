#!/usr/bin/env python3.6
import pandas as pd
import os
import sys
import pickle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
from datetime import datetime
import random
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import statistics
import subprocess
from collections import defaultdict

# This is the method to call for plotting the non-symbiodiniacea bar plots of the non-coral samples
def non_symbiodiniaceae_taxonomy_work_do_plotting():
    taxa_abund_dict_holder_dict = pickle.load(open('taxa_abund_dict_holder_dict.pickle', 'rb'))
    taxa_abund_dict_holder_dict_group = pickle.load(open('taxa_abund_dict_holder_dict_group.pickle', 'rb'))
    taxa_abund_dict_holder_dict_match_values = pickle.load(open('taxa_abund_dict_holder_dict_match_values.pickle', 'rb'))

    master_taxa_abund_dict = pickle.load(open('master_taxa_abund_dict.pickle', 'rb'))
    master_taxa_abund_dict_group = pickle.load(open('master_taxa_abund_dict_group.pickle', 'rb'))

    symbiodinium_list = pickle.load(open('symbiodinium_list.pickle', 'rb'))

    sorted_taxa_list = [a[0] for a in sorted(master_taxa_abund_dict.items(), key=lambda x: x[1], reverse=True)]

    info_df = generate_info_df_for_samples()

    # SETUP AXES
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(10, 6))

    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    gs = plt.GridSpec(5, 4, figure=fig, height_ratios=[1, 1, 1, 0.3, 1])

    # The easiest way to do the labeling of the axes is to create an axis that is each of the Gridspec cells
    # this way we can just apply labels and titles to these
    # first make the axes for the main 3 x 3
    labeling_ax = []
    for row_ind in range(3):
        temp_list = []
        for col_ind in range(3):
            ax = plt.subplot(gs[row_ind, col_ind])
            remove_axes_but_allow_labels(ax)
            temp_list.append(ax)
        labeling_ax.append(temp_list)
    OAaxx = plt.subplot(gs[0, 3])
    remove_axes_but_allow_labels(OAaxx)


    # within each of the GrdiSpec subplots we will make a subplotspec which is three plots on one row

    ax_list = []

    grid_spec_subplot_list = []
    for row_ind in range(3):
        for col_ind in range(4):
            if col_ind == 3 and row_ind == 0:
                # this is the axis for the OA samples
                ax = plt.subplot(gs[row_ind, col_ind])
                ax_list.append(ax)
            elif col_ind !=3:
                # put in the main data 3 plots
                temp_grid_spec_subplot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[row_ind, col_ind])
                grid_spec_subplot_list.append(temp_grid_spec_subplot)
                for i in range(2):
                    # NB this might be a 2d array, lets see.
                    ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                    ax_list.append(ax)
                    fig.add_subplot(ax)


    # now do the invisible row that will give us the space we want
    ax_space = plt.subplot(gs[3, :])
    remove_axes_but_allow_labels(ax_space)
    # now split up the final row to put the legend in. One for DIVs and one for TYPEs
    # temp_grid_spec_subplot_leg = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[4, :])
    leg_axes = []
    for i in range(2):
        ax = plt.subplot(gs[4, :-1])
        # ax = plt.Subplot(fig, temp_grid_spec_subplot_leg[i])
        leg_axes.append(ax)
        # fig.add_subplot(ax)

    max_n_cols = 4
    max_n_rows = 7
    num_leg_cells = max_n_cols * max_n_rows
    if os.path.isfile('taxa_colour_dict.pickle'):
        colour_dict_taxa = pickle.load(open('taxa_colour_dict.pickle', 'rb'))
    else:
        colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                              create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000, num_cols=len(master_taxa_abund_dict.keys()),
                                                 time_out_iterations=10000)]

        grey_palette_taxa = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']



        colour_dict_taxa = {}
        for i in range(len(sorted_taxa_list)):
            if i < num_leg_cells:
                colour_dict_taxa[sorted_taxa_list[i]] = colour_palette_pas[i]
            else:
                grey_index = i % len(grey_palette_taxa)
                colour_dict_taxa[sorted_taxa_list[i]] = grey_palette_taxa[grey_index]
        pickle.dump(colour_dict_taxa, open('taxa_colour_dict.pickle', 'wb'))

    ax_count = 0
    extra_ax_count = 0
    for site in ['SITE01', 'SITE02', 'SITE03']:
        for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
            for spp in ['CSW', 'SURFACE']:
                ax = ax_list[ax_count]
                patches_list = []
                ind = 0
                colour_list = []

                # for each set of location, site and spp, we basically want to get a list of the samples
                # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                # to a sample ID using the smp_name_to_smp_id_dict.

                # get sample_names that fit the requirements
                sample_names_of_set = info_df.loc[
                    (info_df['location'] == location) &
                    (info_df['site'] == site) &
                    (info_df['spp_water'] == spp)
                    ].index.values.tolist()

                # convert these to sample IDs
                # The sample names in symportal are actually the full file names version rather than
                # the shorter versions in the info_df. As such we should we will have to do a conversion here
                full_sample_names = [
                    '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for smp_name in
                    sample_names_of_set]


                num_smp_in_this_subplot = len(full_sample_names)
                x_tick_label_list = []

                plot_sample_on_ax_taxa(colour_dict_taxa, colour_list, full_sample_names, ind, patches_list,
                                       sorted_taxa_list, taxa_abund_dict_holder_dict, x_tick_label_list)

                if site=='SITE03':
                    if spp == 'CSW':
                        x_tick_label_list = ['CSW', 'CSW']
                    elif spp == 'SURFACE':
                        x_tick_label_list = ['SURFACE']
                    apply_patches_to_ax_taxa(ax, colour_list, num_smp_in_this_subplot, patches_list, spp,
                                             max_num_smpls_in_subplot=3, x_tick_label_list=x_tick_label_list)
                else:
                    apply_patches_to_ax_taxa(ax, colour_list, num_smp_in_this_subplot, patches_list, spp, max_num_smpls_in_subplot=3)

                ax_count += 1

                # now here we should check to see if we should be plotting the OA samples.
                # this should happen just after we plotted ISLAND 15, site 2, SURFACE
                if spp == 'SURFACE' and location == 'ISLAND15' and site =='SITE01':
                    # get sample_names that fit the requirements
                    sample_names_of_set = info_df.loc[info_df['spp_water'] == 'PLANKTON'].index.values.tolist()

                    # convert these to sample IDs
                    # The sample names in symportal are actually the full file names version rather than
                    # the shorter versions in the info_df. As such we should we will have to do a conversion here
                    full_sample_names = [
                        '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for
                        smp_name in sample_names_of_set]

                    num_smp_in_this_subplot = len(full_sample_names)
                    x_tick_label_list = []
                    patches_list = []
                    ind = 0
                    colour_list = []
                    ax = ax_list[ax_count]
                    plot_sample_on_ax_taxa(colour_dict_taxa, colour_list, full_sample_names, ind, patches_list,
                                           sorted_taxa_list, taxa_abund_dict_holder_dict, x_tick_label_list)

                    apply_patches_to_ax_taxa(ax, colour_list, num_smp_in_this_subplot, patches_list, spp,
                                             max_num_smpls_in_subplot=num_smp_in_this_subplot,
                                             x_tick_label_list=x_tick_label_list)

                    ax_count += 1

    #now plot the leg axes
    # plot_div_legend(colour_dict_div=colour_dict_taxa, leg_axes=leg_axes, max_n_cols_div=max_n_cols, max_n_rows_div=max_n_rows, num_leg_cells_div=num_leg_cells, ordered_list_of_seqs=sorted_taxa_list)
    plot_type_legend(colour_dict_type=colour_dict_taxa, leg_axes=leg_axes, max_n_cols_type=max_n_cols, max_n_rows_type=max_n_rows, num_leg_cells_type=num_leg_cells, sorted_type_prof_names_by_local_abund=sorted_taxa_list, string_cut_off=18)
    # plot_type_legend(colour_dict_type, leg_axes, max_n_cols_type, max_n_rows_type, num_leg_cells_type,
    #                  sorted_type_prof_names_by_local_abund)

    for site_ind, site_label in enumerate(['SITE01', 'SITE02', 'SITE03']):
        for loc_ind, loc_label in enumerate(['ISLAND06', 'ISLAND10', 'ISLAND15']):
            if site_ind == 0:
                labeling_ax[site_ind][loc_ind].set_title(loc_label)
            if loc_ind == 0:
                labeling_ax[site_ind][loc_ind].set_ylabel(site_label, fontsize='large')
            # labeling_ax[0].set_xlabel('porites', fontsize='medium')
    OAaxx.set_title('OA samples')

    date_time_str = str(datetime.now()).replace(' ', '_').replace(':', '-')


    fig_output_base = '{}/{}'.format(os.getcwd(), date_time_str)
    sys.stdout.write('\nsaving as .svg\n')
    plt.savefig('{}_tara_init_results_non_symbiodiniaceae_bar_plots.svg'.format(fig_output_base))
    sys.stdout.write('\nsaving as .png\n')
    plt.savefig('{}_tara_init_results_non_symbiodiniaceae_bar_plots.png'.format(fig_output_base))


    apples = 'asdf'

def apply_patches_to_ax_taxa(ax, colour_list, num_smp_in_this_subplot, patches_list, spp, max_num_smpls_in_subplot, x_tick_label_list=None):
    # We can try making a custom colour map
    # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
    this_cmap = ListedColormap(colour_list)
    # here we should have a list of Rectangle patches
    # now create the PatchCollection object from the patches_list
    patches_collection = PatchCollection(patches_list, cmap=this_cmap)
    patches_collection.set_array(np.arange(len(patches_list)))
    # if n_subplots is only 1 then we can refer directly to the axarr object
    # else we will need ot reference the correct set of axes with i
    # Add the pathces to the axes
    ax.add_collection(patches_collection)
    ax.autoscale_view()
    ax.figure.canvas.draw()
    # also format the axes.
    # make it so that the x axes is constant length

    ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
    ax.set_ylim(0, 1)
    if x_tick_label_list:
        ax.set_xticks(range(num_smp_in_this_subplot))
        ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

        remove_axes_but_allow_labels(ax, x_tick_label_list)
    else:
        remove_axes_but_allow_labels(ax)
    # as well as getting rid of the top and right axis splines
    # I'd also like to restrict the bottom spine to where there are samples plotted but also
    # maintain the width of the samples
    # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
    # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
    # ax.spines['bottom'].set_visible(False)
    ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

def plot_sample_on_ax_taxa(colour_dict_taxa, colour_list, full_sample_names, ind, patches_list, sorted_taxa_list,
                           taxa_abund_dict_holder_dict, x_tick_label_list):
    for smple_name_to_plot in full_sample_names:
        # General plotting
        sys.stdout.write('\rPlotting sample: {}'.format(smple_name_to_plot))
        x_tick_label_list.append(smple_name_to_plot.split('_')[0])
        # for each sample we will start at 0 for the y and then add the height of each bar to this

        bottom_div = 0
        # for each sequence, create a rect patch
        # the rect will be 1 in width and centered about the ind value.
        for taxa in list(sorted_taxa_list):
            # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
            rel_abund_div = taxa_abund_dict_holder_dict[smple_name_to_plot][taxa]
            if rel_abund_div > 0:
                patches_list.append(
                    Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_taxa[taxa]))
                # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                colour_list.append(colour_dict_taxa[taxa])
                bottom_div += rel_abund_div

        ind += 1

def non_symbiodiniaceae_taxonomy_work_do_blasting():
    # the aim of this metod will be to see what we actually have in the non-Symbidiniaceae sequences
    # I have changed SymPortal so that it now puts out the throwaway non-symbiodiniaceae sequences
    # on a sample by sample basis. This should allow us to go through and look at the taxonomy
    # of whats in them.
    # I think the basic principle will be to do the surface and csw samples in the same layout as the corals but
    # maybe with piecharts or stacks
    # and then to do the OA samples seperately.

    info_df = generate_info_df_for_samples()

    # lets write out the blast directory instructions in the cwd
    # write out the blast locaiton instructions to the current dir.
    # blastnPath = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/bin/blastn'
    ntDBDir = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/'

    # Print out the .ncbirc file in the cwd that shows blast were to find the DBs in question
    ncbirc = ['[BLAST]', 'BLASTDB={}'.format(ntDBDir)]
    with open('{}/.ncbirc'.format(os.getcwd()), 'w') as f:
        for line in ncbirc:
            f.write('{}\n'.format(line))

    # have a list that will hold sequences that have been called symbiodinium just so that we can see
    # if they are mis annotates
    symbiodinium_list = []

    # a dictionary that will hold the total for the binomials
    master_taxa_abund_dict = defaultdict(int)
    # also have a dict for the taxa group
    master_taxa_abund_dict_group = defaultdict(int)

    # also a dict that will hold the sample dicts
    taxa_abund_dict_holder_dict = {}
    taxa_abund_dict_holder_dict_group = {}
    # this dict will hold the average for the match and coverage percentage
    taxa_abund_dict_holder_dict_match_values = {}

    # first do the 6 OA samples
    # get sample_names that fit the requirements
    sample_names_of_set = info_df.loc[info_df['spp_water'] == 'PLANKTON'].index.values.tolist()

    # need to conver these short names to the full length names that will be written out in the
    # throw_away directory
    full_sample_names = [
        '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for
        smp_name in
        sample_names_of_set]

    for sample_name in full_sample_names:
        blast_sample_taxonomy(master_taxa_abund_dict, master_taxa_abund_dict_group, sample_name,
                              taxa_abund_dict_holder_dict, taxa_abund_dict_holder_dict_group, taxa_abund_dict_holder_dict_match_values, symbiodinium_list)

    for site in ['SITE01', 'SITE02', 'SITE03']:
        for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
            for spp in ['CSW', 'SURFACE', 'PLANKTON']:

                # for each set of location, site and spp, we basically want to get a list of the samples
                # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                # to a sample ID using the smp_name_to_smp_id_dict.

                # get sample_names that fit the requirements
                sample_names_of_set = info_df.loc[
                    (info_df['location'] == location) &
                    (info_df['site'] == site) &
                    (info_df['spp_water'] == spp)
                    ].index.values.tolist()

                # need to conver these short names to the full length names that will be written out in the
                # throw_away directory
                full_sample_names = [
                    '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for smp_name in
                    sample_names_of_set]

                for sample_name in full_sample_names:
                    blast_sample_taxonomy(master_taxa_abund_dict, master_taxa_abund_dict_group, sample_name,
                                          taxa_abund_dict_holder_dict, taxa_abund_dict_holder_dict_group, taxa_abund_dict_holder_dict_match_values,
                                          symbiodinium_list)



    # here we should pickle the totally sweet dictionaries that we have made
    pickle.dump(taxa_abund_dict_holder_dict, open('taxa_abund_dict_holder_dict.pickle', 'wb'))
    pickle.dump(taxa_abund_dict_holder_dict_group, open('taxa_abund_dict_holder_dict_group.pickle', 'wb'))
    pickle.dump(taxa_abund_dict_holder_dict_match_values, open('taxa_abund_dict_holder_dict_match_values.pickle', 'wb'))

    pickle.dump(master_taxa_abund_dict, open('master_taxa_abund_dict.pickle', 'wb'))
    pickle.dump(master_taxa_abund_dict_group, open('master_taxa_abund_dict_group.pickle', 'wb'))

    pickle.dump(symbiodinium_list, open('symbiodinium_list.pickle', 'wb'))

    return

def blast_sample_taxonomy(master_taxa_abund_dict, master_taxa_abund_dict_group, sample_name, taxa_abund_dict_holder_dict, taxa_abund_dict_holder_dict_group, taxa_abund_dict_holder_dict_match_values, symbiodinium_list):
    print('Blasting {}'.format(sample_name))
    # these are now names that we can iterate through and look up the fasta and name files accordingly
    # get the fasta
    sample_dir = '{}/throw_awayseqs/{}'.format(os.getcwd(), sample_name)
    path_to_fasta = '{}/{}_throw_away_seqs.fasta'.format(sample_dir, sample_name)
    path_to_name = '{}/{}_throw_away_seqs.name'.format(sample_dir, sample_name)

    with open(path_to_fasta, 'r') as f:
        fasta_file = [line.rstrip() for line in f]

    with open(path_to_name, 'r') as f:
        name_file = [line.rstrip() for line in f]

    # make a name dict that simply returns the abundance of the sequence
    seq_name_to_abundance_dict = {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    # make a fasta dict of seq name to seq
    fasta_dict = {fasta_file[i][1:] : fasta_file[i+1] for i in range(0, len(fasta_file), 2)}

    # here we have the fasta_file that we are going to blast.
    outputpath = '{}/blast.out'.format(sample_dir)
    outputFmt = '6 qseqid sseqid staxids sscinames sblastnames evalue pident qcovs'
    # Run the blast command
    # we need to make sure that we are supplying the absolute paths to the inputs and outputs as
    # we are going to be running the acutal blast command the from the cwd rather than the
    # samples directory
    completedProcess = subprocess.run(
        ['blastn', '-out', outputpath, '-outfmt', outputFmt, '-query', path_to_fasta, '-db', 'nt',
         '-max_target_seqs',
         '10', '-num_threads', '20'])
    # Read in the blast output file and perform the taxa counts
    with open(outputpath, 'r') as f:
        blast_output = [line.rstrip() for line in f]
    # make a default dict that will be binomial to abundance
    # for each line in the blast_output
    # create a dict from the blast ouput where key is the sequences name and the results go in
    blast_out_dict = defaultdict(list)
    for line in blast_output:
        components = line.split('\t')
        blast_out_dict[components[0]].append(components[1:])
    # now we can go key by key in the blast_out_dict and collect counts for what binomials we can
    taxa_count_dict = defaultdict(int)
    taxa_count_dict_group = defaultdict(int)
    taxa_count_dict_match_values = []
    total_seqs = sum(seq_name_to_abundance_dict.values())
    for seq_key, blast_list in blast_out_dict.items():
        specific_binomial_found = False
        # the blast_list needs sorting according to it e-value
        # go through the list and get a list of the evalues or put 0 if not found
        # create list of tuples where first value is the evalue and second the list of items
        # then sort this and work through the new list below
        new_tup_list = []
        for sub_list in blast_list:
            try:
                evalue = sub_list[4].split('-')[-1]
                new_tup_list.append((evalue, sub_list))
            except:
                continue
        # now sort the new_typ list
        sorted_list_of_lists = [a[1] for a in sorted(new_tup_list, key=lambda x: x[0], reverse=True)]
        for result_list in sorted_list_of_lists:
            binomial = result_list[2]
            group = result_list[3]
            ident = float(result_list[5])
            cov = float(result_list[6])
            if binomial != 'uncultured eukaryote' and binomial != 'N/A':
                # then we should use this binomial
                if binomial.split(' ')[0] in ['symbiodinium', 'Symbiodinium', 'Symbiodiniaceae', 'symbiodiniaceae', 'Cladocopium',
                                'cladocopium', 'Fugacium', 'fugacium', 'Breviolum', 'breviolum',
                                'Durusdinium', 'durussinium', 'Gerakladium', 'gerakladium', 'Effrenium', 'effrenium']:
                    symbiodinium_list.extend(['>{}'.format(seq_key), fasta_dict[seq_key]])

                rel_abund_of_seq = seq_name_to_abundance_dict[seq_key] / total_seqs
                taxa_count_dict[binomial.split(' ')[0]] += rel_abund_of_seq
                taxa_count_dict_group[group] += rel_abund_of_seq
                # we will multiply the identify and coverage values by the relative abundance of the sequence
                # when summed this should give us an accurate average
                rel_ident = ident * rel_abund_of_seq
                rel_cov = cov * rel_abund_of_seq
                taxa_count_dict_match_values.append((rel_ident, rel_cov))
                master_taxa_abund_dict[binomial.split(' ')[0]] += rel_abund_of_seq
                master_taxa_abund_dict_group[group] += rel_abund_of_seq

                specific_binomial_found = True
                break
        # if we get here then we didn't find a binomial that wasn't uncultured eukaryote
        if not specific_binomial_found:
            taxa_count_dict['uncultured eukaryote'] += seq_name_to_abundance_dict[seq_key] / total_seqs
            taxa_count_dict_group['uncultured eukaryote'] += seq_name_to_abundance_dict[seq_key] / total_seqs
            master_taxa_abund_dict['uncultured eukaryote'] += seq_name_to_abundance_dict[seq_key] / total_seqs
            master_taxa_abund_dict_group['uncultured eukaryote'] += seq_name_to_abundance_dict[seq_key] / total_seqs
    taxa_abund_dict_holder_dict[sample_name] = taxa_count_dict
    taxa_abund_dict_holder_dict_group[sample_name] = taxa_count_dict_group
    taxa_abund_dict_holder_dict_match_values[sample_name] = taxa_count_dict_match_values
    apples = 'asdf'



# This is the method to call for plotting the coral and csw and surface samples DIVs and ITS2 type profiles
def figure_making_bar_plots_corals():
    '''The aim of this will be to make a figure that splits up the samples into species, site and island
    so that we can get a better idea about how the symbiodiniaceae changes across the locations and species.'''

    # I think we should be able to hack some of the code used in the symportal repo for this.
    # we esentially just want to re-jig the symportal plotting outputs so that we have the same plotting order of samples
    # which are based on the types and then we can just plot subsamples of these according to the order of
    # islands, then site then species. I think we can plot 3 levels of lines below each plot that can indicate
    # which island, site and species they are from.

    info_df = generate_info_df_for_samples()

    # this dictionary will hold the details of which corals are related to the csw samples
    # sadly I put together this list by hand so if we come to scale this up you will likely have to code the extraction
    # of these details. Shouldn't be too much work.
    # key is the short name of the sample eg. IW0000XXX or CO000XXXX. value is 'INDIVIDUAL1' or 'INDIVIDUAL10'
    with open('individual_info', 'r') as f:
        indi_info_list = [line.rstrip() for line in f]

    indi_info = {line.split(' ')[0]:line.split(' ')[1] for line in indi_info_list}

    # we want to read in the table that contains both the coral and non-coral sequences so that we can
    # plot the csw and the surface water samples along side the corals
    path_to_tab_delim_rel_count_DIV_coral_non_coral_standalone = '2018-10-21_08-59-37.620726.DIVs.relative.txt'


    path_to_tab_delim_rel_count_type = '34_init_tara_standalone_all_samps_151018_2018-10-21_08-45-56.507454.profiles.relative.txt'

    output_directory = '/home/humebc/projects/tara/initial_its2_processing'

    generate_stacked_bar_data_submission(path_to_tab_delim_rel_count_DIV_coral_non_coral_standalone, path_to_tab_delim_rel_count_type, output_directory, info_df, individual_info=indi_info, time_date_str=None)

    return

def generate_stacked_bar_data_submission(path_to_tab_delim_count_DIV, path_to_tab_delim_count_type, output_directory, info_df, individual_info, time_date_str=None):
    print('Generating stacked bar data submission')
    # /Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/outputs/non_analysis/35.DIVs.relative.txt

    # Here we will generate our standard stacked bar output.
    # We should take into account that we don't know how many samples will be coming through.
    # I think we should aim for a standard width figure but which can get deeper if there are many samples.
    # I.e. we should split up very large sets of samples into multiple plots to keep interpretability
    # as high as possible.

    # I think it would be cool to see the its2 type plotted below the sequences for each of the samples
    # we can also hack the code for this from SP

    # read in the SymPortal relative abundance output
    smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df_div = process_div_df(path_to_tab_delim_count_DIV)

    # now read in a process the type df too.

    colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund, max_n_cols_type, max_n_rows_type, num_leg_cells_type,  = process_type_df(path_to_tab_delim_count_type)

    colour_dict_div, max_n_cols_div, max_n_rows_div, num_leg_cells_div, ordered_list_of_seqs = get_div_colour_dict_and_ordered_list_of_seqs(
        sp_output_df_div)

    # the ordered_list_of_seqs can also be used for the plotting order

    # we should consider doing a plot per clade but for the time being lets start by doing a single plot that will
    # contain all of the clades

    # if we are plotting this in companion with an ITS2 type profile output then we will be passed a
    # sample_order_list. It is very useful to have the ITS2 type profile output figure and the seq figure
    # in the same sample order for direct comparison


    ordered_sample_list = sp_output_df_type.index.values.tolist()
    # pickle.dump(ordered_sample_list, open('ordered_sample_list_for_18s_work.pickle', 'wb'))
    # let's reorder the columns and rows of the sp_output_df according to the sequence sample and sequence
    # order so that plotting the data is easier

    sp_output_df_div = sp_output_df_div[ordered_list_of_seqs]
    # we need to get rid of this slicing so that we don't cut out the csw and surface samples
    # sp_output_df_div = sp_output_df_div.reindex(ordered_sample_list)

    # At this stage we are ready to plot
    # The three following links show how we should be able to construct a list of matplotlib
    # patches (Rectangles in this case) and add these patches to a PatchCollection before finally
    # adding this patch collection to the ax using ax.add_collection().
    # https://matplotlib.org/api/_as_gen/matplotlib.patches.Rectangle.html
    # https://matplotlib.org/examples/api/patch_collection.html
    # https://matplotlib.org/users/artists.html
    # I hope that this will be quicker than using the bar helper sequence by sequence as we normally do
    # It turns out that the colour parameters are ignored from the individual patches when using

    # SETUP AXES
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(14, 10))

    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    gs = plt.GridSpec(5, 6, figure=fig, height_ratios=[1,1,1,0.2,1], width_ratios=[1,0.2,1,0.2,1,0.2])
    # within each of the GrdiSpec subplots we will make a subplotspec which is three plots on one row

    ax_list = []
    extra_ax_list = []
    grid_spec_subplot_list = []
    for row_ind in range(3):
        for col_ind in range(0, 6, 2):
            # put in the main data 3 plots
            temp_grid_spec_subplot = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=gs[row_ind,col_ind])
            grid_spec_subplot_list.append(temp_grid_spec_subplot)
            for i in range(3):
                #NB this might be a 2d array, lets see.
                ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                ax_list.append(ax)
                fig.add_subplot(ax)
            # now put in the csw and surface plots which will be 1,2
            temp_grid_spec_subplot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[row_ind, col_ind + 1])
            grid_spec_subplot_list.append(temp_grid_spec_subplot)
            for i in range(2):
                # NB this might be a 2d array, lets see.
                ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                # we will not add
                extra_ax_list.append(ax)
                fig.add_subplot(ax)

    # now do the invisible row that will give us the space we want
    ax_space = plt.subplot(gs[3, :])
    remove_axes_but_allow_labels(ax_space)
    # now split up the final row to put the legend in. One for DIVs and one for TYPEs
    temp_grid_spec_subplot_leg = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[4, :])
    leg_axes = []
    for i in range(2):
        ax = plt.Subplot(fig, temp_grid_spec_subplot_leg[i])
        leg_axes.append(ax)
        fig.add_subplot(ax)


    # we will leave one subplot empty for making the legend in at the end
    plot_data_axes(ax_list, extra_ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list, smp_id_to_smp_name_dict,
                   smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, individual_info)

    # PLOT DIV LEGEND
    plot_div_legend(colour_dict_div, leg_axes, max_n_cols_div, max_n_rows_div, num_leg_cells_div, ordered_list_of_seqs)

    # PLOT TYPE LEGEND
    plot_type_legend(colour_dict_type, leg_axes, max_n_cols_type, max_n_rows_type, num_leg_cells_type,
                     sorted_type_prof_names_by_local_abund)

    # add the labels and text here so that we don't have to debug through all of the plotting each time
    add_labels(ax_list, leg_axes, extra_ax_list)

    date_time_str = str(datetime.now()).replace(' ', '_').replace(':', '-')

    # plt.tight_layout()
    fig_output_base = '{}/{}'.format(os.getcwd(), date_time_str)
    sys.stdout.write('\nsaving as .svg\n')
    plt.savefig('{}_tara_init_results_coral_bar_plot_with_csw.svg'.format(fig_output_base))
    sys.stdout.write('\nsaving as .png\n')
    plt.savefig('{}_tara_init_results_coral_bar_plot_with_csw.png'.format(fig_output_base))
    # plt.show()
    return

def add_labels(ax_list, leg_axes, extra_ax_list):
    ax_list[1].set_title('ISLAND06')
    ax_list[4].set_title('ISLAND10')
    ax_list[7].set_title('ISLAND15')

    ax_list[0].set_ylabel('SITE 1', fontsize='x-large')
    ax_list[9].set_ylabel('SITE 2', fontsize='x-large')
    ax_list[18].set_ylabel('SITE 3', fontsize='x-large')

    ax_list[18].set_xlabel('porites', fontsize='medium')
    ax_list[19].set_xlabel('pocillopora', fontsize='medium')
    ax_list[20].set_xlabel('millepora', fontsize='medium')

    ax_list[21].set_xlabel('porites', fontsize='medium')
    ax_list[22].set_xlabel('pocillopora', fontsize='medium')
    ax_list[23].set_xlabel('millepora', fontsize='medium')

    ax_list[24].set_xlabel('porites', fontsize='medium')
    ax_list[25].set_xlabel('pocillopora', fontsize='medium')
    ax_list[26].set_xlabel('millepora', fontsize='medium')

    leg_axes[0].set_xlabel('sequence')
    leg_axes[1].set_xlabel('ITS2 type profile')

    extra_ax_list[12].set_xlabel('CSW', horizontalalignment='center', rotation='vertical')
    extra_ax_list[13].set_xlabel('surface', horizontalalignment='center', rotation='vertical')
    extra_ax_list[14].set_xlabel('CSW', horizontalalignment='center', rotation='vertical')
    extra_ax_list[15].set_xlabel('surface', horizontalalignment='center', rotation='vertical')
    extra_ax_list[16].set_xlabel('CSW', horizontalalignment='center', rotation='vertical')
    extra_ax_list[17].set_xlabel('surface', horizontalalignment='center', rotation='vertical')

def plot_type_legend(colour_dict_type, leg_axes, max_n_cols_type, max_n_rows_type, num_leg_cells_type,
                     sorted_type_prof_names_by_local_abund, string_cut_off=10):
    # Since the matplotlib legends are pretty rubbish when made automatically, I vote that we make our own axes
    # all in favour... Ok.
    # Let's plot the boxes and text that are going to make up the legend in another subplot that we will put underneath
    # the one we currenty have. So.. we will add a subplot when we initially create the figure. We will make the axis
    # 100 by 100 just to make our coordinate easy to work with. We can get rid of all of the axes lines and ticks
    # The type names are generally quite long so we will cut the type legends down to 4 x 8
    # we should start plotting in the top left working right and then down
    # until we have completed 100 sequences.
    # Y axis coordinates
    # we will allow a buffer of 0.5 of the legend box's height between each legend box.
    # as such the coordinates of each y will be in increments of 100 / (1.5 * num rows)
    # the depth of the Rectangle for the legend box will be 2/3 * the above.
    y_coord_increments = 100 / (max_n_rows_type)
    leg_box_depth = 2 / 3 * y_coord_increments
    # X axis coordinates
    # for the x axis we will work in sets of three columns were the first col will be for the box
    # and the second and third cols will be for the text
    # as such the x coordinates will be in increments of 100 / (3 * numcols) starting with 0
    # the width of the legend Rectangle will be the above number * 1/6 (I am making this smaller for the types).
    x_coord_increments = 100 / max_n_cols_type
    leg_box_width = x_coord_increments / 6
    # go column by column
    # we can now calculate the actual number of columns and rows we are going to need.
    if len(sorted_type_prof_names_by_local_abund) < num_leg_cells_type:
        if len(sorted_type_prof_names_by_local_abund) % max_n_cols_type != 0:
            n_rows_type = int(len(sorted_type_prof_names_by_local_abund) / max_n_cols_type) + 1
        else:
            n_rows_type = int(len(sorted_type_prof_names_by_local_abund) / max_n_cols_type)
        last_row_len = len(sorted_type_prof_names_by_local_abund) % max_n_cols_type
    else:
        n_rows_type = max_n_rows_type
        last_row_len = max_n_cols_type
    its2_profile_count = 0
    # Once we know the number of rows, we can also adjust the y axis limits
    leg_axes[1].set_xlim(0, 100)
    # axarr[-1].set_ylim(0, 100)
    leg_axes[1].set_ylim(0, ((n_rows_type - 1) * y_coord_increments) + leg_box_depth)
    leg_axes[1].invert_yaxis()
    # If there are more sequences than there are rows x cols then we need to make sure that we are only going
    # to plot the first row x cols number of sequences.
    sys.stdout.write(
        '\nGenerating figure legend for {} most common sequences\n'.format(str(max_n_rows_type * max_n_cols_type)))

    for row_increment in range(min(n_rows_type, max_n_rows_type)):
        # if not in the last row then do a full set of columns
        if row_increment + 1 != n_rows_type:
            for col_increment in range(max_n_cols_type):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[1].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_type[
                                                    sorted_type_prof_names_by_local_abund[its2_profile_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                # lets limit the name to 15 characters and '...'
                if len(sorted_type_prof_names_by_local_abund[its2_profile_count]) > string_cut_off:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count][
                                      :string_cut_off] + '...'
                else:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count]
                leg_axes[1].text(text_x, text_y, text_for_legend,
                                 verticalalignment='center',
                                 fontsize=8)

                # increase the sequence count
                its2_profile_count += 1
        # else just do up to the number of last_row_cols
        else:
            for col_increment in range(last_row_len):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[1].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_type[
                                                    sorted_type_prof_names_by_local_abund[its2_profile_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                # lets limit the name to 15 characters and '...'
                if len(sorted_type_prof_names_by_local_abund[its2_profile_count]) > string_cut_off:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count][
                                      :string_cut_off] + '...'
                else:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count]
                leg_axes[1].text(text_x, text_y, text_for_legend,
                                 verticalalignment='center',
                                 fontsize=8)

                # Increase the sequences count
                its2_profile_count += 1
    remove_axes_but_allow_labels(leg_axes[1])

def plot_div_legend(colour_dict_div, leg_axes, max_n_cols_div, max_n_rows_div, num_leg_cells_div, ordered_list_of_seqs):
    # Since the matplotlib legends are pretty rubbish when made automatically, I vote that we make our own axes
    # all in favour... Ok.
    # Let's plot the boxes and text that are going to make up the legend in another subplot that we will put underneath
    # the one we currenty have. So.. we will add a subplot when we initially create the figure. We will make the axis
    # 100 by 100 just to make our coordinate easy to work with. We can get rid of all of the axes lines and ticks
    # lets aim to plot a 10 by 10 legend max
    # we should start plotting in the top left working right and then down
    # until we have completed 100 sequences.
    # Y axis coordinates
    # we will allow a buffer of 0.5 of the legend box's height between each legend box.
    # as such the coordinates of each y will be in increments of 100 / (1.5 * num rows)
    # the depth of the Rectangle for the legend box will be 2/3 * the above.
    y_coord_increments = 100 / (max_n_rows_div)
    leg_box_depth = 2 / 3 * y_coord_increments
    # X axis coordinates
    # for the x axis we will work in sets of three columns were the first col will be for the box
    # and the second and third cols will be for the text
    # as such the x coordinates will be in increments of 100 / (3 * numcols) starting with 0
    # the width of the legend Rectangle will be the above number * 1/3.
    x_coord_increments = 100 / max_n_cols_div
    leg_box_width = x_coord_increments / 3
    # go column by column
    # we can now calculate the actual number of columns and rows we are going to need.
    if len(ordered_list_of_seqs) < num_leg_cells_div:
        if len(ordered_list_of_seqs) % max_n_cols_div != 0:
            n_rows_div = int(len(ordered_list_of_seqs) / max_n_cols_div) + 1
        else:
            n_rows_div = int(len(ordered_list_of_seqs) / max_n_cols_div)
        last_row_len = len(ordered_list_of_seqs) % max_n_cols_div
    else:
        n_rows_div = max_n_rows_div
        last_row_len = max_n_cols_div
    sequence_count = 0
    # Once we know the number of rows, we can also adjust the y axis limits
    leg_axes[0].set_xlim(0, 100)
    # axarr[-1].set_ylim(0, 100)
    leg_axes[0].set_ylim(0, ((n_rows_div - 1) * y_coord_increments) + leg_box_depth)
    leg_axes[0].invert_yaxis()
    # If there are more sequences than there are rows x cols then we need to make sure that we are only going
    # to plot the first row x cols number of sequences.
    sys.stdout.write(
        '\nGenerating figure legend for {} most common sequences\n'.format(str(max_n_rows_div * max_n_cols_div)))
    for row_increment in range(min(n_rows_div, max_n_rows_div)):
        # if not in the last row then do a full set of columns
        if row_increment + 1 != n_rows_div:
            for col_increment in range(max_n_cols_div):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[0].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_div[ordered_list_of_seqs[sequence_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                leg_axes[0].text(text_x, text_y, ordered_list_of_seqs[sequence_count], verticalalignment='center',
                                 fontsize=8)

                # increase the sequence count
                sequence_count += 1
        # else just do up to the number of last_row_cols
        else:
            for col_increment in range(last_row_len):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[0].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_div[ordered_list_of_seqs[sequence_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                leg_axes[0].text(text_x, text_y, ordered_list_of_seqs[sequence_count], verticalalignment='center',
                                 fontsize=8)

                # Increase the sequences count
                sequence_count += 1
    remove_axes_but_allow_labels(leg_axes[0])

def plot_data_axes(ax_list, extra_ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list, smp_id_to_smp_name_dict,
                   smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, indi_indo):
    ax_count = 0
    extra_ax_count = 0
    for site in ['SITE01', 'SITE02', 'SITE03']:
        for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
            for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                ax = ax_list[ax_count]
                patches_list = []
                ind = 0
                colour_list = []

                # for each set of location, site and spp, we basically want to get a list of the samples
                # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                # to a sample ID using the smp_name_to_smp_id_dict.

                # get sample_names that fit the requirements
                sample_names_of_set = info_df.loc[
                    (info_df['location'] == location) &
                    (info_df['site'] == site) &
                    (info_df['spp_water'] == spp)
                    ].index.values.tolist()

                # temporarily remove CO0002044, CO0002041
                # from the above list
                # if 'CO0002044' or 'CO0002041' in sample_names_of_set:
                #     sample_names_of_set = [name for name in sample_names_of_set if
                #                            name not in ['CO0002044', 'CO0002041']]

                # convert these to sample IDs
                # The sample names in symportal are actually the full file names version rather than
                # the shorter versions in the info_df. As such we should we will have to do a conversion here
                full_sample_names = [
                    '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for smp_name in
                    sample_names_of_set]
                try:
                    smple_ids_of_set = [smp_name_to_smp_id_dict[smp_name] for smp_name in full_sample_names]
                except:
                    apples = 'asdf'
                # now we want to plot in the order of the ordered_sample_list
                ordered_smple_ids_of_set = [smpl_id for smpl_id in ordered_sample_list if smpl_id in smple_ids_of_set]

                coral_csw_x_val_tup_list=None
                if spp == 'POCILLOPORA':
                    # here we need to work out if any of these are csw associated samples
                    # if so then we need to follow through and get the id of them
                    coral_csw_id_one = None
                    coral_csw_id_ten = None
                    for smp_name in sample_names_of_set:
                        if smp_name in indi_indo.keys():
                            # use the same logic below to get the id of this sample
                            temp_full_name = '_'.join(
                                info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])
                            # then this is one of the csw associated samples
                            if indi_indo[smp_name] == 'INDIVIDUAL1':
                                coral_csw_id_one = smp_name_to_smp_id_dict[temp_full_name]
                            elif indi_indo[smp_name] == 'INDIVIDUAL10':
                                coral_csw_id_ten = smp_name_to_smp_id_dict[temp_full_name]

                    # now we need to get the x values of where to draw the line to and from for the csw associated corals
                    # we can work this out by seeing what the index of the sample ids are in the ordered_smple list
                    # as the width we use is 1, the x values will then be the index -+ 0.5
                    coral_csw_x_val_tup_list = []
                    if coral_csw_id_one:
                        if coral_csw_id_one in ordered_smple_ids_of_set:
                            temp_index_of_sample_in_list = ordered_smple_ids_of_set.index(coral_csw_id_one)
                            coral_csw_x_value_tup = (temp_index_of_sample_in_list -0.5, temp_index_of_sample_in_list + 0.5, 'brown')
                            coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)
                    if coral_csw_id_ten:
                        if coral_csw_id_ten in ordered_smple_ids_of_set:
                            temp_index_of_sample_in_list = ordered_smple_ids_of_set.index(coral_csw_id_ten)
                            coral_csw_x_value_tup = (temp_index_of_sample_in_list -0.5, temp_index_of_sample_in_list + 0.5, 'gray')
                            coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)

                num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
                x_tick_label_list = []

                for smple_id_to_plot in ordered_smple_ids_of_set:


                    # General plotting
                    sys.stdout.write('\rPlotting sample: {}'.format(smple_id_to_plot))
                    x_tick_label_list.append(smp_id_to_smp_name_dict[smple_id_to_plot].split('_')[0])
                    # for each sample we will start at 0 for the y and then add the height of each bar to this

                    # PLOT DIVs
                    plot_div_over_type(colour_dict_div, colour_list, ind, patches_list, smple_id_to_plot,
                                       sp_output_df_div)

                    # PLOT type
                    plot_type_under_div(colour_dict_type, colour_list, ind, patches_list, smple_id_to_plot,
                                        sp_output_df_type)
                    ind += 1


                paint_rect_to_axes_div_and_type(ax=ax, colour_list=colour_list,
                                                num_smp_in_this_subplot=num_smp_in_this_subplot,
                                                patches_list=patches_list,
                                                coral_csw_x_val_tup_list=coral_csw_x_val_tup_list,
                                                x_tick_label_list=x_tick_label_list,
                                                max_num_smpls_in_subplot=10)

                ax_count += 1

            # PLOT csw and surface

            # first identify which the csw sample is
            csw_samples = info_df.loc[
                (info_df['location'] == location) &
                (info_df['site'] == site) &
                (info_df['spp_water'] == 'CSW')
                ].index.values.tolist()

            # convert these to sample IDs
            # The sample names in symportal are actually the full file names version rather than
            # the shorter versions in the info_df. As such we should we will have to do a conversion here
            full_sample_names_csw = [
                '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for
                smp_name in csw_samples]
            smple_ids_of_set_csw = [smp_name_to_smp_id_dict[smp_name] for smp_name in full_sample_names_csw]

            # here we need to work out which sample is which individual for the csw association to corals

            coral_csw_id_one = None
            coral_csw_id_ten = None
            for smp_name in csw_samples:
                if smp_name in indi_indo.keys():
                    # use the same logic below to get the id of this sample
                    temp_full_name = '_'.join(
                        info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])
                    # then this is one of the csw associated samples
                    if indi_indo[smp_name] == 'INDIVIDUAL1':
                        coral_csw_id_one = smp_name_to_smp_id_dict[temp_full_name]
                    elif indi_indo[smp_name] == 'INDIVIDUAL10':
                        coral_csw_id_ten = smp_name_to_smp_id_dict[temp_full_name]

            # now we need to get the x values of where to draw the line to and from for the csw associated corals
            # we can work this out by seeing what the index of the sample ids are in the ordered_smple list
            # as the width we use is 1, the x values will then be the index -+ 0.5
            coral_csw_x_val_tup_list = []

            if coral_csw_id_one and coral_csw_id_ten:
                # then both exist
                # check to see that they are in the order 'one' then 'ten'. if not, reverse
                if smple_ids_of_set_csw.index(coral_csw_id_one) == 1:
                    # then we need to reverse
                    smple_ids_of_set_csw = list(reversed(smple_ids_of_set_csw))
                # now add the tuples to the coral_csw_x_val_tup_list
                for i, colour in enumerate(['brown', 'gray']):
                    coral_csw_x_value_tup = (i - 0.5, i + 0.5, colour)
                    coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)


            num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
            x_tick_label_list = []

            colour_list = []
            ind = 0
            patches_list = []
            for smple_id_to_plot in smple_ids_of_set_csw:
                bottom_div = 0
                # for each sequence, create a rect patch
                # the rect will be 1 in width and centered about the ind value.
                for seq in list(sp_output_df_div):
                    # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                    rel_abund_div = sp_output_df_div.loc[smple_id_to_plot, seq]
                    if rel_abund_div > 0:
                        patches_list.append(
                            Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
                        # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                        colour_list.append(colour_dict_div[seq])
                        bottom_div += rel_abund_div
                ind += 1

            paint_rect_to_axes_div_and_type(ax=extra_ax_list[extra_ax_count], colour_list=colour_list,
                                            num_smp_in_this_subplot=2, coral_csw_x_val_tup_list=coral_csw_x_val_tup_list,
                                            patches_list=patches_list, max_num_smpls_in_subplot=2)
            extra_ax_count += 1


            # now get the surface sample
            surface_sample = info_df.loc[
                (info_df['location'] == location) &
                (info_df['site'] == site) &
                (info_df['spp_water'] == 'SURFACE')
                ].index.values.tolist()[0]

            # convert these to sample IDs
            # The sample names in symportal are actually the full file names version rather than
            # the shorter versions in the info_df. As such we should we will have to do a conversion here
            full_sample_name_surface = '_'.join(
                info_df.loc[surface_sample]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])

            smple_id_of_surface_sample = smp_name_to_smp_id_dict[full_sample_name_surface]

            colour_list = []
            bottom_div = 0
            ind = 0
            patches_list = []
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            for seq in list(sp_output_df_div):
                # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                rel_abund_div = sp_output_df_div.loc[smple_id_of_surface_sample, seq]
                if rel_abund_div > 0:
                    patches_list.append(
                        Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
                    # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                    colour_list.append(colour_dict_div[seq])
                    bottom_div += rel_abund_div


            paint_rect_to_axes_div_and_type(ax=extra_ax_list[extra_ax_count], colour_list=colour_list,
                                            num_smp_in_this_subplot=1,
                                            patches_list=patches_list, max_num_smpls_in_subplot=2)
            extra_ax_count += 1

def remove_axes_but_allow_labels(ax, x_tick_label_list=None):
    ax.set_frame_on(False)
    if not x_tick_label_list:
        ax.set_xticks([])
    ax.set_yticks([])

def paint_rect_to_axes_div_and_type(ax, colour_list, num_smp_in_this_subplot,  patches_list, coral_csw_x_val_tup_list=None, x_tick_label_list=None,  max_num_smpls_in_subplot=10):
    # We can try making a custom colour map
    # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
    this_cmap = ListedColormap(colour_list)
    # here we should have a list of Rectangle patches
    # now create the PatchCollection object from the patches_list
    patches_collection = PatchCollection(patches_list, cmap=this_cmap)
    patches_collection.set_array(np.arange(len(patches_list)))
    # if n_subplots is only 1 then we can refer directly to the axarr object
    # else we will need ot reference the correct set of axes with i
    # Add the pathces to the axes
    ax.add_collection(patches_collection)
    ax.autoscale_view()
    ax.figure.canvas.draw()
    # also format the axes.
    # make it so that the x axes is constant length
    ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
    ax.set_ylim(-0.2, 1)
    # ax.set_xticks(range(num_smp_in_this_subplot))
    # ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

    remove_axes_but_allow_labels(ax)

    # as well as getting rid of the top and right axis splines
    # I'd also like to restrict the bottom spine to where there are samples plotted but also
    # maintain the width of the samples
    # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
    # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
    # ax.spines['bottom'].set_visible(False)
    ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

    # once we have added the black line we should add the grey and brown line that will associate the
    # coral samples to the csw samples
    if coral_csw_x_val_tup_list:
        for coral_csw_x_value_tup in coral_csw_x_val_tup_list:
            ax.add_line(Line2D((coral_csw_x_value_tup[0], coral_csw_x_value_tup[1]), (0, 0), linewidth=2, color=coral_csw_x_value_tup[2]))

def plot_div_over_type(colour_dict_div, colour_list, ind, patches_list, smple_id_to_plot, sp_output_df_div):
    bottom_div = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.
    for seq in list(sp_output_df_div):
        # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
        rel_abund_div = sp_output_df_div.loc[smple_id_to_plot, seq]
        if rel_abund_div > 0:
            patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            colour_list.append(colour_dict_div[seq])
            bottom_div += rel_abund_div

def plot_type_under_div(colour_dict_type, colour_list, ind, patches_list, smple_id_to_plot, sp_output_df_type):
    # the idea of the type is to put it as a reflection below the y=0 line
    # as such we should just want to make everything negative
    bottom_type = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.
    # we want to plot the rects so that they add to 1. As such we want to divide
    # each value by the total for that sample.
    tot_for_sample = sp_output_df_type.loc[smple_id_to_plot].sum()
    for its2_profile in list(sp_output_df_type):
        rel_abund = sp_output_df_type.loc[smple_id_to_plot, its2_profile]
        if rel_abund > 0:
            depth = -0.2 * (rel_abund / tot_for_sample)
            patches_list.append(
                Rectangle((ind - 0.5, bottom_type), 1, depth,
                          color=colour_dict_type[its2_profile]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            colour_list.append(colour_dict_type[its2_profile])
            bottom_type += depth

def get_div_colour_dict_and_ordered_list_of_seqs(sp_output_df_div):
    colour_palette_div = get_colour_list()
    grey_palette_div = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    # get a list of the sequences in order of their abundance and use this list to create the colour dict
    # the abundances can be got by simply summing up the columns making sure to ommit the last columns
    abundance_dict = {}
    for col in list(sp_output_df_div):
        abundance_dict[col] = sum(sp_output_df_div[col])
    # get the names of the sequences sorted according to their totalled abundance
    ordered_list_of_seqs = [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    # If we aer only going to have a legend that is cols x rows as shown below, then we should only use
    # that many colours in the plotting.
    max_n_cols = 8
    max_n_rows = 7
    num_leg_cells = max_n_cols * max_n_rows
    colour_dict_div = {}
    for i in range(len(ordered_list_of_seqs)):
        if i < num_leg_cells:
            colour_dict_div[ordered_list_of_seqs[i]] = colour_palette_div[i]
        else:
            grey_index = i % len(grey_palette_div)
            colour_dict_div[ordered_list_of_seqs[i]] = grey_palette_div[grey_index]
    return colour_dict_div, max_n_cols, max_n_rows, num_leg_cells, ordered_list_of_seqs

def process_type_df(path_to_tab_delim_count_type):
    sp_output_df_type = pd.read_csv(path_to_tab_delim_count_type, sep='\t', lineterminator='\n',
                                    skiprows=[0, 1, 2, 3, 5],
                                    header=None)
    # get a list of tups that are the seq names and the abundances zipped together
    type_profile_to_abund_tup_list = [(name, int(abund)) for name, abund in
                                      zip(sp_output_df_type.iloc[1][2:].values.tolist(),
                                          sp_output_df_type.iloc[0][2:].values.tolist())]
    # convert the names that are numbers into int strings rather than float strings.
    int_temp_list = []
    for name_abund_tup in type_profile_to_abund_tup_list:
        try:
            int_temp_list.append((str(int(name_abund_tup[0])), int(name_abund_tup[1])))
        except:
            int_temp_list.append((name_abund_tup[0], int(name_abund_tup[1])))
    type_profile_to_abund_tup_list = int_temp_list
    # need to drop the rows that contain the sequence accession and species descriptions
    for i, row_name in enumerate(sp_output_df_type.iloc[:, 0]):
        if 'Sequence accession' in row_name:
            # then we want to drop all rows from here until the end
            index_to_drop_from = i
            break
    sp_output_df_type = sp_output_df_type.iloc[:index_to_drop_from]
    # now drop the sample name columns
    sp_output_df_type.drop(columns=1, inplace=True)
    # make headers
    sp_output_df_type.columns = ['sample_id'] + [a[0] for a in type_profile_to_abund_tup_list]
    # now drop the local abund row and promote the its2_type_prof names to columns headers.
    sp_output_df_type.drop(index=[0, 1], inplace=True)
    sp_output_df_type = sp_output_df_type.set_index(keys='sample_id', drop=True).astype('float')

    max_n_cols = 4
    max_n_rows = 7
    num_leg_cells = max_n_cols * max_n_rows

    # we will use the col headers as the its2 type profile order for plotting but we
    # we should colour according to the abundance of the its2 type profiles
    # as we don't want to run out of colours by the time we get to profiles that are very abundant.
    # The sorted_type_prof_names_by_local_abund object has the names of the its2 type profile in order of abundance
    # we will use the index order as the order of samples to plot
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    sorted_type_prof_names_by_local_abund = [a[0] for a in
                                             sorted(type_profile_to_abund_tup_list, key=lambda x: x[1], reverse=True)]

    # we should plot sample by sample and its2 type by its2 type in the order of the output
    # the problem with doing he convert_to_pastel is that the colours become very similar
    # colour_palette = convert_to_pastel(get_colour_list())
    # Rather, I will attempt to generate a quick set of colours that are pastel and have a minimum distance
    # rule for any colours that are generated from each other.
    # let's do this for 50 colours to start with and see how long it takes.
    # turns out it is very quick. Easily quick enough to do dynamically.
    # When working with pastel colours (i.e. mixing with 255,255,255 it is probably best to work with a smaller dist cutoff
    if os.path.isfile('{}/colour_dict_type.pickle'.format(os.getcwd())):
        colour_dict_type = pickle.load(open('{}/colour_dict_type.pickle'.format(os.getcwd()), 'rb'))
    else:
        colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                              create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=2000, num_cols=30,
                                                 time_out_iterations=100000)]
        # todo, I want to have a constant set of pastelles here so that I can reporduce this figure and also
        # use aspects of it in other figures and have agreements. To do this I will pickle out this pastelle pallette
        # to reuse.

        # # The below 3d scatter produces a 3d scatter plot to examine the spread of the colours created
        # from mpl_toolkits.mplot3d import Axes3D
        # colour_palette = create_colour_list(sq_dist_cutoff=5000)
        # hex_pal = ['#%02x%02x%02x' % rgb_tup for rgb_tup in colour_palette]
        # colcoords = [list(a) for a in zip(*colour_palette)]
        # print(colcoords)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter(colcoords[0], colcoords[1], colcoords[2], c=hex_pal, marker='o')
        # colour_palette = get_colour_list()
        grey_palette_type = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']


        colour_dict_type = {}
        for i in range(len(sorted_type_prof_names_by_local_abund)):
            if i < num_leg_cells:
                colour_dict_type[sorted_type_prof_names_by_local_abund[i]] = colour_palette_pas[i]
            else:
                grey_index = i % len(grey_palette_type)
                colour_dict_type[sorted_type_prof_names_by_local_abund[i]] = grey_palette_type[grey_index]
        pickle.dump(colour_dict_type, open('{}/colour_dict_type.pickle'.format(os.getcwd()), 'wb'))
    return colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund, max_n_cols, max_n_rows, num_leg_cells

def process_div_df(path_to_tab_delim_count_DIV):
    sp_output_df = pd.read_csv(path_to_tab_delim_count_DIV, sep='\t', lineterminator='\n', header=0, index_col=0)

    # In order to be able to drop the DIV row at the end and the meta information rows, we should
    # drop all rows that are after the DIV column. We will pass in an index value to the .drop
    # that is called here. To do this we need to work out which index we are working with
    index_values_as_list = sp_output_df.index.values.tolist()
    for i in range(-1, -(len(index_values_as_list)), -1):
        if index_values_as_list[i].startswith('DIV'):
            # then this is the index (in negative notation) that we need to cut from
            meta_index_to_cut_from = i
            break
    sp_output_df = sp_output_df.iloc[:meta_index_to_cut_from]

    # create sample id to sample name dict
    smp_id_to_smp_name_dict = {ID: '_'.join(nm.split('_')[:3]) for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}
    smp_name_to_smp_id_dict = {'_'.join(nm.split('_')[:3]): ID for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}

    # now lets drop the QC columns from the SP output df and also drop the clade summation columns
    # we will be left with just clumns for each one of the sequences found in the samples
    sp_output_df.drop(columns=['sample_name', 'noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                               'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                               'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                               'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                               'post_taxa_id_absolute_non_symbiodinium_seqs',
                               'post_taxa_id_unique_non_symbiodinium_seqs',
                               'size_screening_violation_absolute', 'size_screening_violation_unique',
                               'post_med_absolute', 'post_med_unique'
                               ], inplace=True)
    sp_output_df = sp_output_df.astype('float')
    return smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df


# This is the method to call for plotting the QC quantitative plot
def generate_qc_summary_figure():
    info_df = generate_info_df_for_samples()
    # #  switch off the 2044 and 2041 exclusion when they have been run
    # info_df.drop(['CO0002044', 'CO0002041'], axis=0, inplace=True)

    # we want to read in the table that contains both the coral and non-coral sequences so that we can
    # plot the csw and the surface water samples along side the corals
    path_to_tab_delim_rel_count_DIV_coral_non_coral_standalone = '2018-10-21_08-59-37.620726.DIVs.absolute.txt'

    # read in the SymPortal output
    sp_output_df = pd.read_csv(path_to_tab_delim_rel_count_DIV_coral_non_coral_standalone, sep='\t', lineterminator='\n')

    sp_output_df = sp_output_df.iloc[:-4]
    # get rid of the bottom

    # generate a new index from the sample names


    new_index = [smple_name.split('_')[0] for smple_name in sp_output_df['sample_name'].values.tolist()]
    sp_output_df.index = new_index




    # The SP output contains the QC info columns between the DIVs and the no_name ITS2 columns.
    # lets put the QC info columns into a seperate df.
    QC_info_df = sp_output_df[['sample_name', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                               'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                               'post_taxa_id_absolute_non_symbiodinium_seqs',
                               'post_taxa_id_unique_non_symbiodinium_seqs',
                               'size_screening_violation_absolute', 'size_screening_violation_unique',
                               'post_med_absolute', 'post_med_unique']]

    f, axarr = plt.subplots(3, 1, figsize=(7, 7))
    # counter to reference which set of axes we are plotting on
    axarr_index = 0
    # y_axis_labels = ['raw_contigs', 'post_qc', 'Symbiodinium', 'non-Symbiodinium', 'post-MED', 'post-MED / pre-MED']
    x_axis_labels = ['raw_contigs', 'non-Symbiodiniaceae', 'Symbiodiniaceae']

    for sub_plot_type in [('raw_contigs',),
                          ('post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs'),
                          ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs')]:

        # The grids were confusing when there were two axes
        # axarr[axarr_index].grid(b=True, which='major', axis='y')
        # for each of the sub plots we will want to grab the absolute and unique counts and plot these
        # for each of the sample types.
        # go environment type by environment type

        # we will create some x axis indicies to arranage where we will be ploting
        # we can be smart with these later on and create some nice spacing layouts but for the time
        # being lets just get things plotted. Let's have one idices for each sample type and work
        # relatively from there.
        ind = range(4)
        ind_index = 0


        if sub_plot_type[0] != 'raw_contigs':
            ax2 = axarr[axarr_index].twinx()
            # ax2.set_yscale('symlog')

            # axarr[axarr_index].set_yscale('symlog')
        else:
            axarr[axarr_index].set_xlabel(x_axis_labels[axarr_index])
            # axarr[axarr_index].set_yscale('symlog')

        # we will convert the sed_close and sed_far to simply sed
        env_types_list = ['CORAL', 'CSW', 'SURFACE', 'PLANKTON']

        for env_type in env_types_list:

            if sub_plot_type[0] == 'raw_contigs':
                # here we will plot just the raw_contigs
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                if env_type == 'CORAL':
                    env_info_df = info_df[(info_df['spp_water'] == 'PORITES') | (info_df['spp_water'] == 'MILLEPORA') | (info_df['spp_water'] == 'POCILLOPORA')]
                else:
                    env_info_df = info_df[info_df['spp_water'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                # axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'PLANKTON':
                    axarr[axarr_index].spines['left'].set_color(c='blue')
                    # axarr[axarr_index].set_ylabel('total sequences', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')
                    axarr[axarr_index].spines['right'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    axarr[axarr_index].spines['top'].set_visible(False)

                    # set the ticks
                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    # set the xaxis title
                    axarr[axarr_index].set_xlabel('raw_contigs', fontsize='large')
                    axarr[axarr_index].set_ylim(10000, 100000)


                    # axarr[axarr_index].set_ylim((0, 1200000))

                ind_index += 1
            elif sub_plot_type[0] != 'med_ratio':
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                if env_type == 'CORAL':
                    env_info_df = info_df[(info_df['spp_water'] == 'PORITES') | (info_df['spp_water'] == 'MILLEPORA') | (info_df['spp_water'] == 'POCILLOPORA')]
                else:
                    env_info_df = info_df[info_df['spp_water'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                # axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                # if env_type == 'CORAL':
                #     axarr[axarr_index].set_ylabel('', color='b')
                #     axarr[axarr_index].tick_params('y', colors='b')


                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = list(env_QC_info_df.loc[:, sub_plot_type[1]])
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                # ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'PLANKTON':
                    if sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        # the legend work for the SYmbiodiniaceae
                        axarr[axarr_index].spines['left'].set_color(c='blue')
                        ax2.spines['left'].set_color(c='blue')
                        axarr[axarr_index].tick_params('y', colors='b')

                        axarr[axarr_index].set_ylim(0, 25000)
                        ax2.set_ylim(0, 400)
                        axarr[axarr_index].spines['right'].set_color(c='red')
                        ax2.spines['right'].set_color(c='red')

                        ax2.tick_params('y', colors='r')

                        axarr[axarr_index].spines['top'].set_visible(False)
                        ax2.spines['top'].set_visible(False)

                        # axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False,
                        #                                labelbottom=False)
                        # ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        axarr[2].set_xticks([0.1875, 1.1875, 2.1875, 3.1875])
                        axarr[2].set_xticklabels(env_types_list, fontsize='large')
                        axarr[axarr_index].set_xlabel(x_axis_labels[axarr_index], fontsize='large')

                    else:
                        # the legend work for non-symbiodiniaceae
                        # axarr[axarr_index].set_ylabel('total sequences', color='b')
                        axarr[axarr_index].spines['left'].set_color(c='blue')
                        ax2.spines['left'].set_color(c='blue')
                        axarr[axarr_index].tick_params('y', colors='b')

                        axarr[axarr_index].set_ylim(0, 25000)
                        ax2.set_ylim(0, 400)
                        axarr[axarr_index].spines['right'].set_color(c='red')
                        ax2.spines['right'].set_color(c='red')
                        # ax2.set_ylabel('distinct sequences', color='r')
                        ax2.tick_params('y', colors='r')

                        axarr[axarr_index].spines['top'].set_visible(False)
                        ax2.spines['top'].set_visible(False)


                        axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        axarr[axarr_index].set_xlabel(x_axis_labels[axarr_index], fontsize='large')



                ind_index += 1

        axarr_index += 1

    apples = 'asdf'
    f.text(0.02, 0.55, 'total sequences', va='center', ha='center', rotation='vertical', color='b', fontsize='large')
    f.text(1 - 0.02, 0.40, 'distinct sequences', ha='center', va='center', rotation='vertical', color='r', fontsize='large')
    # f.text(0.07, 0.18, 'ratio', va='center', rotation='vertical', color='b')
    # f.text(1 - 0.05, 0.18, 'ratio', va='center', rotation='vertical', color='r')


    # plt.tight_layout()
    f.savefig('qc_stats.svg')
    f.savefig('qc_stats.png')
    return


# These are methods used in the production of the dataframes that are used in the plotting methods
# They work on the SymPortal Outputs or the directory structure of the Genoscope ftp server
def write_out_stats_and_reorder_files():

    # The purpose of this will be to produce a couple of very basic stats of how many samples
    # we have for the different types and then to write out the files into two directories.
    # We will write out the coral samples into one directory and the water samples into another.
    # We do this as the coral samples will go all the way through the SP analysis and have predicted.
    # The water samples can just be submitted to the db and we can get the sequences through QC. They
    # will still need further processing but this will be more tertiary analysis.


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


# These are general methods used in throughout to get colours
def get_colour_list():
    colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                  "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                  "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                  "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                  "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                  "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                  "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                  "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                  "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                  "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                  "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                  "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                  "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                  "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                  "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                  "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                  "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                  "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                  "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                  "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                  "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                  "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                  "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                  "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                  "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                  "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                  "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                  "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list

def create_colour_list(sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000, avoid_black_and_white=True):
    new_colours = []
    min_dist = []
    attempt = 0
    while len(new_colours) < num_cols:
        attempt += 1
        # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
        if attempt > time_out_iterations:
            sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                     'and have not been able to find a colour that fits into your defined colour space.\n'
                     'Please lower the number of colours you are trying to find, '
                     'the minimum distance between them, or both.'.format(attempt))
        if mix_col:
            r = int((random.randint(0, 255) + mix_col[0]) /2)
            g = int((random.randint(0, 255) + mix_col[1]) /2)
            b = int((random.randint(0, 255) + mix_col[2]) /2)
        else:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        # now check to see whether the new colour is within a given distance
        # if the avoids are true also
        good_dist = True
        if sq_dist_cutoff:
            dist_list = []
            for i in range(len(new_colours)):
                distance = (new_colours[i][0] - r)**2 + (new_colours[i][1] - g)**2 + (new_colours[i][2] - b)**2
                dist_list.append(distance)
                if distance < sq_dist_cutoff:
                    good_dist = False
                    break
            # now check against black and white
            d_to_black = (r - 0)**2 + (g - 0)**2 + (b - 0)**2
            d_to_white = (r - 255)**2 + (g - 255)**2 + (b - 255)**2
            if avoid_black_and_white:
                if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                    good_dist = False
            if dist_list:
                min_dist.append(min(dist_list))
        if good_dist:
            new_colours.append((r,g,b))
            attempt = 0

    return new_colours


figure_making_bar_plots_corals()
