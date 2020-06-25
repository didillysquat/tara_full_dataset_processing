#!/usr/bin/env python3
"""
Script for plotting up the PCoAs of the distance plots and plotting up the line
graphs that show how the mantel tests change in reaction to the line graphs.

In general the order of these plots will be quite important as they are dependent on each other.
Preceeding plots will define proceeding as we will learn parameters from one that will be used in the
next. This narative should be written into the figure legends.

List of the plots that we want

PCoA plots:
1 - PCoA using the sequences as they are (i.e. with the most abundant sequences in place.)
2 - PCoA using with the majority sequences removed (but with the secondary samples still in place).
3 - PCOA using with the majority seuqences and secondary samples removed. This should likely be
plotted using the parameters inferred from the line plots.



Line plots:
1 - three row plot
Three rows to this plot that correspond to the three 'approaches' of the paper.
A - In the first row we will plot normalisation abundance against the persons correlation.
We can annotate individual point to show significance of the results. We should also likely annotate the
number of samples that are being compared.
As lines we will have each of the species, distance method and SNP w/wo = 2x2x2 = 8.
We also want to see the effect of normalisation method, but perhaps we will plot this in a separate plot.
We also want to plot the same but for minimum distinct seqs. This will probably have to go in a separate plot.

B - The second row will be similar to the first but looking at samples_at_least_threshold. So this will
be on the X.

C - The third row will again be similar format to the 2 rows above but looking at most_abund_seq_cutoff.

For rows B and C we will hold constant values according to the results of A.

As such. Let's start with the first row (A).

"""

from base_18s import EighteenSBase
import os
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from clustering_18s import ComputeClassificationAgreement
from multiprocessing import Queue, Process

class MSPlots(EighteenSBase):
    def __init__(self):
        super().__init__()

    def plot_three_row(self, classifications=False):
        """
        The three row plot.
        We are quite advanced with this figure now. In the end it looks like we may split this
        into three figures that can each likely go into the suppelementary figures.
        
        FIG ONE can be norm abundance testing and norm ethod testing.
        redo the plotting of the norm abundance using the best cutoff arguments.
        You'll need to write this in in distance.py and then write in the plotting here.
        The effect of normalisation can also be done on these parameters. These will be computed
        as a result of changing distance.py. Scale the X of these down to ~20000, as there is little change.
        see if with these new parameters being used the unifrac settles down a bit and whether it is
        worth testing a greater normalisation abundance for unifrac (hopefully not).
        
        FIG TWO
        This will be the line plots. Do this as a 2x2 with rows as the samples_at_least_threshold and 
        most_abund_seq_cutoff. Columns will be the coeficient and p_value.
        put the vlines on here that represent the maximums/areas of interest. Maybe make these a shaded vertical
        area using an alpha.

        Fig three
        The contours of coefficient for each of the species and dist methods.
        increase the values up to 300 for the most_abund_seq_cutoff (these are calculating now).
        Then make into a 2 x 2 and put a scale for each.

        then it will be time to move towards seeing the level of agreement we're getting.
        """
        
        tr = ThreeRow(parent=self)
        if not classifications:
            tr.plot()
        else:
            tr.plot_classifications()

class ThreeRow:
    def __init__(self, parent):
        self.parent = parent
        self.genera = ['Pocillopora', 'Porites']
        # Let's start with the first plot quick and dirty and then we can add the others
        # and refactorize.
        
        # plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_fig.png'), dpi=300)
        self.line_style_dict = {'rai':'-', 'pwr':'--', 'braycurtis':'-', 'unifrac':'--', True:'-', False:'--'}
        self.line_color_dict = {'Pocillopora':'black', 'Porites':'red'}

    def plot(self):
        """
        For the record, the results string format:
        self.unique_string = f'{self.genus}_{self.remove_maj_seq}_{self.exclude_secondary_seq_samples}_' \
        f'{self.exclude_no_use_samples}_{self.use_replicates}_' \
        f'{snp_distance_type}_{self.dist_method_18S}_' \
        f'{self.approach}_{self.normalisation_abundance}_{self.normalisation_method}_' \
        f'{self.only_snp_samples}_{self.samples_at_least_threshold}_' \
        f'{self.most_abund_seq_cutoff}_{self.min_num_distinct_seqs_per_sample}'
        self.result_path = os.path.join('/home/humebc/projects/tara/tara_full_dataset_processing/18s/output', f'{self.unique_string}_mantel_result.txt')
        """
        self.norm_fig, self.norm_ax = plt.subplots(1,2, figsize=(8,4))
        self.norm_ax[0].set_ylabel("Pearson's correlation coefficient")
        self.norm_ax[1].set_ylabel("Pearson's correlation coefficient")
        self._plot_first_row()
        plt.tight_layout()
        plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_norm_fig.png'), dpi=1200)

        self.cutoff_fig, self.cutoff_ax = plt.subplots(2,2, figsize=(8,8))
        self._plot_second_row()
        self._plot_third_row()
        self.cutoff_ax[0,0].set_ylabel("Pearson's correlation coefficient")
        self.cutoff_ax[1,0].set_ylabel("Pearson's correlation coefficient")
        self.cutoff_ax[0,1].set_ylabel('p-value')
        self.cutoff_ax[1,1].set_ylabel('p-value')
        self.cutoff_ax[0,0].set_xlabel('minimum in sample cutoff')
        self.cutoff_ax[1,0].set_xlabel('most abundant sequence cutoff')
        self.cutoff_ax[0,1].set_xlabel('minimum in sample cutoff')
        self.cutoff_ax[1,1].set_xlabel('most abundant sequence cutoff')
        self.cutoff_ax[0,1].legend(loc='upper left', fontsize='x-small')
        plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_cutoff_fig.png'), dpi=1200)
        
        self.contour_fig, self.contour_ax = plt.subplots(2,2, figsize=(8,8))
        self._plot_fourth_row()
        self.contour_ax[0,0].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[1,0].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[0,1].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[1,1].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[0,0].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[1,0].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[0,1].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[1,1].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        # 00
        # hlines
        self.contour_ax[0,0].plot([0,300], [0.08,0.08], 'r-', linewidth=0.5)
        # vlines
        y_max = self.contour_ax[0,0].get_ylim()[1]
        for i in range(100,300,50):
            self.contour_ax[0,0].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #01
        self.contour_ax[0,1].plot([0,300], [0.08,0.08], 'r-', linewidth=0.5)
        y_max = self.contour_ax[0,1].get_ylim()[1]
        for i in range(100,225,25):
            self.contour_ax[0,1].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #10
        self.contour_ax[1,0].plot([0,300], [.02,.02], 'r-', linewidth=0.5)
        self.contour_ax[1,0].plot([0,300], [.66,.66], 'r-', linewidth=0.5)
        y_max = self.contour_ax[1,0].get_ylim()[1]
        for i in range(25,150,25):
            self.contour_ax[1,0].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #11
        self.contour_ax[1,1].plot([75,75], [0,self.contour_ax[1,1].get_ylim()[1]], 'r-', linewidth=0.5)
        for i in np.arange(0,0.6,0.1):
            self.contour_ax[1,1].plot([0,300], [i,i], 'r-', linewidth=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_contour_fig.png'), dpi=1200)
        
        print('DONE')
        self.foo = 'bar'

    def plot_classifications(self):
        """
        For the record, the results string format:
        self.unique_string = f'{self.genus}_{self.remove_maj_seq}_{self.exclude_secondary_seq_samples}_' \
        f'{self.exclude_no_use_samples}_{self.use_replicates}_' \
        f'{snp_distance_type}_{self.dist_method_18S}_' \
        f'{self.approach}_{self.normalisation_abundance}_{self.normalisation_method}_' \
        f'{self.only_snp_samples}_{self.samples_at_least_threshold}_' \
        f'{self.most_abund_seq_cutoff}_{self.min_num_distinct_seqs_per_sample}'
        self.result_path = os.path.join('/home/humebc/projects/tara/tara_full_dataset_processing/18s/output', f'{self.unique_string}_mantel_result.txt')
        """
        # self.norm_fig, self.norm_ax = plt.subplots(1,2, figsize=(8,4))
        # self.norm_ax[0].set_ylabel("Pearson's correlation coefficient")
        # self.norm_ax[1].set_ylabel("Pearson's correlation coefficient")
        # self._plot_first_row()
        # plt.tight_layout()
        # plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_norm_fig.png'), dpi=1200)

        # self.cutoff_fig, self.cutoff_ax = plt.subplots(2,2, figsize=(8,8))
        # self._plot_second_row()
        # self._plot_third_row()
        # self.cutoff_ax[0,0].set_ylabel("Pearson's correlation coefficient")
        # self.cutoff_ax[1,0].set_ylabel("Pearson's correlation coefficient")
        # self.cutoff_ax[0,1].set_ylabel('p-value')
        # self.cutoff_ax[1,1].set_ylabel('p-value')
        # self.cutoff_ax[0,0].set_xlabel('minimum in sample cutoff')
        # self.cutoff_ax[1,0].set_xlabel('most abundant sequence cutoff')
        # self.cutoff_ax[0,1].set_xlabel('minimum in sample cutoff')
        # self.cutoff_ax[1,1].set_xlabel('most abundant sequence cutoff')
        # self.cutoff_ax[0,1].legend(loc='upper left', fontsize='x-small')
        # plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_cutoff_fig.png'), dpi=1200)
        
        self.contour_fig, self.contour_ax = plt.subplots(2,2, figsize=(8,8))
        self._plot_fourth_row_classification()
        self.contour_ax[0,0].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[1,0].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[0,1].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[1,1].set_ylabel('minimum in sample cutoff', fontsize='x-small')
        self.contour_ax[0,0].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[1,0].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[0,1].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        self.contour_ax[1,1].set_xlabel('most abundant sequence cutoff', fontsize='x-small')
        # 00
        # hlines
        self.contour_ax[0,0].plot([0,300], [0.08,0.08], 'r-', linewidth=0.5)
        # vlines
        y_max = self.contour_ax[0,0].get_ylim()[1]
        for i in range(100,300,50):
            self.contour_ax[0,0].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #01
        self.contour_ax[0,1].plot([0,300], [0.08,0.08], 'r-', linewidth=0.5)
        y_max = self.contour_ax[0,1].get_ylim()[1]
        for i in range(100,225,25):
            self.contour_ax[0,1].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #10
        self.contour_ax[1,0].plot([0,300], [.02,.02], 'r-', linewidth=0.5)
        self.contour_ax[1,0].plot([0,300], [.66,.66], 'r-', linewidth=0.5)
        y_max = self.contour_ax[1,0].get_ylim()[1]
        for i in range(25,150,25):
            self.contour_ax[1,0].plot([i,i], [0,y_max], 'r-', linewidth=0.5)

        #11
        self.contour_ax[1,1].plot([75,75], [0,self.contour_ax[1,1].get_ylim()[1]], 'r-', linewidth=0.5)
        for i in np.arange(0,0.6,0.1):
            self.contour_ax[1,1].plot([0,300], [i,i], 'r-', linewidth=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_contour_fig_classification.png'), dpi=1200)
        
        print('DONE')
        self.foo = 'bar'

    def _plot_first_row(self):
        # RESULT This plot shows us that UniFrac is bascially going mental
        # The results basically look random.
        # RESULT For Bray curtis we see a slight improvement
        # up to about the 10 000 point. So I think this looks like a sensible cutoff to work with.
        # RESULT The pwr vs rai makes very little difference so I would say that we can work with either
        # The only difference being that one comes from below and one from above. We can work with either moving forwards
        # RESULT For the with and without SNP samples there is almost no difference. I think we can put a very positive
        # spin on this. This means that the additional samples are not effecting how they are being resolved.
        # plot out an example pcoa plot of this to show that the poosition of samples does not change much
        # actually not sure that this is possible but maybe investigate.
        # As for displaying these facts, Porites and Pocillopora separate nicely as there appears to be a far stronger
        # correlation between porties (almost double). This may be due to the lack of structuring in Pocillopora.
        # This means that it works well to have both genera in each of the plots. I think we can do a plot for
        # distance method, normalisation method and for only_smp_samples.
        # given that the unifrac basically doesn't work with the 0 cutoff, there's not much point testing it using
        # this value. So Rather, we should test it using one of the cutoffs that will be determined in row 2 I guess.
        for g in self.genera:
            for m in ['unifrac', 'braycurtis']:
                self._plot_line_first_row(ax=self.norm_ax[0], genus=g, label=f'{g}_{m}', color=self.line_color_dict[g], linestyle=self.line_style_dict[m], normalisation_method='rai', distance_method=m, snp_only=False)

            for n_m in ['pwr', 'rai']:
                self._plot_line_first_row(ax=self.norm_ax[1], genus=g, label=f'{g}_{n_m}', color=self.line_color_dict[g], linestyle=self.line_style_dict[n_m],normalisation_method=n_m, distance_method='braycurtis', snp_only=False)

        for ax in self.norm_ax:
            ax.set_xlabel('normalisation abundance')
            ax.set_xlim(0,10000)
            ax.set_ylim(0,0.25)
            ax.legend(loc='lower right', fontsize='x-small')
            # NB the only way that having SNP samples can make a difference is during the Unifrac and Unifrac isn't really working here
            # Otherwise, if you think about it, exactly the same pairwise comparisons are going to be considered for the braycurtis and
            # then we're already stripping down to only the SNP comparison in the mantel.
            # We want to be looking at the SNP/noSNP for the cluster assignment.
            # for only_snp_samples in [True, False]: # We are computing this now.
            #     self._plot_line_first_row(ax=self.ax[0][2], genus=g, color=self.line_color_dict[g], linestyle=self.line_style_dict[only_snp_samples],normalisation_method='rai', distance_method='unifrac', snp_only=only_snp_samples)
        
        self.norm_ax[0].set_title('distance_method')
        self.norm_ax[1].set_title('normalisation_method')

    def _plot_second_row(self):
        # RESULT This shows us that the effect of the samples_at_least_threshold is dependent on the genus
        # being investigated. For Porites it has a negative effect.
        # However, for Pocillopora, it hav a positive effect.
        # The distance method i.e. unifrac or braycurtis has little effect.
        for g in self.genera:
            for m in ['unifrac', 'braycurtis']:
                if m == 'unifrac':
                    norm_abund = 1000
                else:
                    norm_abund = 10000
                self._plot_line_second_row(ax=self.cutoff_ax[0][0], plot_type='coef', genus=g, label=f'{g}_{m}', color=self.line_color_dict[g], normalisation_abundance=norm_abund, linestyle=self.line_style_dict[m], normalisation_method='pwr', distance_method=m, snp_only=False)
                self._plot_line_second_row(ax=self.cutoff_ax[0][1], plot_type='p_val', genus=g, label=f'{g}_{m}', color=self.line_color_dict[g], normalisation_abundance=norm_abund, linestyle=self.line_style_dict[m], normalisation_method='pwr', distance_method=m, snp_only=False)

    def _plot_third_row(self):
        # RESULTS This shows us that once again, the effect is very genus dependent
        # For Porites, using this threshold argubly has some benefit to a small degree but questionable.
        # However for Pocillopora it appears to have little or no effect.
        # test some combinations of these two factors to see if we find some surprising results.
        for g in self.genera:
            for m in ['unifrac', 'braycurtis']:
                if m == 'unifrac':
                    norm_abund = 1000
                else:
                    norm_abund = 10000
                self._plot_line_third_row(ax=self.cutoff_ax[1][0], genus=g, plot_type='coef', color=self.line_color_dict[g], normalisation_abundance=norm_abund, linestyle=self.line_style_dict[m], normalisation_method='pwr', distance_method=m, snp_only=False)
                self._plot_line_third_row(ax=self.cutoff_ax[1][1], genus=g, plot_type='p_val', color=self.line_color_dict[g], normalisation_abundance=norm_abund, linestyle=self.line_style_dict[m], normalisation_method='pwr', distance_method=m, snp_only=False)
     
    def _plot_fourth_row(self):
        # RESULT. This has worked pretty well. It shows us that we can basically take a linear combination
        # of the results from the individual testing of the samples_at_least_threshold and the most_abund_seq_cutoff
        # and predict the result with the combination. Similar to the individual results it shows us that
        # the approach needs to be species specific. Pocillopora shoud have a samples_at_least_threshold of
        # approximately 0.1 (we should draw this onto the plot as a vline on the individual and contour.
        # And that the Porites should use a most_abundant_seq_cutoff of approximately 60 or 70 (draw on as vline).
        # We should also extend the values for which we have results in the contours upto 300 so that we can compare
        # better to the individual plots. scale the individual plots down to 300.
        # We can additionally now in theory work on the crosses of the samples_at_least_threshold and the 
        # most_abund_seq_cutoff.
        # We will need to look at a countour for each variable combination I guess
        # so one for each of the genera and dist methods
        for g_i, g in enumerate(self.genera):
            for m_i, m in enumerate(['unifrac', 'braycurtis']):
                if m == 'unifrac':
                    norm_abund = 1000
                else:
                    norm_abund = 10000
                contour = self._plot_countour(ax=self.contour_ax[g_i,m_i], genus=g, normalisation_abundance=norm_abund, normalisation_method='pwr', distance_method=m, snp_only=False)
                self.contour_ax[g_i,m_i].set_title(f'{g}_{m}', fontsize='x-small')
                cbar = self.contour_fig.colorbar(contour, ax=self.contour_ax[g_i,m_i])
                cbar.ax.set_ylabel("Pearson's correlation coefficient")

    def _plot_fourth_row_classification(self):
        """
        A variation on the _plot_fourth_row. This will look at classification agreement rather than mantel based Pearson's
        """
        for g_i, g in enumerate(self.genera):
            for m_i, m in enumerate(['unifrac', 'braycurtis']):
                for i_i, i in enumerate(['original_three', 'ten_plus_one']):
                    if m == 'unifrac':
                        norm_abund = 1000
                    else:
                        norm_abund = 10000
                    # If unifrac need to incorporate the island_list string diretly when looking for the distance matrix.
                    # If bray curtis then we will work with a single base distance matrix and remove samples accordding to which
                    # island list we will be working with.
                    contour = self._plot_countour_classification(ax=self.contour_ax[g_i,m_i], genus=g, normalisation_abundance=norm_abund, normalisation_method='pwr', distance_method=m, snp_only=False, island_list=i)
                    self.contour_ax[g_i,m_i].set_title(f'{g}_{m}', fontsize='x-small')
                    cbar = self.contour_fig.colorbar(contour, ax=self.contour_ax[g_i,m_i])
                    cbar.ax.set_ylabel("Pearson's correlation coefficient")
       
    def _plot_line_first_row(self, ax, genus, linestyle, color, label, normalisation_method='rai', distance_method='braycurtis', snp_only=False):
        
        # We are testing the normalisation values using optimised parameters that are genus dependent.
        # As such, the file we search for needs to be adjusted per genus
        if genus == 'Pocillopora':
            samples_at_least_threshold = str(0.08)
            most_abund_seq_cutoff = str(150)
        else: # genus == 'Porites':
            samples_at_least_threshold = str(0.02)
            most_abund_seq_cutoff = str(75)
        results_dict = {}
        for result_file in [_ for _ in os.listdir(self.parent.output_dir_18s) if _.endswith('_mantel_result.txt')]:
            if result_file.startswith(f'{genus}_True_True_True_False_biallelic_{distance_method}_dist'):
                # inbetween these to conditions is the nomalisation_abundance
                if result_file.endswith(f'{normalisation_method}_{snp_only}_{samples_at_least_threshold}_{most_abund_seq_cutoff}_3_mantel_result.txt'):
                # if result_file.endswith(f'{n_m}_False_0_0_3_mantel_result.txt'):
                    # Then this is a set of points for plotting
                    # We want to get them in order
                    # the normalisation depth is the 8th item
                    # create a dict of normalisation depth to tuple.
                    # Where tuple is p_value and correlation coef.
                    with open(os.path.join(self.parent.output_dir_18s, result_file), 'r') as f:
                            (cor_coef, p_val) = [float(_) for _ in f.read().rstrip().lstrip().split('\t')]
                    norm_value = int(result_file.split('_')[8])
                    if norm_value not in results_dict:
                        results_dict[norm_value] = (cor_coef, p_val)
                    else:
                        raise RuntimeError('Dict already contains norm_value')
        
        # Here we have a results dict ready for plotting for one of the g/m/snp combos.
        # For starters plot up the line plain. Then changge line characters. then annotate with p value
        sorted_norm = sorted(results_dict.keys())
        ax.plot(
            sorted_norm, 
            [results_dict[_][0] for _ in sorted_norm],
            label=f'{label}', color=color, linestyle=linestyle, linewidth=0.5)
        # ax.legend(fontsize=2, loc='best')
        foo = 'bar'

    def _plot_line_second_row(self, ax, plot_type, genus, label, linestyle, color, normalisation_abundance, normalisation_method='pwr', distance_method='braycurtis', snp_only=False):
        # We want to find all files in the 18s output directory where:
        # genus == g
        results_dict = {}
        for result_file in [_ for _ in os.listdir(self.parent.output_dir_18s) if _.endswith('_mantel_result.txt')]:
            if result_file.startswith(f'{genus}_True_True_True_False_biallelic_{distance_method}_dist_{normalisation_abundance}_{normalisation_method}_{snp_only}'):
                # inbetween these to conditions is the nomalisation_abundance
                if result_file.endswith(f'_0_3_mantel_result.txt'):
                    # Then this is a set of points for plotting
                    # We want to get them in order
                    # the samples_at_least_threshold is the 8th item
                    # create a dict of samples_at_least_threshold to tuple.
                    # Where tuple is p_value and correlation coef.
                    with open(os.path.join(self.parent.output_dir_18s, result_file), 'r') as f:
                            (cor_coef, p_val) = [float(_) for _ in f.read().rstrip().lstrip().split('\t')]
                    norm_value = float(result_file.split('_')[11])
                    if norm_value not in results_dict:
                        results_dict[norm_value] = (cor_coef, p_val)
                    else:
                        if norm_value != 0:
                            raise RuntimeError('Dict already contains norm_value')
                        else:
                            # We have _0_ and _0.0_ because of how the distances were calculated.
                            # They are actually different but both give p_val >> 0.05 and coef ~ 0.
                            pass
        
        # Here we have a results dict ready for plotting for one of the g/m/snp combos.
        # For starters plot up the line plain. Then changge line characters. then annotate with p value
        sorted_norm = sorted(results_dict.keys())
        # ax2=ax.twinx()
        if plot_type == 'coef':
            ax.plot(
                sorted_norm, 
                [results_dict[_][0] for _ in sorted_norm],
                label=f'{label}', 
                color=color,
                linestyle=linestyle,
                linewidth=0.5)
        else:
            ax.plot(
                sorted_norm, 
                [results_dict[_][1] for _ in sorted_norm],
                label=f'{label}', 
                color=color,
                linestyle=linestyle,
                linewidth=0.5)
        # ax.legend(fontsize=2, loc='best')
        foo = 'bar'

    def _plot_line_third_row(self, ax, plot_type, genus, linestyle, color, normalisation_abundance, normalisation_method='pwr', distance_method='braycurtis', snp_only=False):
        results_dict = {}
        for result_file in [_ for _ in os.listdir(self.parent.output_dir_18s) if _.endswith('_mantel_result.txt')]:
            if result_file.startswith(f'{genus}_True_True_True_False_biallelic_{distance_method}_dist_{normalisation_abundance}_{normalisation_method}_{snp_only}_0_'):
                # inbetween these to conditions is the nomalisation_abundance
                if result_file.endswith(f'_3_mantel_result.txt'):
                    # Then this is a set of points for plotting
                    # We want to get them in order
                    # the samples_at_least_threshold is the 8th item
                    # create a dict of samples_at_least_threshold to tuple.
                    # Where tuple is p_value and correlation coef.
                    with open(os.path.join(self.parent.output_dir_18s, result_file), 'r') as f:
                            (cor_coef, p_val) = [float(_) for _ in f.read().rstrip().lstrip().split('\t')]
                    norm_value = int(result_file.split('_')[12])
                    if norm_value not in results_dict:
                        results_dict[norm_value] = (cor_coef, p_val)
                    else:
                        if norm_value != 0:
                            raise RuntimeError('Dict already contains norm_value')
                        else:
                            # We have _0_ and _0.0_ because of how the distances were calculated.
                            # They are actually different but both give p_val >> 0.05 and coef ~ 0.
                            pass
        
        # Here we have a results dict ready for plotting for one of the g/m/snp combos.
        # For starters plot up the line plain. Then changge line characters. then annotate with p value
        sorted_norm = sorted(results_dict.keys())
        # ax2=ax.twinx()
        if plot_type == 'coef':
            ax.plot(
                sorted_norm, 
                [results_dict[_][0] for _ in sorted_norm],
                label=f'{genus}_{normalisation_method}_{distance_method}_{snp_only}', 
                color=color,
                linestyle=linestyle,
                linewidth=0.5)
        else:
            ax.plot(
                sorted_norm, 
                [results_dict[_][1] for _ in sorted_norm],
                label=f'{genus}_{normalisation_method}_{distance_method}_{snp_only}', 
                color=color,
                linestyle=linestyle,
                linewidth=0.5)
        # ax.legend(fontsize=2, loc='best')
        foo = 'bar'

    def _plot_countour(self, ax, genus, distance_method, normalisation_abundance, normalisation_method='pwr', snp_only=False):
        # Plot a contour plot where we have samples_at_least_threshold on the x and most_abund_seq_cutoff on the y
        # and then the coef on the z.
        x_samples_at_least_threshold = []
        y_most_abund_seq_cutoff = []
        z_coef = []
        z_p_val = []
        for result_file in [_ for _ in os.listdir(self.parent.output_dir_18s) if _.endswith('_mantel_result.txt')]:
            if result_file.startswith(f'{genus}_True_True_True_False_biallelic_{distance_method}_dist_{normalisation_abundance}_{normalisation_method}_{snp_only}_'):
                # inbetween these to conditions are the misco and the masco scores
                if result_file.endswith(f'_3_mantel_result.txt'):
                    # Then this is a set of points for plotting
                    
                    with open(os.path.join(self.parent.output_dir_18s, result_file), 'r') as f:
                            (cor_coef, p_val) = [float(_) for _ in f.read().rstrip().lstrip().split('\t')]
                    samples_at_least_threshold = float(result_file.split('_')[11])
                    most_abund_seq_cutoff = int(result_file.split('_')[12])
                    # We need to through out the 75 value most_abund_seq_cutoff as we only have this for a single
                    # samples_at_least_threshold due to the normalistaion testing.
                    # If we leave this 75 in it creates a vertical white stripe at 75.
                    if most_abund_seq_cutoff == 75:
                        continue
                    x_samples_at_least_threshold.append(samples_at_least_threshold)
                    y_most_abund_seq_cutoff.append(most_abund_seq_cutoff)
                    z_coef.append(cor_coef)
                    z_p_val.append(p_val)
                    
        df = pd.DataFrame(columns=sorted([int(_) for _ in set(y_most_abund_seq_cutoff)]), index=sorted(list(set(x_samples_at_least_threshold))))
        
        for x,y,z in zip(x_samples_at_least_threshold, y_most_abund_seq_cutoff, z_coef):
            df.at[x,y] = z
        
        # the samples_at_least_threshold has most_abund_seq_cutoff up to 580, but the combinations
        # were only computed up to 100 so crop to this
        df = df.iloc[:,:df.columns.get_loc(300) + 1]
        contour = ax.contourf(list(df), list(df.index), df.to_numpy())
        return contour

    def _plot_countour_classification(self, ax, genus, distance_method, normalisation_abundance, island_list, normalisation_method='pwr', snp_only=False, num_proc=100):
        """
        This is a modification for _plot_contour to work with classification agreement rather than Pearsons mantel agreement
        """
        # Plot a contour plot where we have samples_at_least_threshold on the x and most_abund_seq_cutoff on the y
        # and then the agreement on the z.
        x_samples_at_least_threshold = []
        y_most_abund_seq_cutoff = []
        z_coef = []
        z_p_val = []
        in_q = Queue()
        out_q = Queue()
        in_q_len = 0
        for dist_file in [_ for _ in os.listdir(self.parent.output_dir_18s) if _.endswith('.dist.gz')]:
            if dist_file.startswith(f'{genus}_True_True_True_False_biallelic_{distance_method}_dist_{normalisation_abundance}_{normalisation_method}_{snp_only}_'):
                # inbetween these two conditions are the misco and the masco scores
                if distance_method == 'unifrac':
                    # Then we need to search for the island_list in the dist name
                    if dist_file.endswith(f'_3_{island_list}.dist.gz'):
                        # Then this is a set of points for plotting
                        # Then this is a distance matrix that we want to check the classification agreement for
                        results_file  = dist_file.replace('.dist.gz', '_classification_result.txt')
                        results_path = os.path.join(self.parent.output_dir_18s, results_file)
                        #TODO Implement multiprocessing here
                        # If results path exists, then add the results to the x,y and z
                        # else, add it to the input queue for multiprocessing.
                        samples_at_least_threshold = float(dist_file.split('_')[11])
                        most_abund_seq_cutoff = int(dist_file.split('_')[12])
                        if most_abund_seq_cutoff == 75:
                                continue
                        if os.path.exists(results_path):
                            # Then we already have the results computed and we can simply read them in and plot them up
                            # We have not completed this yet
                            with open(results_path, 'r') as f:
                                max_agreement = float(f.readline().split(',')[0])
                            # We need to throw out the 75 value most_abund_seq_cutoff as we only have this for a single
                            # samples_at_least_threshold due to the normalistaion testing.
                            # If we leave this 75 in it creates a vertical white stripe at 75.
                            x_samples_at_least_threshold.append(samples_at_least_threshold)
                            y_most_abund_seq_cutoff.append(most_abund_seq_cutoff)
                            z_coef.append(max_agreement)
                        else:
                            # Then we need to compute the classificaiton agreement
                            # We can make a call here to a class of clustering.
                            distance_matrix_path = os.path.join(self.parent.output_dir_18s, dist_file)
                            in_q.put(
                                (distance_method, distance_matrix_path, island_list, genus, 
                                samples_at_least_threshold, most_abund_seq_cutoff, 
                                self.parent.output_dir, self.parent.cache_dir_18s, self.parent.input_dir_18s, 
                                self.parent.output_dir_18s, self.parent.fastq_info_df_path, self.parent.temp_dir_18s))
                            in_q_len += 1
                elif distance_method == 'braycurtis':
                    if dist_file.endswith(f'_3.dist.gz'):
                        # Then this is a set of points for plotting
                        # Then this is a distance matrix that we want to check the classification agreement for
                        # For the bray curtis we'll need to now include the island_list into the string
                        results_file  = dist_file.replace('.dist.gz', '_{island_list}_classification_result.txt')
                        results_path = os.path.join(self.parent.output_dir_18s, result_file)
                        #TODO Implement multiprocessing here
                        # If results path exists, then add the results to the x,y and z
                        # else, add it to the input queue for multiprocessing.
                        samples_at_least_threshold = float(dist_file.split('_')[11])
                        most_abund_seq_cutoff = int(dist_file.split('_')[12])
                        if most_abund_seq_cutoff == 75:
                                continue
                        if os.path.exists(results_path):
                            # Then we already have the results computed and we can simply read them in and plot them up
                            # We have not completed this yet
                            with open(results_path, 'r') as f:
                                max_agreement = float(f.readline().split(',')[0])
                            # We need to throw out the 75 value most_abund_seq_cutoff as we only have this for a single
                            # samples_at_least_threshold due to the normalistaion testing.
                            # If we leave this 75 in it creates a vertical white stripe at 75.
                            x_samples_at_least_threshold.append(samples_at_least_threshold)
                            y_most_abund_seq_cutoff.append(most_abund_seq_cutoff)
                            z_coef.append(max_agreement)
                        else:
                            # Then we need to compute the classificaiton agreement
                            # We can make a call here to a class of clustering.
                            distance_matrix_path = os.path.join(self.parent.output_dir_18s, dist_file)
                            in_q.put(
                                (distance_method, distance_matrix_path, island_list, genus, 
                                samples_at_least_threshold, most_abund_seq_cutoff, 
                                self.parent.output_dir, self.parent.cache_dir_18s, self.parent.input_dir_18s, 
                                self.parent.output_dir_18s, self.parent.fastq_info_df_path, self.parent.temp_dir_18s))
                            in_q_len += 1
                else:
                    raise NotImplementedError          

        # Here we have the x, y and z populated with the results that have already been calculated and we have
        # an MP in_q populated with tuples that represent results that still need to be processed.
        # We need to MP process them here.
        for n in range(num_proc):
            in_q.put('STOP')
        
        all_processes = []
        for p in range(num_proc):
            p = Process(target=self.execute_classification_agreement, args=(in_q, out_q))
            all_processes.append(p)
            p.start()

        # Now process the out_q in the main thread
        done_count = 0
        prog_count = 0
        while done_count < num_proc:
            item = out_q.get()
            if item == 'ALL_DONE':
                done_count += 1
            else:
                prog_count += 1
                print(f'{prog_count}/{in_q_len} dist analyses complete')
                max_agreement, misco, masco = item
                x_samples_at_least_threshold.append(misco)
                y_most_abund_seq_cutoff.append(masco)
                z_coef.append(max_agreement)

        for p in all_processes:
            p.join()

        df = pd.DataFrame(columns=sorted([int(_) for _ in set(y_most_abund_seq_cutoff)]), index=sorted(list(set(x_samples_at_least_threshold))))
        
        for x,y,z in zip(x_samples_at_least_threshold, y_most_abund_seq_cutoff, z_coef):
            df.at[x,y] = z
        
        # the samples_at_least_threshold has most_abund_seq_cutoff up to 580, but the combinations
        # were only computed up to 100 so crop to this
        df = df.iloc[:,:df.columns.get_loc(300) + 1]
        contour = ax.contourf(list(df), list(df.index), df.to_numpy())
        return contour

    @staticmethod
    def execute_classification_agreement(in_q, out_q):
        for class_args in iter(in_q.get, 'STOP'):
            (distance_method, distance_matrix_path, island_list, genus, 
            samples_at_least_threshold, most_abund_seq_cutoff, output_dir, 
            cache_dir_18s, input_dir_18s, output_dir_18s, 
            fastq_info_df_path, temp_dir_18s) = class_args
            
            max_agreement = ComputeClassificationAgreement(
                distance_method=distance_method, distance_matrix_path=distance_matrix_path, 
                island_list=island_list, genus=genus,
                output_dir=output_dir, cache_dir_18s=cache_dir_18s, input_dir_18s=input_dir_18s,
                output_dir_18s=output_dir_18s, fastq_info_df_path=fastq_info_df_path, temp_dir_18s=temp_dir_18s).compute_classficiation_agreement()

            # Once the max_agreement calculation is complete we will want to put a tuple into the out_q
            # that represents the max_agreement, the misco and the masco scores, so that these can be added
            # in the main thread to the plotting coordinates list
            out_q.put((max_agreement, samples_at_least_threshold, most_abund_seq_cutoff))
            print(f'max_agreement={max_agreement:.2f} for misco of {samples_at_least_threshold} and masco of {most_abund_seq_cutoff}')
        out_q.put('ALL_DONE')
        

# MSPlots().plot_three_row()
MSPlots().plot_three_row(classifications=True)