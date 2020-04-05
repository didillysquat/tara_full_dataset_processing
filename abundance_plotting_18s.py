"""We need to do some investigation into the absolute number of sequences and unique number of sequences
found in samples to better quality control the distance metrics"""

from base_18s import EighteenSBase
import compress_pickle
import os
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter, defaultdict
import sys

class EighteenSAbundancePlotting(EighteenSBase):
    def __init__(self):
        super().__init__()
        self._add_additional_info_to_info_df()
        self.absolute_consolidated_abundance_dict = compress_pickle.load(os.path.join(self.cache_dir, 'consolidated_df_dict_output_tables.p.bz'))
        self.fig, self.ax = plt.subplots(6,1)

    def _add_additional_info_to_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz')):
            self.info_df = compress_pickle.load(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))
        else:
            print('updating info_df')
            most_abund_coral_genus_df_list = []
            most_abund_seq_of_coral_genus_df_list = []
            for sample_name in self.info_df.index:
                sys.stdout.write(f'\r{sample_name}')
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
                coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
                most_abund_coral_genus_df_list.append(self._identify_most_abund_coral_genus(rel_all_seq_abundance_dict, coral_annotation_dict))
                consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                most_abund_seq_of_coral_genus_df_list.append(sorted([_ for _ in consolidated_host_seqs_abund_dict.items()], key=lambda x: x[1], reverse=True)[0][0])
            self.info_df['most_abund_coral_genus'] = most_abund_coral_genus_df_list
            self.info_df['most_abund_seq_of_coral_genus'] = most_abund_seq_of_coral_genus_df_list
            compress_pickle.dump(self.info_df, os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))
            print()

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

    def _get_hard_and_rel_sub_dicts(self, sample_names):
        if os.path.isfile(os.path.join(self.cache_dir, 'hard_sub_sample_dict.p.bz')):
            if os.path.isfile(os.path.join(self.cache_dir, 'rel_sub_sample_dict.p.bz')):
                return compress_pickle.load(os.path.join(self.cache_dir, 'hard_sub_sample_dict.p.bz')), compress_pickle.load(os.path.join(self.cache_dir, 'rel_sub_sample_dict.p.bz'))

        hard_sub_sample_dict = {}
        rel_sub_sample_dict = {}
        count = 0
        tot_samples = len(sample_names)
        for sample_name in sample_names:
            count += 1
            sys.stdout.write(f'\r{sample_name}: {count}/{tot_samples}')
            abund_list = self.absolute_consolidated_abundance_dict[sample_name]
            if sum(abund_list) < 10000:
                continue

            # Make a redundant list of the seqs
            non_z = []
            for i, abund in enumerate(abund_list):
                if abund > 0:
                    non_z.append(i)

            redundant_list = []
            # prob_list = []
            tot = sum(abund_list)
            for i in non_z:
                seq = self.ordered_seq_names[i]
                abund = abund_list[i]
                # prob = abund/tot
                redundant_list.extend([seq for _ in range(abund)])
                # prob_list.extend([prob for _ in range(abund)])

            hard_sub_sample_list = np.random.choice(redundant_list, 10000, replace=False)
            hard_abunds_dict = dict(Counter(hard_sub_sample_list))
            hard_sub_sample_dict[sample_name] = hard_abunds_dict

            # For soft
            norm_abund_dict = {self.ordered_seq_names[i]: int((abund_list[i] / tot) * 100) for i in non_z if int((abund_list[i] / tot) * 10000) > 0}
            rel_sub_sample_dict[sample_name] = norm_abund_dict

        compress_pickle.dump(hard_sub_sample_dict, os.path.join(self.cache_dir, 'hard_sub_sample_dict.p.bz'))
        compress_pickle.dump(rel_sub_sample_dict, os.path.join(self.cache_dir, 'rel_sub_sample_dict.p.bz'))
        return hard_sub_sample_dict, rel_sub_sample_dict

    def plot_master_scatter(self):
        sample_names = list(self.absolute_consolidated_abundance_dict.keys())
        host_only_master_seq_info_dict = compress_pickle.load(os.path.join(self.cache_dir, 'host_only_master_seq_info_dict.p.bz'))
        self.ordered_seq_names = [tup[0] for tup in sorted([_ for _ in host_only_master_seq_info_dict.items()],
                                  key=lambda x: x[1],
                                  reverse=True)]

        hard_sub_sample_dict, rel_sub_sample_dict = self._get_hard_and_rel_sub_dicts(sample_names)
        plot_type = 'seq_loss'
        if plot_type == 'hard_soft_ab':
            # Total number of seqs
            x = [sum(self.absolute_consolidated_abundance_dict[k]) for k in sample_names]
            # Unique number of seqs
            y = [len([_ for _ in v if _ > 0]) for v in [self.absolute_consolidated_abundance_dict[k] for k in sample_names]]
            self.ax[0].scatter(x=x, y=y, c='black', s=2)
            self.ax[0].set_xlim((self.ax[0].get_xlim()[0], 500000))

            hard_sample_names = hard_sub_sample_dict.keys()
            # Hard x total
            x_hard = [sum(hard_sub_sample_dict[k]) for k in hard_sample_names]
            # Hard y unique
            y_hard = [len([_ for _ in v if _ > 0]) for v in [hard_sub_sample_dict[k] for k in hard_sample_names]]

            # Rel x total
            x_rel = [sum(rel_sub_sample_dict[k]) for k in hard_sample_names]
            # Rel y unique
            y_rel = [len(rel_sub_sample_dict[k]) for k in hard_sample_names]

            # Plot hard
            self.ax[1].scatter(x=list(range(len(y_hard))), y=y_hard, c='black', s=2)
            # Plor rel
            self.ax[2].scatter(x=list(range(len(y_rel))), y=y_rel, c='black', s=2)
            foo = 'bar'
        elif plot_type == 'seq_loss':
            # Here I think it will be useful to plot the number of unique sequences that are lost
            # due to the rel_abund cutoff
            try:
                start_seqs = compress_pickle.load(os.path.join(self.cache_dir, 'start_seqs.p.bz'))
                norm_seqs = compress_pickle.load(os.path.join(self.cache_dir, 'norm_seqs.p.bz'))
                no_norm_seqs_50 = compress_pickle.load(os.path.join(self.cache_dir, 'no_norm_seqs_50.p.bz'))
                no_norm_seqs_100 = compress_pickle.load(os.path.join(self.cache_dir, 'no_norm_seqs_100.p.bz'))
                samples_to_plot = compress_pickle.load(os.path.join(self.cache_dir, 'samples_to_plot.p.bz'))
                seq_to_total_abund_dict = compress_pickle.load(
                    os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
            except:
                start_seqs = []
                norm_seqs = []
                no_norm_seqs_50 = []
                no_norm_seqs_100 = []
                samples_to_plot = []
                seq_to_total_abund_dict = compress_pickle.load(
                    os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
                counter = 0
                total_sample_name = len(sample_names)
                seq_50 = [seq for seq, val in seq_to_total_abund_dict.items() if val > 50]
                seq_100 = [seq for seq, val in seq_to_total_abund_dict.items() if val > 100]
                for sample_name in sample_names:
                    counter += 1
                    sys.stdout.write(f'\r{sample_name} {counter}/{total_sample_name}')
                    abs_abunds = self.absolute_consolidated_abundance_dict[sample_name]
                    non_z = [i for i, val in enumerate(abs_abunds) if val > 0]
                    start_seqs_sample = [self.ordered_seq_names[i] for i in non_z]

                    try:
                        norm_seqs.append(list(rel_sub_sample_dict[sample_name].keys()))
                    except KeyError:
                        continue
                    start_seqs.append(start_seqs_sample)
                    samples_to_plot.append(sample_name)
                    no_norm_seqs_50.append([seq for seq in start_seqs_sample if seq in seq_50])
                    no_norm_seqs_100.append([seq for seq in start_seqs_sample if seq in seq_100])


                compress_pickle.dump(start_seqs, os.path.join(self.cache_dir, 'start_seqs.p.bz'))
                compress_pickle.dump(norm_seqs, os.path.join(self.cache_dir, 'norm_seqs.p.bz'))
                compress_pickle.dump(no_norm_seqs_50, os.path.join(self.cache_dir, 'no_norm_seqs_50.p.bz'))
                compress_pickle.dump(no_norm_seqs_100, os.path.join(self.cache_dir, 'no_norm_seqs_100.p.bz'))
                compress_pickle.dump(samples_to_plot, os.path.join(self.cache_dir, 'samples_to_plot.p.bz'))

            # Make an island color dict
            color_list = get_colour_list()
            island_c_dict = {
                island: color for island, color in zip(
                self.info_df['island'].unique(),
                color_list[:len(self.info_df['island'].unique())])}
            c_dict_island = {
                sample_name: island_c_dict[self.info_df.at[sample_name, 'island']] for
                sample_name in self.info_df.index
            }
            # Create a color dictionary that is binary for primary or secondary sequence sample
            #
            maj_seq_c_dict = {maj_seq: color for maj_seq, color in
                              zip(self.info_df['most_abund_seq_of_coral_genus'].unique(),
                                  color_list[:len(self.info_df['most_abund_seq_of_coral_genus'].unique())])}
            c_dict_maj_seq = {sample_name: maj_seq_c_dict[self.info_df.at[sample_name, 'most_abund_seq_of_coral_genus']]
                    for sample_name in
                    self.info_df.index}
            start_abunds = [len(_) for _ in start_seqs]
            lost_through_norm = [len(pre) - len(post) for pre, post in zip(start_seqs, norm_seqs)]
            lost_through_threshold_50 = [len(pre) - len(thresh) for pre, thresh in zip(start_seqs, no_norm_seqs_50)]
            lost_through_threshold_100 = [len(pre) - len(thresh) for pre, thresh in zip(start_seqs, no_norm_seqs_100)]
            self.ax[0].scatter(x=start_abunds, y=lost_through_norm,
                               c=[c_dict_island[sample_name] for sample_name in samples_to_plot], s=2)
            self.ax[1].scatter(x=start_abunds, y=lost_through_threshold_50,
                               c=[c_dict_island[sample_name] for sample_name in samples_to_plot], s=2)
            self.ax[2].scatter(x=start_abunds, y=lost_through_threshold_100,
                               c=[c_dict_island[sample_name] for sample_name in samples_to_plot], s=2)
            #  plot up the absolute abundance against the unique seqs
            self.ax[3].scatter(x=[sum(self.absolute_consolidated_abundance_dict[sample_name]) for sample_name in samples_to_plot], y=start_abunds,
                               c=[c_dict_island[sample_name] for sample_name in samples_to_plot], s=2)
            # # I want to know which islands these pink and purple ones are
            # count_dict = defaultdict(int)
            # count_dict_high = defaultdict(int)
            # for i, sample_name in enumerate(samples_to_plot):
            #     if start_abunds[i] < 100:
            #         island = self.info_df.at[sample_name, 'island']
            #         abs_abund = self.absolute_consolidated_abundance_dict[sample_name]
            #         abs_abund_minor = sum(abs_abund) - max(abs_abund)
            #         print(f'{sample_name}: {island} {start_abunds[i]} {abs_abund_minor}')
            #         count_dict[island] += 1
            #     if start_abunds[i] > 575:
            #         island = self.info_df.at[sample_name, 'island']
            #         count_dict_high[island] += 1
            # print(count_dict)
            # print(count_dict_high)

            # Plot the unique number of sequence

            # Plot the number of sequences remaining with the threshold
            cutoff = 0
            cutoff_list = [_ for _ in seq_to_total_abund_dict.items() if _[1] > cutoff]
            x_thresh = []
            y_thresh = []
            print()
            while cutoff_list:
                sys.stdout.write(f'\r{cutoff}')
                x_thresh.append(cutoff)
                y_thresh.append(len(cutoff_list))
                cutoff += 1
                cutoff_list = [_ for _ in seq_to_total_abund_dict.items() if _[1] > cutoff]

            self.ax[4].scatter(x=x_thresh, y=y_thresh, c='black', s=2)
            print()

            # # Plot the number of total sequnces left in host samples after the most abundant seq removed
            # remaining = []
            # for v in self.absolute_consolidated_abundance_dict.values():
            #     remaining.append(sum(v)-max(v))
            # self.ax[5].scatter(x=range(len(self.absolute_consolidated_abundance_dict)), y=remaining, c='black', s=2)
            # self.ax[5].set_ylim(0,10000)


            # Want to know the number of unique sequencs that are
            foo = 'bar'

            # I want a historgram of the number of background seqs in samples
            hist_list = []
            low_hist_list = []
            high_hist_list = []
            for i, sample_name in enumerate(samples_to_plot):
                hist_list.append(start_abunds[i])
                island = self.info_df.at[sample_name, 'island']
                if island in ['I11', 'I08']:
                    low_hist_list.append(start_abunds[i])
                elif island == 'I15':
                    high_hist_list.append(start_abunds[i])

            # This shows us that there are esentially two distributions and I bet that the first lower
            # abund distribution is the low sequencing islands Islands 8 and 11. I also bet that the highstuff
            # is island 15. I will test this by plotting up histograms of this color different.
            # Yeah, so this is basically correct. The vast majority of samples have more than 100 sequences.

            fig, ax = plt.subplots(1, 1)
            ax.hist(hist_list, bins=50)
            ax.hist(low_hist_list, bins=10, color='black')
            ax.hist(high_hist_list, bins=10, color='red')
            foo = 'bar'

            # I also want to look at what the distribution (in terms of abundance) of these back ground sequences
            # is. I
            print()
            cont = []
            for i, sample_name in enumerate(samples_to_plot):
                sys.stdout.write(f'\r{sample_name}')
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                # sort the list and
                cont.append([consolidated_host_seqs_abund_dict[_] for _ in sorted(consolidated_host_seqs_abund_dict, key=consolidated_host_seqs_abund_dict.get, reverse=True)])

            average_list = []
            list_index = 0
            print()
            while True:
                sys.stdout.write(f'\r{list_index}')
                temp_list = []
                for abund_list in cont:
                    try:
                        temp_list.append(abund_list[list_index])
                    except IndexError:
                        pass
                if not temp_list:
                    break
                else:
                    average_list.append(sum(temp_list)/len(temp_list))
                    list_index += 1
            foo = 'bar'

            # TODO plot up the stacked bars according to
            fig, ax = plt.subplots(1, 1)
            ax.scatter(x=range(len(average_list)), y=average_list, c='black', s=2)
            ax.set_ylim(-0.0003, 0.003)
            foo = 'bar'



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

# TODO other things to plot that will be useful.
# We can have a look to see how many sequences samples are losing and correlate this
# with the total number of unique sequences they started with

EighteenSAbundancePlotting().plot_master_scatter()