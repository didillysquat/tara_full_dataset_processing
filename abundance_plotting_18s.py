"""We need to do some investigation into the absolute number of sequences and unique number of sequences
found in samples to better quality control the distance metrics"""

from base_18s import EighteenSBase
import compress_pickle
import os
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import sys

class EighteenSAbundancePlotting(EighteenSBase):
    def __init__(self):
        super().__init__()
        self.absolute_consolidated_abundance_dict = compress_pickle.load(os.path.join(self.cache_dir, 'consolidated_df_dict_output_tables.p.bz'))
        self.fig, self.ax = plt.subplots(3,1)

    def plot_master_scatter(self):
        sample_names = list(self.absolute_consolidated_abundance_dict.keys())
        host_only_master_seq_info_dict = compress_pickle.load(os.path.join(self.cache_dir, 'host_only_master_seq_info_dict.p.bz'))
        ordered_seq_names = [tup[0] for tup in sorted([_ for _ in host_only_master_seq_info_dict.items()],
                                  key=lambda x: x[1],
                                  reverse=True)]
        # Total number of seqs
        x = [sum(self.absolute_consolidated_abundance_dict[k]) for k in sample_names]
        # Unique number of seqs
        y = [len([_ for _ in v if _ > 0]) for v in [self.absolute_consolidated_abundance_dict[k] for k in sample_names]]
        self.ax[0].scatter(x=x, y=y, c='black', s=2)
        self.ax[0].set_xlim((self.ax[0].get_xlim()[0], 500000))


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
                seq = ordered_seq_names[i]
                abund = abund_list[i]
                # prob = abund/tot
                redundant_list.extend([seq for _ in range(abund)])
                # prob_list.extend([prob for _ in range(abund)])

            hard_sub_sample_list = np.random.choice(redundant_list, 10000, replace=False)
            hard_abunds = list(Counter(hard_sub_sample_list).values())
            hard_sub_sample_dict[sample_name] = hard_abunds

            # For soft
            norm_abunds = [int((abund_list[i]/tot)*10000) for i in non_z]
            rel_sub_sample_dict[sample_name] = norm_abunds

        sample_names = hard_sub_sample_dict.keys()
        # Hard x total
        x_hard = [sum(hard_sub_sample_dict[k]) for k in sample_names]
        # Hard y unique
        y_hard = [len([_ for _ in v if _ > 0]) for v in [hard_sub_sample_dict[k] for k in sample_names]]

        # Rel x total
        x_rel = [sum(rel_sub_sample_dict[k]) for k in sample_names]
        # Rel y unique
        y_rel = [len([_ for _ in v if _ > 0]) for v in [rel_sub_sample_dict[k] for k in sample_names]]

        # Plot hard
        self.ax[1].scatter(x=list(range(len(y_hard))), y=y_hard, c='black', s=2)
        # Plor rel
        self.ax[2].scatter(x=list(range(len(y_rel))), y=y_rel, c='black', s=2)
        foo = 'bar'


EighteenSAbundancePlotting().plot_master_scatter()