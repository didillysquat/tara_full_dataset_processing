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

class MSPlots(EighteenSBase):
    def __init__(self):
        super().__init__()

    def plot_three_row(self):
        """
        The three row plot.
        """
        tr = ThreeRow(parent=self).plot

class ThreeRow:
    def __init__(self, parent):
        self.parent = parent
        # Let's start with the first plot quick and dirty and then we can add the others
        # and refactorize.
        self.fig, self.ax = plt.subplots(1,1, figsize=(5,5))
        plt.savefig(os.path.join(self.parent.eighteens_dir, 'temp_fig.png'))
        self.foo = 'bar'

    def plot(self):
        foo = 'bar'

MSPlots().plot_three_row()