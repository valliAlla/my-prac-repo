!pip install bioinfokit
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
matplotlib.style.use('ggplot')
%matplotlib inline
 x = pd.read_csv("/Users/allanagalalithavalli/Documents/Internship/BOK KO2.csv")
 ids= x[((x["log(FC)"]<=-0.5) | (x["log(FC)"]>=0.5)) & ((x["log(P)"]<=-1) | (x["log(P)"]>=1)) & (x["p value"]>=0.05)]["Symbol"].tolist()
 from bioinfokit import visuz
visuz.volcano(table="/Users/allanagalalithavalli/Documents/Internship/BOK KO2.csv", lfc="log(FC)",lfc_thr=0.5, pv="p value", geneid="Symbol", genenames= tuple(ids),valpha=0.4)
visuz.gene_exp.volcano(x, lfc="log(FC)", pv="p value",
                       lfc_thr=0.5, valpha=0.4, geneid="Symbol", genenames=tuple(ids))
class Volcano(object):
def __init__(self, ratio, p_val, label=None, s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5):
assert len(ratio) == len(p_val)
        self.df = pd.DataFrame({"ratio": ratio, "p_val": p_val})
        if label is not None:
            self.df["label"] = label
        self.s_curve_y_axis_overplot = s_curve_y_axis_overplot
        self.p_val_cutoff = self.get_p_val_cutoff()
        self.ratio_cutoff = self.get_ratio_cutoff()
        self.df["s_val"] = self.df.apply(self.calc_s_from_row, axis=1)
        self.ratio_for_s = pd.Series(np.linspace(self.df["ratio"].min() - s_curve_x_axis_overplot,
                                                 self.df["ratio"].max() + s_curve_x_axis_overplot,
                                                 num=1000))
        self.p_for_s_larger_1 = self.ratio_for_s.apply(self.calc_p_for_s_equals_1)

def get_p_val_cutoff(self):
"""
        p_val_cutoff = 0.05
        pc = 3.5 + median(p_val(50% lowest log2_ratios)) --> is what Jan uses for whatever reason ???
        -log10_pval of 2.0 --> pval of 0.01
        """
quant = self.df["ratio"].quantile(0.5)
        return 1.3
 def get_ratio_cutoff(self):
        """
        log2_ratio_cutoff = 2.0 
        ratio_cutoff_high = 2 + median(ratio(50% lowest log10_p_values))
        ratio_cutoff_low = 0.5 - median(ratio(50% lowest log10_p_values))        
        """
  quant = self.df["p_val"].quantile(0.5)
        median_ = self.df.loc[self.df["p_val"] < quant, "ratio"].median()
        
        ratio_cutoff_high = median_  #1 + median_
        ratio_cutoff_low = -1* median_  #-1 + median_
        return ratio_cutoff_low, ratio_cutoff_high

    def calc_s_from_row(self, row):
        p_val = row["p_val"]
        ratio = row["ratio"]
        return self.calc_s(p_val, ratio)

    def calc_s(self, p_val, ratio):
        """
        so the algorithmn for finding stuff with s > 1 is:
        discard stuff below the ratio_cutoff
        discard stuff below the p-val cutoff
        do the calcuation for the stuff above BOTH cutoffs and accept all with s > 1
        s = (p_val - p_val_cutoff) * (ratio - ratio_cutoff)
        :param p_val: Float(-log10 p-value)
        :param ratio: Float(log2 ratio)
        :return: Float
        """
        ratio_cutoff_low, ratio_cutoff_high = self.ratio_cutoff
        if ratio > 0:
            ratio_delta = ratio - ratio_cutoff_high
            if ratio_delta < 0:
                return 0
        elif ratio < 0:
            ratio_delta = ratio - ratio_cutoff_low
            if ratio_delta > 0:
                return 0
        ratio_delta = abs(ratio_delta)
        p_val_delta = p_val - self.p_val_cutoff
        if p_val_delta < 0:
            return 0
        return p_val_delta * ratio_delta

    def calc_p_for_s_equals_1(self, ratio):
        """
        :param ratio: Float(log2 ratio)
        :return: Float
        """
        ratio_cutoff_low, ratio_cutoff_high = self.ratio_cutoff
        ratio_delta_high = ratio - ratio_cutoff_high
        ratio_delta_low = ratio - ratio_cutoff_low

        if ratio > ratio_cutoff_high:
            return (1.0 / ratio) + self.p_val_cutoff
        elif ratio < ratio_cutoff_low:
            return (1.0 / (ratio * -1)) + self.p_val_cutoff
        else:
            return np.nan

    def get_fig(self, title="Volcano plot", s_value_cutoff=1.0):
        fig, ax1 = plt.subplots(figsize=(12, 12))
        ax1.set_title(title, fontsize=22, fontweight="bold")

        cond_p_val = self.p_for_s_larger_1 <= (self.df["p_val"].max() + self.s_curve_y_axis_overplot)
        cond_pos = self.ratio_for_s > 0
        x1 = self.ratio_for_s[cond_pos & cond_p_val]
        y1 = self.p_for_s_larger_1[cond_pos & cond_p_val]

        cond_neg = self.ratio_for_s < 0
        x2 = self.ratio_for_s[cond_neg & cond_p_val]
        y2 = self.p_for_s_larger_1[cond_neg & cond_p_val]

        x3 = self.df["ratio"]
        y3 = self.df["p_val"]

        ax1.plot(x1, y1, 'r-', x2, y2, 'r-', x3, y3, 'b.', markersize=10)  # , alpha=0.7)
        ax1.set_xlabel('log2(ratio)', fontsize=20)
        ax1.set_ylabel('-log10(p-value)', fontsize=20)

        if "label" in self.df.columns:
            cond = self.df["s_val"] >= s_value_cutoff
            for index_, row in self.df[cond].iterrows():
                label = row["label"]
                x_coord = row["ratio"]
                y_coord = row["p_val"]
                if x_coord < 0:
                    ax1.plot(x_coord,y_coord,'r.',markersize=10)
                    ax1.annotate(label,xy=(x_coord, y_coord),xycoords='data',xytext=(5,5),
                        textcoords='offset points',arrowprops=dict(arrowstyle="-"),fontsize=12,color='r')
                else:
                    ax1.plot(x_coord,y_coord,'g.',markersize=10)
                    ax1.annotate(label, xy=(x_coord, y_coord), xycoords='data', xytext=(5, 5),
                    textcoords='offset points', arrowprops=dict(arrowstyle="-"), fontsize=12,color="g")



                #ax1.annotate(label, xy=(x_coord, y_coord), xycoords='data', xytext=(5, 5),
                    #textcoords='offset points', arrowprops=dict(arrowstyle="-"), fontsize=12)

        return fig, ax1
       volc = Volcano(x["log(FC)"], x["log(P)"], x["Symbol"],
               s_curve_x_axis_overplot=2, s_curve_y_axis_overplot=.5)
fig, axs = volc.get_fig()
axs.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
