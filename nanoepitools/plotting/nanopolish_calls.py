import matplotlib.pyplot as plt
import numpy as np


def plot_met_hist_types(metcall, thres=2.5, typecol='type', bins=200,
                        alpha=0.5, typedict=None, typecolor=None,
                        title='Methylation log likelihood ratios'):
    types = set(metcall[typecol])
    if typedict is None:
        typedict = {t:t for t in types}
    if typecol is None:
        typecolor = {t:None for t in types}
    for t in types:
        llr_type = metcall['log_lik_ratio'].loc[metcall[typecol] == t]
        llr_type = np.clip(llr_type, -20, 20)
        plt.hist(met_llr, bins=bins, alpha=alpha, color=typecolor[t],
                 label=typedict[t])
    plt.xlabel("Methylation log-likelihood ratio")
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()