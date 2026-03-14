import matplotlib.pyplot as plt

def plotdefaults():
    plot_params = {'axes.linewidth': 1,
                   'xtick.labelsize': 'medium',
                   'ytick.labelsize': 'medium',
                   'xtick.major.pad': 8,
                   'xtick.major.size': 12,
                   'xtick.minor.size': 6,
                   'ytick.major.size': 12,
                   'ytick.minor.size': 6,
                   'xtick.direction': 'in',
                   'ytick.direction': 'in',
                   'patch.linewidth': 0,
                   'font.size': 12}
    plt.rcParams.update(plot_params)
