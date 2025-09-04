import matplotlib.pyplot as plt


plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 14})


def plot_levels(calculated, experimental):
    fig, ax = plt.subplots()
    unit = "MeV"
    if unit == "keV":
        unit_conversion = 1
    else:
        unit_conversion = 1 / 1000
    width = 2
    space = 1
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    for level in calculated.keys():
        try:
            ax.hlines(y=experimental[level] * unit_conversion, xmin=0, xmax=width, linewidth=2, color='k')
            plt.text(x=width, y=experimental[level] * unit_conversion, s="$" + level + "$")
            ax.hlines(y=calculated[level] * unit_conversion, xmin=space+width, xmax=space+width*2, linewidth=2, color='k')
            plt.text(x=space+width*2, y=calculated[level] * unit_conversion, s="$" + level + "$")
        except:
            pass
    plt.ylabel("Energy ({})".format(unit))
    plt.text(width/3, y=-500*unit_conversion, s="Expt.")
    plt.text(width + space + width/3, y=-500*unit_conversion, s="IBM")
    plt.xlim([-space/2, width*2+space*2])
    plt.ylim([-900 * unit_conversion, max(list(calculated.values()) + list(experimental.values()))*unit_conversion + 100*unit_conversion])

