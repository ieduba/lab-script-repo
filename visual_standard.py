import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler

plt.rcdefaults()
# Could use Style Guide Instead of Custom
thickness = 3
fsize = 18
mpl.rcParams['lines.linewidth'] = thickness
mpl.rcParams['lines.linestyle'] = '-'
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.titlesize'] = fsize
mpl.rcParams['axes.labelsize'] = fsize

mpl.rcParams['xtick.labelsize'] = fsize#-8
mpl.rcParams['ytick.labelsize'] = fsize#-8

# Settings ensure that textboxes are recognizable in illustrator once exported
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


#Set Border Width
mpl.rcParams['axes.linewidth'] = 2

#Tick Mark Settings
mpl.rcParams['xtick.major.size'] = 3*thickness
mpl.rcParams['xtick.major.width'] = thickness
mpl.rcParams['ytick.major.size'] = 3*thickness
mpl.rcParams['ytick.major.width'] = thickness

# Fonts
mpl.rcParams['font.sans-serif'] = 'Helvetica'
# mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['axes.prop_cycle'] = cycler(color=['#0072B2', '#D55E00', '#009E73', '#CC79A7','darkgrey', '#56B4E9','#E69F00','#F0E442']) # can add black: '#000000'
