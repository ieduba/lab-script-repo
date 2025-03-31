import numpy as np
from scipy.signal import savgol_filter
from matplotlib import pyplot as plt

for cond in ['HBR1-5ends', 'GBR1-5ends', 'HBR1-wholeread', 'GBR1-wholeread']:
	table = np.loadtxt(f'{cond}-metagene.tab', dtype=str, delimiter='\t', skiprows=2)
	smoothed = {}
	for row in range(len(table)):
		name = table[row,0]
		data = table[row,2:]
		smoothed[name] = savgol_filter(data, 3, 2)
		plt.plot(range(len(data)),data,label=name)

	plt.legend()
	plt.savefig(f'{cond}-smoothed-metagene.pdf', format='pdf', bbox_inches='tight')
	plt.close
