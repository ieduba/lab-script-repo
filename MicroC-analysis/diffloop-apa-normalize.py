import numpy as np
import matplotlib.pyplot as plt

#looptypes = ['scr-combo-differential_loops1','scr-combo-differential_loops2','dH1-combo-differential_loops1','dH1-combo-differential_loops2']
looptypes = ['scrscr','scrlow','lowlow','lowscr', 'scrcommon', 'lowcommon']

maxval = 0

apa = {}
for looptype in looptypes:
	apa[looptype] = np.loadtxt(f'{looptype}-APA.txt', dtype = str, delimiter = ',')
	for row in range(len(apa[looptype])):
		apa[looptype][row][0] = apa[looptype][row][0].lstrip('[')
		apa[looptype][row][-1] = apa[looptype][row][-1].rstrip(']')
	apa[looptype] = apa[looptype].astype(float)
	if np.max(apa[looptype]) > maxval:
		maxval = np.max(apa[looptype])

for looptype in looptypes:
	fig, ax = plt.subplots()
	heatmap = ax.imshow(apa[looptype], cmap = 'Reds', vmin = 0, vmax = maxval*0.5, extent=[-10,10,-10,10])
	bar = plt.colorbar(heatmap)
	plt.savefig(f'{looptype}-scaled_APA.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()


