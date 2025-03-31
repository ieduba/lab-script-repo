import numpy as np
import matplotlib.pyplot as plt

looptypes = ['scr-scrdiffloop','scr-lowdiffloop','low-scrdiffloop','low-lowdiffloop','scrscr-rep','scrlow-rep','lowlow-rep','lowscr-rep','scrcommon-rep','lowcommon-rep','scr-CTCF','scr-noCTCF','low-CTCF','low-noCTCF','scr-upgenes','scr-nochangegenes','low-upgenes','low-nochangegenes']

print('looptype\tcenter/LL')

apa = {}
enrich = {}
for looptype in looptypes:
	apa[looptype] = np.loadtxt(f'{looptype}-APA.txt', dtype = str, delimiter = ',')
	for row in range(len(apa[looptype])):
		apa[looptype][row][0] = apa[looptype][row][0].lstrip('[')
		apa[looptype][row][-1] = apa[looptype][row][-1].rstrip(']')
	apa[looptype] = apa[looptype].astype(float)
	enrich[looptype] = apa[looptype][10][10]/apa[looptype][20][0]

	print(f'{looptype}\t{enrich[looptype]}')

