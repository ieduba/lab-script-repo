import numpy as np
from matplotlib import pyplot as plt
from sklearn.mixture import GaussianMixture

def model(counts, ID):
	frag_lens = range(len(counts))
	all_frags = []
	for i in range(0, len(frag_lens)):
		temp = np.repeat(frag_lens[i], int(counts[i]))*10
		all_frags = np.concatenate((all_frags, temp))

	frags_subset_idx = np.where(all_frags<1200)
	frags_subset = all_frags[frags_subset_idx]

	# Set up the dataset.
	X = frags_subset.reshape(-1,1)
	X_train = X

	min_range = 0
	max_range = 1199

	# fit models with 1-10 components
	N = np.arange(1, 11) #undo, make 1-11
	models = [None for i in range(len(N))]

	for i in range(len(N)):
		models[i] = GaussianMixture(N[i]).fit(X_train)

	M_best = models[5]# Manually select 6 components models
	weights = M_best.weights_
	peaks = M_best.means_

	fig, ax = plt.subplots()
	x = np.linspace(min_range, max_range, 1200)
	logprob = M_best.score_samples(x.reshape(-1, 1))
	responsibilities = M_best.predict_proba(x.reshape(-1, 1))
	pdf = np.exp(logprob)
	pdf_individual = responsibilities * pdf[:, np.newaxis]

	ax.hist(X_train, 120, density=True, histtype='stepfilled', alpha=0.4)
	ax.plot(x, pdf, '-k')
	ax.plot(x, pdf_individual, '--k')
	ax.text(0.04, 0.96, "Best-fit Mixture", ha='left', va='top', transform=ax.transAxes)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$p(x)$')

	plt.savefig(f'{ID}-GMMfit.pdf')
	
	return(weights, peaks)
