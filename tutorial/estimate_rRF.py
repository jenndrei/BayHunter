import numpy as np
import matplotlib.pyplot as plt

from BayHunter import utils

#
# ----------------------------------------------------------- RF and parameters
#
# load receiver function and define parameters
rfx, rfy = np.loadtxt('observed/st3_prf.dat').T  # synth RF, Gauss factor a=1
rfa = 1  # a
dt = 0.2  # RF sampling rate
draws = 50000
# more draws (= larger sample) show more precise and smoother results,
# but are slower in computation.

rrfs = [0.75, 0.85, 0.95, 0.97, 0.98, 0.99]  # test r_RFs

pars = {'rfx': rfx, 'rfy': rfy, 'rfa': rfa,
        'a': rfa, 'dt': dt, 'rrfs': rrfs,
        'draws': draws}


#
# ------------------------------------------ visualize 'raw' data and estimates
#

fig = utils.plot_rrf_estimate(pars=pars)
fig.savefig('st3_rrf_estimate.pdf', bbox_inches='tight')


#
# -------------------------------------- return values for costum visualization
#

# update test r_RFs for more continuous estimates
pars['rrfs'] = np.linspace(0.9, 0.999, 25)
pars['draws'] = 2000
# quicker computation by less draws (= smaller sample size from distribution),
# but also larger uncertainties of results. To compensate for the smaller
# sample, draw more samples by looping over a number.
# However, I recommend to increase the number of draws (e.g., 50 000)
# and wait a bit longer...
# As always, test different settings to find your best r_RF... estimate

plt.close()
fig, ax = plt.subplots()

for sample in range(10):
    rrf, a = utils.rrf_estimate(pars=pars)  # return values
    ax.plot(rrf, a, color='k', marker='x', ls='')

ax.axhline(rfa, color='gray', label='reference')
ax.set_xlabel('$r_{RF}$')
ax.set_ylabel('Gauss factor a')
ax.grid(color='lightgray')
ax.legend(loc=1)
fig.savefig('rrf-a_rel.pdf', bbox_inches='tight')
