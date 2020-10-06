# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import zmq
import pickle
import numpy as np
import os.path as op
from configobj import ConfigObj
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

rstate = np.random.RandomState(333)


class SerializingSocket(zmq.Socket):
    """A class with some extra serialization methods
    send_array sends numpy arrays with metadata necessary
    for reconstructing the array on the other side (dtype, shape).
    """
    def send_array(self, arr, flags=0, copy=True, track=False):
        md = dict(
            dtype=str(arr.dtype),
            shape=arr.shape,
        )
        self.send_json(md, flags | zmq.SNDMORE)
        return self.send(arr, flags, copy=copy, track=track)

    def recv_array(self, flags=0, copy=True, track=False):
        md = self.recv_json(flags=flags)
        msg = self.recv(flags=flags, copy=copy, track=track)
        arr = np.frombuffer(msg, dtype=md['dtype'])
        return arr.reshape(md['shape'])


class SerializingContext(zmq.Context):
    _socket_class = SerializingSocket


def string_decode(section):
    keywords = ['station', 'savepath']

    for key in section:
        if key in keywords:
            continue
        try:
            section[key] = eval(section[key])
        except:
            for i, value in enumerate(section[key]):
                section[key][i] = eval(value)
    return section


def load_params(initfile):
    config = ConfigObj(initfile)
    keywords = ['station', 'savepath']
    params = []
    for configsection in config.sections:
        if configsection == 'datapaths':
            continue
        section = config[configsection]
        section = string_decode(section)
        params.append(section)
    return params


def load_params_user(initfile, station, slowness=7):
    import linecache
    config = ConfigObj(initfile)

    paths = {}
    for key in config['datapaths']:
        if key.split('.')[-1] == 'bin':
            file = config['datapaths'][key] % (station, slowness)
        else:
            file = config['datapaths'][key] % station

        if op.exists(file):
            newkey = key.split('_')[-1]
            paths[newkey] = file

            # only for receiver functions
            if key.split('.')[-1] == 'bin':
                slow = float(linecache.getline(file, 2).strip().replace('#', ''))
                paths['slowness.bin'] = slow

            if key.split('.')[-1] == 'stack':
                slow = float(linecache.getline(file, 2).strip().replace('#', ''))
                paths['slowness.stack'] = slow

    modelpriors = string_decode(config['modelpriors'])
    initparams = string_decode(config['initparams'])
    initparams['station'] = station
    initparams['savepath'] = initparams['savepath'] % (station, '%.2f')
    return paths, modelpriors, initparams


def save_baywatch_config(targets, path='.', priors=dict(), initparams=dict(),
                         refmodel=dict()):
    """
    Saves a configfile that you will need if using BayWatch.

    targets: JointTarget instance fed into the inversion
    path: where to save the configfile
    priors: used for inversion
    refmodel: reference model / values to be illustrated
    """
    configfile = op.join(path, 'baywatch.pkl')
    data = {}

    for target in targets.targets:
        target.get_covariance = None

    data['targets'] = targets.targets
    data['priors'] = priors
    data['initparams'] = initparams
    data['refmodel'] = refmodel

    with open(configfile, 'wb') as f:
        pickle.dump(data, f)


def save_config(targets, configfile, priors=dict(), initparams=dict()):
    """
    Conveniently saves a configfile that you can easily use to view the data
    and parameters used for inversion. This configfile (.pkl) will also be used
    for PlotFromStorage plotting methods after the inversion. With this you can
    redo the plots with the correct data used for the inversion.

    targets: JointTarget instance from inversion
    configfile: outfile name
    priors, initparams: parameter dictionaries important for plotting,
    contains e.g. prior distributions, noise params, iterations etc.
    """
    data = {}
    refs = []

    for target in targets.targets:
        target.get_covariance = None
        ref = target.ref
        refs.append(ref)

    data['targets'] = targets.targets
    data['targetrefs'] = refs
    data['priors'] = priors
    data['initparams'] = initparams

    with open(configfile, 'wb') as f:
        pickle.dump(data, f)


def read_config(configfile):
    try:  # python2
        with open(configfile, 'rb') as f:
            data = pickle.load(f)
    except:  # python3
        with open(configfile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')

    return data


def get_path(name):
    fn = op.join(op.dirname(__file__), 'defaults', name)
    if not op.exists(fn):
        raise OSError('%s does not exist!' % name)
    return fn


# following functions for estimate r_RF
def _compute_gaussnoise(size, corr=0.85, sigma=0.0125, draws=1):
    """Gaussian correlated noise - use for RF if Gauss filter applied."""
    # from BayHunter SynthObs
    idx = np.fromfunction(lambda i, j: (abs((i+j) - 2*i)),
                          (size, size))
    rmatrix = corr**(idx**2)

    Ce = sigma**2 * rmatrix
    data_noise = rstate.multivariate_normal(np.zeros(size), Ce, draws)

    return np.concatenate(data_noise)


def compute_spectrum(y, Fs):
    """
    Return (normed) single-sided amplitude spectrum of y(t).
    y: receiver function amplitudes
    Fs: Frequency sampling (df) in Hz.
    """
    y = y - np.mean(y)

    n = y.size  # length of the signal
    n_half = int(n/2.)
    k = np.arange(n)
    T = n/Fs
    frq = k/T
    frq = frq[:n_half]  # frequency range

    Y = np.fft.fft(y)/n  # fft computing and normalization
    Y = Y[:n_half]
    Y = abs(Y)
    Y = Y/Y.max()

    return frq, Y


def gauss_fct(a, x):
    """Return Gaussian curve."""
    return np.exp(-(x*2*np.pi)**2 / (4*a**2))


def _min_fct(a, x, y):
    """Minimizing function for Gaussian bell and to be regressed x and y."""
    return gauss_fct(a, x) - y


def _spec_resample(frq, Y):
    """Computes a 2-D histogram of frequencies and spectral energy, and
    resamples the large amount of data to 120 bins.
    """
    bins = 120
    limit = 3  # minimum occurrences per bin
    y_values = np.zeros((bins)) * np.nan  # Y

    hist, xedges, yedges = np.histogram2d(frq, Y, bins=bins)
    xbin = (xedges[:-1] + xedges[1:])/2.
    ybin = (yedges[:-1] + yedges[1:])/2.

    ybin = ybin[::-1]

    histp = hist.T
    histp = histp[::-1]
    for i_y, row in enumerate(histp):
        # going through the rows, starting at top plot (no energy, all frequencies)
        for i_x, occurence in enumerate(row):
            if y_values[i_x] > 0:
                continue
            elif occurence > limit:
                y_values[i_x] = ybin[i_y]

    return xbin, y_values


def plot_rrf_estimate(pars=dict()):
    """ Returns a figure illustrating RF, RF-spectrum, reference Gaussian filter of given
    Gaussfactor [a] and randomly generated noise based on [rrfs] with correlated Gauss 
    factor as displayed in the legend. To estimate a proper rrf for the input RF data,
    choose a rrf-curve that best matches the Gaussian factor used for RF computation.
    As there is still a random factor included (which is reduced by the large amount of 
    drawings), repeat several times to ensure steadiness of your chosen r.
    Please forward all parameters as listed below, otherwise, there might be incompatible
    default values:

    rfx: RF, time in s, default np.linspace(-5, 35, 201)
    rfy: RF, amplitude, default None
    rfa: RF, Gauss factor for label only, default None
    rrfs: rrf to be inspected, default [0.75, 0.85, 0.95]
    a: Gauss factor for reference curve, default 2
    dt: time sampling, default computation from rfx (from rfx default: 0.2)
    draws: number of draws for random noise generation, default 50000

    """
    rfx = pars.get('rfx', np.linspace(-5, 35, 201))
    rfy = pars.get('rfy', None)
    rfa = pars.get('rfa', None)

    rfdt = np.median(np.unique(rfx[1:] - rfx[:-1]))

    rrfs = pars.get('rrfs', [0.75, 0.85, 0.95])
    a = pars.get('a', 2.)  # reference plotting
    dt = pars.get('dt', rfdt)
    df = 1./dt

    fig = plt.figure()

    if rfx is not None and rfy is not None:
        ax_rf = fig.add_subplot(2, 1, 1)
        try:
            label = 'RF, a=%.1f' % rfa
        except Exception:
            label = 'RF'

        ax_rf.plot(rfx, rfy, 'k', lw=1, label=label)
        ax_rf.set_xlabel('Time in s')
        ax_rf.set_ylabel('Amplitude')
        ax_rf.set_xlim(rfx.min(), rfx.max())

        ax_rf.legend(loc=1)

        frq, Y = compute_spectrum(rfy, df)

        ax_p = fig.add_subplot(2, 1, 2)
        ax_p.plot(frq, Y, 'k', lw=1, label='RF-spec', zorder=200)

        print('Time sampling given for RF spectrum: %.3f s' % dt)

    # ----------------------------------------

    else:
        ax_p = fig.add_subplot(1, 1, 1)

    draws = pars.get('draws', 50000)
    sigma = 0.0125  # should not matter
    rrfs = np.array(rrfs)
    rrfs.sort()
    a0 = 1

    Y_all = np.zeros((rrfs.size, int(draws*rfx.size/2.)))

    print('a\trrf')
    for rrf in rrfs:
        rfnoise = _compute_gaussnoise(rfx.size, rrf, sigma, draws=draws)
        frq, Y = compute_spectrum(rfnoise, df)

        # find envelope of spectrum, first resample values
        res_frq, res_Y = _spec_resample(frq, Y)
        res_Y_max = res_Y.max()
        res_Y = res_Y / res_Y.max()

        env_lsq = least_squares(_min_fct, a0, args=(res_frq, res_Y))
        env_a = env_lsq.x[0]
        env_G = gauss_fct(env_a, res_frq)
        label = 'a=%.1f; $r_{RF}$=%.2f' % (env_a, rrf)
        line, = ax_p.plot(res_frq, env_G, lw=1.2, zorder=100, label=label)
        color = line.get_color()

        # ax_p.plot(res_frq, res_Y/res_Y_max, marker='x', zorder=200, lw=1.2,
        #         color=color, label=label)

        ax_p.plot(frq, Y/res_Y_max, lw=0.3, alpha=0.5, color=color)

        print('%.3f\t%.3f' % (env_a, rrf))

    ax_p.set_xlabel('Frequency in Hz')
    ax_p.set_ylabel('Spectral Power')
    ax_p.set_ylim(ymin=0)
    ax_p.set_xlim(frq.min(), frq.max())

    # reference curve for Gaussian 'a' 
    G = gauss_fct(a, res_frq)
    ax_p.plot(res_frq, G, label='a=%.1f' % a, color='k', ls='--', zorder= 200)

    # legend sort by a
    handles, labels = ax_p.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    # ax.legend(handles, labels)
    ax_p.legend(handles[::-1], labels[::-1], loc=2, bbox_to_anchor=(1,1.1))

    fig.subplots_adjust(hspace=0.4)
    return fig


def rrf_estimate(pars=dict()):
    """ Returns rrf and a pairs. See explanation for plot_rrf_estimate

    rfx: RF, time in s, default np.linspace(-5, 35, 201)
    rrfs: rrf to be inspected, default [0.75, 0.85, 0.95]
    dt: time sampling, default computation from rfx (from rfx default: 0.2)
    draws: number of draws for random noise generation, default 50000
    """
    rfx = pars.get('rfx', np.linspace(-5, 35, 201))
    rfdt = np.median(np.unique(rfx[1:] - rfx[:-1]))

    rrfs = pars.get('rrfs', [0.75, 0.85, 0.95])
    dt = pars.get('dt', rfdt)
    df = 1./dt

    draws = pars.get('draws', 50000)
    sigma = 0.0125  # should not matter
    rrfs = np.array(rrfs)
    rrfs.sort()
    a0 = 1

    Y_all = np.zeros((rrfs.size, int(draws*rfx.size/2.)))
    a_est = []
    for rrf in rrfs:
        rfnoise = _compute_gaussnoise(rfx.size, rrf, sigma, draws=draws)
        frq, Y = compute_spectrum(rfnoise, df)

        # find envelope of spectrum, first resample values
        res_frq, res_Y = _spec_resample(frq, Y)
        res_Y_max = res_Y.max()
        res_Y = res_Y / res_Y.max()

        env_lsq = least_squares(_min_fct, a0, args=(res_frq, res_Y))
        env_a = env_lsq.x[0]

        a_est.append(env_a)
        print('%.3f\t%.3f' % (env_a, rrf))

    return rrfs, a_est