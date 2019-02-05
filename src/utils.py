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
    with open(configfile, 'rb') as f:
        data = pickle.load(f)

    return data


def get_path(name):
    fn = op.join(op.dirname(__file__), 'defaults', name)
    if not op.exists(fn):
        raise OSError('%s does not exist!' % name)
    return fn
