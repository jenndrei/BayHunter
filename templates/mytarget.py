# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

from BayHunter.Targets import SingleTarget
import MyForwardModel


class MyOwnTarget(SingleTarget):
    # swd: exponential noise correlation law
    # rf: gaussian if rfnoise_corr is fixed, else exponential
    noiseref = 'swd'

    def __init__(self, x, y):
        ref = 'mydata'
        SingleTarget.__init__(self, x, y, ref)

        # forward your own plugin (instance),
        # otherwise it is None and returns an error
        self.moddata.plugin = MyForwardModel(x, ref)
        self.moddata.xlabel = 'xvalues in unit'
