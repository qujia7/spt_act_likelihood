import numpy as np
import warnings
from scipy.interpolate import interp1d
import os
file_dir = os.path.abspath(os.path.dirname(__file__))
from cobaya.likelihood import Likelihood
import healpy as hp
from scipy.stats import norm




class Omega_mh3(Likelihood):
    # Data type for aggregated chi2 (case sensitive)



    def initialize(self):
        self.norm = norm(loc=0.09635, scale=0.000001)


    def get_requirements(self):
        return {'h': None,'omegam':None}

    def logp(self, **params_values):
        h = self.provider.get_param("h")
        omegam = self.provider.get_param("omegam")
        return self.norm.logpdf(omegam*h**3)
 