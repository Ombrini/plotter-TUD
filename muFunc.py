# chemical potentials from concentrations 
# LFP
import numpy as np

import mpet.geometry as geo
from mpet.config import Config, constants


T = 1
k = constants.k                      # Boltzmann constant, J/(K Li)
Tref = constants.T_ref               # Temp, K
e = constants.e                      # Charge of proton, C
F = constants.F                      # C/mol
c_ref = constants.c_ref
eokT = constants.e / (constants.k * constants.T_ref)

def non_homog_rect_fixed_csurf(y, ybar, B, kappa, ywet):

    N = len(y)
    ytmp = np.empty(N+2, dtype=object)
    ytmp[1:-1] = y
    ytmp[0] = ywet
    ytmp[-1] = ywet
    dxs = 1./N
    curv = np.diff(ytmp, 2)/(dxs**2)
    muR_nh = -kappa*curv + B*(y - ybar)
    return muR_nh

def LiFePO4(y, ybar, resultDir):
    config = Config.from_dicts(resultDir)
    Omga = config["c","Omega_a"][0][0]
    kappa = config["c","kappa"][0][0]
    B = config["c","B"]
    cwet = 0.98

    muR_IS = T*np.log(y/(1-y))
    enthalpyTerm = Omga*(1-2*y)
    muRhomog = muR_IS + enthalpyTerm

    muRtheta = -eokT*3.422
    muRnonHomog = non_homog_rect_fixed_csurf(y, ybar, B, kappa, cwet)
    muR = muRhomog + muRnonHomog
    # actR = np.exp(muR/T)
    # muR += muRtheta + muR_ref
    return muR