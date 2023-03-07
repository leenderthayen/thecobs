import pandas as pd
import numpy as np

from importlib import resources

from thecobs.SpectralFunctions import *
from thecobs.Constants import *

def getExchangeParams(Z):
    with resources.open_text("thecobs.data", "ExchangeData.dat") as fid:
        data = np.genfromtxt(fid)
        index = np.where(data[:, 0] == Z)

        return data[index, 1:][0, 0]

def getExchangeParamsSimkovic(Z):
    with resources.open_text("thecobs.data", "ExchangeDataSimkovic.dat") as fid:
        data = np.genfromtxt(fid)
        index = np.where(data[:, 0] == Z)

        return data[index, 1:][0, 0]

def calculateSpectrum(E, params):
    W = 1 + E/ELECTRON_MASS_KEV

    ph = phase_space(W, **params)
    f = fermi_function(W, **params)
    l0 = finite_size_L0(W, **params)
    u = finite_size_U_fermi(W, **params)
    rc = radiative_correction(W, **params)

    rec = np.ones(len(E))
    recCoul = np.ones(len(E))
    C = np.ones(len(E))

    beta_type = params['Type']

    if beta_type == 'A: Fermi':
        rec = recoil_fermi(W, **params)
        recCoul = recoil_Coulomb_fermi(W, **params)
        C = shape_factor_fermi(W, **params)
    elif beta_type == 'A: Gamow-Teller':
        rec = recoil_gamow_teller(W, **params)
        recCoul = recoil_Coulomb_gamow_teller(W, **params)
        C = shape_factor_gamow_teller(W, **params)
    elif beta_type == 'A: Mixed':
        norm = 1+mixing_ratio**2
        recF = recoil_fermi(W, **params)
        recCoulF = recoil_Coulomb_fermi(W, **params)
        CF = shape_factor_fermi(W, **params)
        recGT = recoil_gamow_teller(W, **params)
        recCoulGT = recoil_Coulomb_gamow_teller(W, **params)
        CGT = shape_factor_gamow_teller(W, **params)
        rec = 1+(recF-1+mixing_ratio**2*(recGT-1))/norm
        recCoul = 1+(recCoulF-1+mixing_ratio**2*(recCoulGT-1))/norm
        C = 1+(CF-1 + mixing_ratio**2*(CGT-1))/norm
    elif 'FU' in beta_type:
        C = shape_factor_unique_forbidden(W, **params)
    s = atomic_screening(W, **params)
    X = atomic_exchange(W, **params)

    sp = ph*f*l0*u*rc*rec*recCoul*C*s*X

    comb = np.stack((E, W, sp, ph, f, l0, rc, C, s, X, u, rec, recCoul), axis=1)

    df = pd.DataFrame(comb, columns = ['Energy', 'W', 'Spectrum', 'PhaseSpace', 'FermiFunction', 'L0', 'RadiativeCorrections', 'ShapeFactor', 'Screening', 'Exchange', 'U', 'Recoil', 'CoulombRecoil'])

    return df
