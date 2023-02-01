import numpy as np
from scipy.special import factorial, factorial2
from scipy.special import gamma as sciGamma

from thecobs.constants import *

def generalizedFermiFunction(W, Z, R, k):
    """Implementation of the generalized Fermi function F_{k-1} according to Behrens et al.

    :param W: Total electron energy in units of its rest mass
    :param Z: Proton number of daughter
    :param R: nuclear radius in natural units
    :param k: absolute value of kappa

    """
    p = np.sqrt(W**2.0-1.)
    y = ALPHA*Z*W/p
    gammak = np.sqrt(k**2.0-(ALPHA*Z)**2.0)
    prefactor = (k*factorial2(2*k-1))**2.0*4**k*(2*p*R)**(2.*(gammak-k))*np.exp(np.pi*y)
    gammas = (np.absolute(sciGamma(gammak+y*1j))/sciGamma(1+2.*gammak))**2.0
    return prefactor*gammas

def lambda_2(W, Z):
    p = np.sqrt(W**2.0-1.)
    y = ALPHA*Z*W/p
    gamma2 = np.sqrt(4-(ALPHA*Z)**2.0)
    gamma1 = np.sqrt(1-(ALPHA*Z)**2.0)
    R = 1.2e-15*(2.5*Z)**(1./3.)/NATURALLENGTH
    prefactor = (2+gamma2)/2/(1+gamma1)/(1+(ALPHA*Z)**2.0/12*(11.577216))
    kinematic = (2*p*R)**(2*(gamma2-gamma1-1))*(np.absolute(sciGamma(gamma2+y*1j))/np.absolute(sciGamma(gamma1+y*1j)))**2.0
    return prefactor*kinematic

def lambda_k(W, Z, k):
    """Coulomb function $\lambda_k$ as per Behrens et al.

    :param W: Total electron energy in units of its rest mass
    :param Z: Proton number of daughter
    :param k: absolute value of kappa

    """
    #return 1.
    gammak = np.sqrt(k**2.0-(ALPHA*Z)**2.0)
    gamma1 = np.sqrt(1.-(ALPHA*Z)**2.0)
    R = 1.2e-15*(2.5*Z)**(1./3.)/NATURALLENGTH
    return generalizedFermiFunction(W, Z, R, k)/generalizedFermiFunction(W, Z, R, 1)*(k+gammak)/(k*(1+gamma1))

def calcQ_k(k, Z, W, beta):
    """Calculate the Q screening ratio according to Buehring (1984)

    :param k: absolute value of kappa
    :param Z: Proton number of daughter
    :param W: Total electron energy in units of its rest mass
    :param beta: screening exponent at r=0

    """
    gamma_k = sqrt(k * k - (ALPHA * Z) ** 2.0)

    up = (ALPHA * Z) ** 2.0 * calcN2(-k, Z, W, beta) + (k + gamma_k) ** 2.0 * calcN2(k, Z, W, beta)

    down = (ALPHA * Z) ** 2.0 * calcNCoul2(-k, Z, W) + (k + gamma_k) ** 2.0 * calcNCoul2(k, Z, W)

    return up / down


def calcN2(kappa, Z, W, beta):
    """Calculate the square of the normalization function according to Buehring (1984)

    :param k: absolute value of kappa
    :param Z: Proton number of daughter
    :param W: Total electron energy in units of its rest mass
    :param beta: screening exponent at r=0
    :param kappa:

    """
    k = abs(kappa)

    Wt = W - ALPHA * Z * beta / 2.0

    p = sqrt(W ** 2.0 - 1.0)

    pt = 0.5 * p + 0.5 * sqrt(p ** 2.0 - 2 * ALPHA * Z * beta * Wt + 0.0 * (1j))

    kappat = kappa * sqrt(1.0 - ALPHA * Z * beta / (1.0 + W) + 0.0 * 1j)

    P = p / beta
    Pt = pt / beta

    gamma_k = sqrt(k * k - (ALPHA * Z) ** 2.0)

    yt = ALPHA * Z * Wt / pt

    term1 = 0.5 * ((kappa - gamma_k) * (kappa * Wt - gamma_k) + 0.5 * (gamma_k / kappa) * (kappat - kappa) ** 2.0)
    R_ab = 1.0
    if (kappa > 0):
        R_ab += 0.25 * k ** 2.0 * (beta / p) ** 2.0

    c1 = 1.0 + 2.0 * gamma_k
    c2 = gamma_k + yt * (1j)
    c3 = gamma_k + 2.0 * Pt * (1j)
    c4 = k + 2.0 * P * (1j)

    Gc1 = fabs(gamma(c1))
    Gc2 = fabs(gamma(c2))

    Gc3 = log(gamma(c3))
    Gc4 = log(gamma(c4))

    Gc34 = exp(Gc3 - Gc4)

    rap2 = (abs(Gc34)) ** 2.0

    if c3.imag > 200.0 and c4.imag > 200.0:
        crap2 = sqrt(c4/c3)*exp(c4-c3)
        lnrap = c3*log(c3)-c4*log(c4)
        rap2 = (abs(crap2*exp(lnrap)))**2.0

    rap1 = (Gc2 / Gc1) ** 2.0

    rap = rap1 * rap2

    term5 = beta ** (2.0 * gamma_k - 2.0 * k) * (2.0 * p) ** (2.0 * k - 1.0)

    return term1 * (rap / R_ab) * term5 + 0.0 * 1j


def calcNCoul2(kappa, Z, W):
    """Calculate the square of the normalization function for a simple point charge, according to Buehring (1984)

    :param kappa: absolute value of kappa
    :param Z: Proton number of daughter
    :param W: Total electron energy in units of its rest mass

    """
    ALPHA = 1.0 / 137.0

    k = abs(kappa)

    gamma_k = sqrt(k * k - (ALPHA * Z) ** 2.0)

    p = sqrt(W ** 2.0 - 1.0)

    y = ALPHA * Z * W / p

    return 0.5 * (kappa - gamma_k) * (kappa * W - gamma_k) * (gamma(1.0 + 2.0 * gamma_k)) ** (-2.0) * (
        np.absolute(gamma(gamma_k + y * 1j))) ** 2.0 * exp(pi * y) * (2.0 * p) ** (
                       2.0 * gamma_k - 1.0) + 0.0 * 1j

def calcLambda_k_scr(k, Z, W, beta):
    """Calculate the screening correction to the $\lambda_k$ Coulomb function

    :param k: absolute value of kappa
    :param Z: Proton number of daughter
    :param W: Total electron energy in units of its rest mass
    :param beta: Screening exponent at r=0

    """
    return calcQ_k(k, Z, W, beta)/calcQ_k(1.0, Z, W, beta)
