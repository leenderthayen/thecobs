import numpy as np
from scipy.special import gamma, loggamma, spence

ALPHA = 1/137.036
PROTON_MASS_C2 = 938.272
NEUTRON_MASS_C2 = 939.565
ELECTRON_MASS_C2 = 0.511
NUCLEON_MASS_C2 = (PROTON_MASS_C2+NEUTRON_MASS_C2)/2
GA = 1.27
GM = 4.706

def phase_space(W, W0):
    return np.sqrt(W**2-1)*(W-W0)**2*W

def fermi_function(Z, W, R):
    if Z == 0:
        return 1

    f = np.zeros(W.shape)
    mask = (W>1)
    W = W[mask]

    g = np.sqrt(1-(ALPHA*Z)**2)
    p = np.sqrt(W**2-1)
    y = ALPHA*Z*W/p

    #We use the traditional Fermi function, i.e. a prefactor 4 instead of 2(1+gamma)
    f[mask] = (4
            *np.power(2*p*R, 2*(g-1))
            *np.exp(np.pi*y)
            *(np.abs(gamma(g+1.j*y))/(gamma(1+2*g)))**2)

    return f

def finite_size_L0(Z, W, R):
    return (1
            -ALPHA*Z*W*R
            +13/60*(ALPHA*Z)**2
            -ALPHA*Z*R/W)

def finite_size_U_fermi(Z, W):
    p = np.sqrt(W**2-1)

    a0 = (-5.6e-5
            -4.94e-5*Z
            +6.23e-8*Z**2)
    a1 = (5.17e-6
            +2.517e-6*Z
            +2e-8*Z**2)
    a2 = (-9.17e-8
            +5.53e-9*Z
            +1.25e-10*Z**2)

    return 1+a0+a1*p+a2*p**2

def sirlin_g(W, W0):
    g = np.zeros(W.shape)
    mask = (W > 1) & (W<W0)
    W = W[mask]

    beta = np.sqrt(W**2-1)/W

    g[mask] = (3*np.log(PROTON_MASS_C2/ELECTRON_MASS_C2)
        -3./4.
        +4./beta*(-1*spence(1-(2*beta/(1+beta))))
        +4*(np.arctanh(beta)/beta-1)*((W0-W)/(3*W)-3/2+np.log(2*(W0-W)))
        +np.arctanh(beta)/beta*(2*(1+beta**2)+(W0-W)**2/(6*W**2)))

    return g

def radiative_correction_o1(W, W0):
    g = sirlin_g(W, W0)

    return ALPHA/2/np.pi*g

def radiative_correction_o2(Z, W, R):
    M = NUCLEON_MASS_C2/ELECTRON_MASS_C2
    L = np.sqrt(10)/R

    gE = 0.5772
    k2 = 0.395

    d1f = (np.log(L/M)
            -k2
            -3/(np.sqrt(10)*np.pi)*L/M
            *(0.5+gE+0.5*np.log(10)+np.log(M/L)))

    d2 = (3/(2*np.sqrt(10)*np.pi)*L/M
            *(1-np.pi/(2*np.sqrt(10))*L/M))

    d3 = (3/(np.sqrt(10)*np.pi)
            *GA
            *GM
            *L/M
            *(gE-1+0.5*np.log(10)+np.log(M/L)+np.pi/(4*np.sqrt(10)*L/M)))

    d01d4 = (np.log(M)
            -2/3*np.log(2*W)
            +35./9
            +np.pi**2/6
            -6*np.log(2))

    return Z*ALPHA**2*(d1f+d2+d3+d01d4)

def radiative_correction_o3(Z, W, W0, R):
    gE = 0.5772
    L = np.sqrt(10)/R

    a = np.pi/3-3/(2*np.pi)
    b = 4/(3*np.pi)*(11/4-gE-np.pi**2/6)
    f = np.log(2*W)-5/6
    g = 0.5*(np.log(R)**2-np.log(2*W)**2)+5/3*np.log(2*R*W)

    return Z**2*ALPHA**3*(a*np.log(L/W)+b*f+4*np.pi/3*g-0.649*np.log(2*W0))

def radiative_correction_L(W0):
    return 1.026725*(1-2*ALPHA/3/np.pi*np.log(2*W0))**2.25

def radiative_correction(Z, W, W0, R):
    o1 = radiative_correction_o1(W, W0)
    o2 = radiative_correction_o2(Z, W, R)
    o3 = radiative_correction_o3(Z, W, W0, R)
    L = radiative_correction_L(W0)

    return ((1+o1-ALPHA/2/np.pi*3*np.log(PROTON_MASS_C2/ELECTRON_MASS_C2/W0))
            *(L+o2+o3))

def radiative_correction_neutrino(Wv, W0):
    r = np.ones(Wv.shape)
    mask = (Wv<W0)

    W = W0-Wv[mask]

    p = np.sqrt(p**2-1)
    beta = p/W

    r[mask] = (1
            +ALPHA/2/np.pi
            *(3*np.log(PROTON_MASS_C2/ELECTRON_MASS_C2)+23/4
                -8/beta*spence(1-(2*beta)/(1+beta))
                +8*(np.arctanh(beta)/beta-1)*np.log(2*W*beta)
                +4*np.arctanh(beta)/beta*((7+3*beta**2)/8-2*np.arctanh(beta))))
    return r

def recoil_fermi(W, W0, A):
    M = A*NUCLEON_MASS_C2/ELECTRON_MASS_C2
    M2 = M**2

    Vr0 = W0 * W0 / 2. / M2 - 11. / 6. / M2
    Vr1 = W0 / 3. / M2
    Vr2 = 2. / M - 4. * W0 / 3. / M2
    Vr3 = 16. / 3. / M2

    result = (1
            +Vr0
            +Vr1/W
            +Vr2*W
            +Vr3*W**2)

    return result

def recoil_gamow_teller(W, W0, A):
    M = A*NUCLEON_MASS_C2/ELECTRON_MASS_C2
    M2 = M**2

    Ar0 = -2. * W0 / 3. / M - W0 * W0 / 6. / M2 - 77. / 18. / M2
    Ar1 = -2. / 3. / M + 7. * W0 / 9. / M2
    Ar2 = 10. / 3. / M - 28. * W0 / 9. / M2
    Ar3 = 88. / 9. / M2

    result = (1
            +Ar0
            +Ar1/W
            +Ar2*W
            +Ar3*W**2)

    return result

def recoil_Coulomb_fermi(Z, W, W0, A, beta_type):
    M = A*NUCLEON_MASS_C2/ELECTRON_MASS_C2
    p = np.sqrt(W**2-1)

    return (1
            -beta_type*ALPHA*Z*np.pi/M/p
            *(1+(W0-W)/(3*W)))

def recoil_Coulomb_gamow_teller(Z, W, W0, A, beta_type):
    M = A*NUCLEON_MASS_C2/ELECTRON_MASS_C2
    p = np.sqrt(W**2-1)

    return (1
            -beta_type*ALPHA*Z*np.pi/M/p
            *(1-1/3*(W0-W)/(3*W)))

def shape_factor_fermi(Z, W, W0, R):
    C0 = (-233/630*(ALPHA*Z)**2
            -1/5*(W0*R)**2
            -6/35*ALPHA*Z*W0*R)
    C1 = (-13/35*ALPHA*Z*R
            +4/15*W0*R**2)
    Cm1 = (2/14*W0*R**2
            +1/70*ALPHA*Z*R)
    C2 = -4/15*R**2

    return 1 + C0 + C1*W + Cm1/W + C2*W**2

def shape_factor_gamow_teller(Z, W, W0, R, A, b, c, d, L):
    M = A*NUCLEON_MASS_C2/ELECTRON_MASS_C2
    bt = Z/abs(Z)

    C0 = (-1/5*(W0*R)**2
            +4/9*R**2*(1-L/20)
            +1/3*W0/M/c*(-bt*b+d)
            +2/5*ALPHA*Z/(M*R*c1)*(bt*2*b+d)
            +2/35*ALPHA*Z*W0*R*(1-L)
            -233/630*(ALPHA*Z)**2)
    C1 = (bt*4/3*b*M/c
            +4/9*W0*R**2*(1-L/10)
            -4/7*ALPHA*Z*R*(1-L/10))
    Cm1 = (-1/3/M/c*(2*bt*b+d)
            -2/45*W0*R**2*(1-L)
            -ALPHA*Z*R/70)
    C2 = -4/9*R**2*(1-L/10)

    return 1 + C0 + C1*W + Cm1/W + C2*W**2

def atomic_screening(Z, W, R, l):
    beta_type = Z/abs(Z)

    S = np.ones(W.shape)
    X = np.ones(W.shape)
    mask = (W>1)

    W = W[mask]

    p = np.sqrt(W**2-1)
    Wt = W - beta_type * 0.5 * ALPHA * (abs(Z)-beta_type)*l

    pt = (0.5*p
            +0.5*np.sqrt(p**2-2*ALPHA*Z*Wt*l+0.j))

    y = ALPHA*Z*W/p
    yt = ALPHA*Z*Wt/pt
    g = np.sqrt(1-(ALPHA*Z)**2)

    S[mask] = (np.abs(gamma(g+1.j*yt))**2/np.abs(gamma(g+1.j*y))**2
           # *np.abs(np.exp(loggamma(g+2.j*pt/l)))**2/np.abs(np.exp(loggamma(1+2.j*p/l)))**2
            *np.abs(gamma(g+1.j*2*pt/l))**2/np.abs(gamma(1+1.j*2*p/l))**2
            *np.exp(-np.pi*y)
            *np.power(2*p/l, 2*(1-g)))

    X[mask] = (1/(1+0.25*(l/p)**2)
            *(1+1/8*(Wt+g)/Wt*(l/p)**2
                +0.5*g**2/(1+np.sqrt(1-ALPHA*Z*l/(W+1)))**2
                *(W-1)/Wt*(l/p)**2
                *(1-1/8*(1-g)/g*(l/p)**2)))

    return X*S


def atomic_mismatch(Z, W, W0, A):
    beta_type = Z/abs(Z)
    dBdZ2 = (44.200 * np.power(Z - beta_type, 0.41)
            +2.3196e-7 * np.power(Z - beta_type, 4.45))
    K = -0.872 + 1.270 * np.power(abs(Z), 0.097) + 9.062e-11 * np.power(abs(Z), 4.5)
    beta = np.sqrt(W**2-1)/W
    l = 1.83E-3 * K * Z / beta
    M = A * NUCLEON_MASS_C2/ELECTRON_MASS_C2
    vR = np.sqrt(1 - M * M / (M * M + (W0 * W0 - 1) / 4.))
    psi2 = 1 + 2 * ALPHA / beta * (np.arctan(1 / l) - l / 2 / (1 + l * l))
    C0 = -ALPHA**3*Z / beta * l / (1 + l * l) / psi2 #CHECK ALPHA power
    C1 = (2*ALPHA**2*Z * vR / beta
            *((0.5 + l * l) / (1 + l * l) - l * np.arctan(1 / l)) / psi2)

    return 1 - 2 / (W0 - W) * (0.5 * dBdZ2 + 2 * (C0 + C1))
