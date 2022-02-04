from context import bsgpy
import bsgpy.spectral_functions as sf

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    W0 = 3
    Z = 10
    R = 0.01
    A = 20
    betaType = 1

    W = np.linspace(1, W0, 1000)

    ps = sf.phase_space(W, W0)
    f = sf.fermi_function(Z, W, R, betaType)
    g = sf.sirlin_g(W, W0)
    rf = sf.recoil_fermi(W, W0, A)
    rc1 = (1+sf.ALPHA/2/np.pi*g)

    plt.figure()
    plt.plot(W, ps)
    plt.plot(W, ps*f)
    plt.plot(W, ps*f*rc1)
    plt.plot(W, ps*f*rc1*rf)

    plt.figure()
    plt.plot(W, f)

    plt.figure()
    plt.plot(W, rc1)

    plt.figure()
    plt.plot(W, rf)

    plt.show(block=True)
