from context import bsgpy
import bsgpy.spectral_functions as sf

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    W0 = 3
    Z = 10
    R = 0.01
    A = 20
    W = np.linspace(1, W0, 1000)
    l = 2/137*Z**0.33

    ps = sf.phase_space(W, W0)
    f = sf.fermi_function(Z, W, R)
    g = sf.sirlin_g(W, W0)
    rf = sf.recoil_fermi(W, W0, A)
    rc1 = (1+sf.ALPHA/2/np.pi*g)
    rc = sf.radiative_correction(Z, W, W0, R)
    print(l)
    s = sf.atomic_screening(Z, W, R, l)

    plt.figure()
    plt.plot(W, ps)
    plt.plot(W, ps*f)
    plt.plot(W, ps*f*rc1)
    plt.plot(W, ps*f*rc1*rf)
    plt.plot(W, ps*f*rc*rf*s)

    plt.figure()
    plt.plot(W, f)

    plt.figure()
    plt.plot(W, s)

    plt.figure()
    plt.plot(W, rc1)
    plt.plot(W, rc)

    plt.figure()
    plt.plot(W, rf)

    plt.show(block=True)
