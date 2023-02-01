from context import bsg
import bsg.spectral_functions as sf

if __name__ == "__main__":
    import numpy as np

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
    s = sf.atomic_screening(Z, W, R, l)

    E = (W-1)*0.511

    try:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(E, ps)
        plt.plot(E, ps*f)
        plt.plot(E, ps*f*rc1)
        plt.plot(E, ps*f*rc1*rf)
        plt.plot(E, ps*f*rc*rf*s)
        plt.xlabel("Electron energy [MeV]")
        plt.ylabel("Spectrum [arb.u.]")

        plt.figure()
#    plt.plot(E, f, label="Fermi function")
        plt.plot(E, s, label="Screening")
        plt.plot(E, rc, label=r"Radiative correction O($\alpha^3Z^2$)")
        plt.legend(loc=0)
        plt.xlabel("Electron energy [MeV]")
        plt.ylabel("Spectral correction factor")

        plt.figure()
        plt.plot(W, s)

        plt.figure()
        plt.plot(W, rc1)
        plt.plot(W, rc)

        plt.figure()
        plt.plot(W, rf)

        plt.show(block=True)
    except:
        pass

    print("Done!")
