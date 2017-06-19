import numpy as np 


karman = 0.41  # Von Karman constant

def manning(wl, n=0.03):
    """
    Manning formulas

    Return:
        Chézy coefficient
    """
    return wl**(1 / 6.) / float(n)

def yangchoi(wd, hv, Cd, n):
    """
    Formula of Yang & Choi (2010) from Augustijn et al (2011)
    """
    g = 9.81
    return np.sqrt(2 * g / (Cd * n * hv)) + 2 * np.sqrt(g * (wd - hv) / wd) / karman * (np.log(wd / hv) - (wd - hv) / wd)


def baptist1(wd, hv, Cd, n, nm=0.03):
    """
    Formula of Baptist (2005)
    """
    
    g = 9.81  # Gravitational acc.
    for h in list(wd):
        Cb = manning(h, n=nm)
        if h > hv:
            yield (np.sqrt(Cb**-2 + Cd * n * hv * (2 * g)**-1))**-1 + np.sqrt(g) / karman * np.log(h / hv)
        else:
            yield (np.sqrt(Cb**-2 + Cd * n * hv * (2 * g)**-1))**-1

