import numpy as np

H0 = 67.6
Omega_m = 0.31
Omega_Lambda = 0.0

def redshift_to_dist(z, type="DCMR", h=H0/100.0, Omega_m=Omega_m, n=1000):
    Omega_r = 4.165E-5 / (h**2)
    Omega_total = Omega_m + Omega_Lambda + Omega_r
    Omega_k = 1 - Omega_total
    a = 1.0
    Tyr = 977.8
    c = 299792.458

    if type == "DCMR":
        az = 1.0 / (1 + 1.0 * z)
        age = 0
        n = 1000
        
        for i in range(n):
            a = az * (i + 0.5) / n
            adot = np.sqrt(Omega_k + (Omega_m / a) + (Omega_r / (a**2)) + (Omega_Lambda*a**2))
            age = age + 1 / adot

        zage = az * age / n
        zage_Gyr = (Tyr / H0) * zage
        DTT = 0.0
        DCMR = 0.0

        for i in range(n):
            a = az + (1 - az) * (i + 0.5) / n
            adot = np.sqrt(Omega_k + (Omega_m / a) + (Omega_r / (a**2)) + (Omega_Lambda*a**2))
            DTT = DTT + 1. / adot
            DCMR = DCMR + 1. / (a * adot)

        DTT = (1 - az) * DTT / n
        DCMR = (1 - az) * DCMR / n
        age = DTT + zage
        age_Gyr = age * (Tyr / H0)
        DTT_Gyr = (Tyr / H0) * DTT
        DCMR_Gyr = (Tyr / H0) * DCMR
        DCMR_Mpc = (c / H0) * DCMR
        return DCMR_Mpc

    return 0