#!/usr/bin/env python3

import argparse

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt


# ignore division-by-zero errors at omega=0
np.seterr(divide='ignore', invalid='ignore')


# Fourier-transformed autocorrelation via multiplication in frequency domain
def compute_autocorrelation(x):
    return np.abs(np.fft.rfft(x - np.mean(x)))**2


# Raman computation from tensor invariants
def compute_spectrum(alpha, args):
    delta_t_md = 1e-15 * args.timestep                  # seconds
    delta_t_alpha = delta_t_md * args.subsampling       # seconds
    kbT = args.temperature * sc.k                       # Joule

    # time derivatives of the polarizability components
    alpha_dot = np.diff(alpha, axis=1) / delta_t_alpha  # 1/seconds

    # frequencies with correction for systematic error in Verlet integration
    omega_hz = np.fft.rfftfreq(len(alpha_dot[0]), delta_t_alpha)                # Hz
    omega_hz = 1 / delta_t_md * np.acos(1 - 1/2 * omega_hz**2 * delta_t_md**2)  # Hz
    omega = omega_hz / (100 * sc.c)                                             # cm⁻¹

    # terms for the tensor invariants a_p and gamma_p (with relevant prefactors)
    terms = np.array([
      (alpha_dot[0] + alpha_dot[1] + alpha_dot[2])/3 * 45.0**0.5,
      (alpha_dot[0] - alpha_dot[1]) * 3.5**0.5,
      (alpha_dot[1] - alpha_dot[2]) * 3.5**0.5,
      (alpha_dot[2] - alpha_dot[0]) * 3.5**0.5,
      alpha_dot[3] * 21.0**0.5,
      alpha_dot[4] * 21.0**0.5,
      alpha_dot[5] * 21.0**0.5,
    ])

    # Fourier-transformed autocorrelations of these terms
    autocorrelations = list(compute_autocorrelation(t) for t in terms)

    # computation of the Raman spectrum
    spectrum = np.sum(autocorrelations, axis=0)
    spectrum *= 1/45 * (args.omega_in - omega)**4 / omega

    # computation of the prefactor
    if args.include_prefactor:
        spectrum *= 1 / (1 - np.exp(-sc.hbar*omega_hz/kbT))

    # correction of the spectrum for systematic error in numerical differentiation
    spectrum *= (np.sin(omega_hz * delta_t_alpha) / (omega_hz * delta_t_alpha))**2

    return omega, spectrum


# Lorentzian or Gaussian smoothing
def smoothen_spectrum(omega, spectrum, smoothing, sigma):
    if smoothing == "none" or sigma == 0:
        return omega, spectrum

    omega_new = np.arange(0, np.max(omega), sigma / 10)
    spectrum_new = np.zeros(len(omega_new))

    for i in range(1, len(omega)):
        if smoothing == "gaussian":
            spectrum_new += spectrum[i] * np.exp(-(omega[i] - omega_new)**2 / (sigma**2))
        elif smoothing == "lorentzian":
            spectrum_new += spectrum[i] / ((omega[i] - omega_new)**2 + sigma**2)
        else:
            raise ValueError(smoothing)

    return omega_new, spectrum_new


if __name__ == "__main__":
    parser = argparse.ArgumentParser("./md_raman.py")
    parser.add_argument('inputfile', type=str)
    parser.add_argument('timestep', type=float)                         # femtoseconds
    parser.add_argument('subsampling', type=float)                      # dimensionless
    parser.add_argument('--temperature', type=float, default=300)       # Kelvin
    parser.add_argument('--omega-in', type=float, default=100/514e-9)   # cm⁻¹
    parser.add_argument('--smoothing', type=str, default="lorentzian")
    parser.add_argument('--sigma', type=float, default=10)              # cm⁻¹
    parser.add_argument('--include-prefactor', action="store_true")
    args = parser.parse_args()

    alpha = np.loadtxt(args.inputfile, unpack=True)

    omega, spectrum = compute_spectrum(alpha, args)
    omega, spectrum = smoothen_spectrum(omega, spectrum, args.smoothing, args.sigma)

    fig, ax = plt.subplots(1, 1)
    ax.plot(omega, spectrum)
    ax.set_xlabel(r"$\omega \, \mathrm{[cm^{-1}]}$")
    ax.set_ylabel(r"$I(\omega) \, \mathrm{[AU]}$")
    ax.set_xlim(xmin=0, xmax=omega[-1])
    ax.set_ylim(ymin=0)
    ax.set_yticks([])
    ax.minorticks_on()
    plt.tight_layout()
    plt.savefig("raman_spectrum.png", dpi=150)

    np.savetxt("raman_spectrum.dat", np.array([omega, spectrum]).T)
