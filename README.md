# Computing Raman Spectra from MD Simulations

This is a tool for computing Raman spectra of solids or molecules based on polarizablity timeseries from on Molecular Dynamics (MD) simulations, developed by the [TheoFEM](https://theofem.de/) research group at TUM.
This repo contains the code and an example for a Raman spectrum of alpha-SiO2 at 300K.

## Usage

The script can be called in the following way:

```
./md_raman.py INPUTFILE TIMESTEP SUBSAMPLING
```

Here, `INPUTFILE` is a file that contains a polarizability timeseries (see below for format details), `TIMESTEP` is the timestep of the underlying MD simulation in femtoseconds, and `SUBSAMPLING` is the subsampling factor for the polarizability timeseries.
For example, if the polarizability was computed for every 10th timestep in an MD simulation, then the subsampling factor would be `10`.

The following additional options can be specified:

```
  --temperature         Temperature in Kelvin (default: 300).
  --omega-in            Incoming laser frequency in cm⁻¹ (default: 514nm converted to cm⁻¹).
  --smoothing           Smoothing type: gaussian, lorentzian, or none (default: lorentzian).
  --sigma SIGMA         Smoothing width in cm⁻¹ (default: 10).
  --include-prefactor   Inclusion of the factor 1/(1-e^(hbar*omega/kbT)) (default: false).
```

## Input Format

The input file should contain a timeseries of polarizability components in the following format:

```
     alpha_xx      alpha_yy      alpha_zz      alpha_xy      alpha_xz      alpha_zx
```

As an example, the first few lines of the file might look like this:

```
     +2.534301     +2.548085     +2.577736     +0.008830     -0.005494     -0.005261
     +2.533335     +2.549526     +2.574396     +0.006313     +0.000704     -0.003617
     +2.565470     +2.573721     +2.618057     -0.001479     +0.003312     -0.000895
     +2.602988     +2.614486     +2.660008     +0.002718     -0.009674     +0.007126
     +2.618783     +2.643407     +2.668517     +0.006070     -0.007723     +0.017191
     +2.578381     +2.609073     +2.644925     +0.006837     +0.016772     +0.003504
     +2.526333     +2.562766     +2.594085     +0.013118     +0.024493     -0.007330
     +2.498738     +2.532230     +2.550043     +0.014374     +0.014851     -0.004920
     +2.493625     +2.514986     +2.537445     +0.010511     +0.001420     -0.003488
     +2.517608     +2.526416     +2.555534     +0.013001     -0.012135     -0.014414
```

## Citation

D.A. Egger, M. Grumet, T. Bučko, "Machine Learning Accelerates Raman Computations from Molecular Dynamics for Materials Science", in preparation.
