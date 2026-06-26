<p align="center">
  <img width="120" src="https://github.com/user-attachments/assets/30171cc6-7d53-48a8-ab5e-027f61a6d743" />
</p>
<h1 align="center">FD-BPM</h1>
<p align="center"><em>Library originally written in 2021</em></p> 

![Language](https://img.shields.io/badge/MATLAB-R2018a%2B-orange)
![Solver](https://img.shields.io/badge/MEX-C%20%2B%20OpenMP-blue)

`FD-BPM` is a **Finite-Difference Beam Propagation Method** solver for the paraxial
wave equation in (2+1)D, tailored for **waveguide simulations**. You define a refractive-index
`n(x, y, z)` from a set of straight/curved waveguides, launch an input field, and it
marches the field along `z`, showing the intensity, phase and index live as it propagates. The
per-step linear solve runs in a C MEX (`thomas.c`) parallelised with OpenMP exploiting a natural parallelisation of this system.

If you have used ```FD-BPM``` in a scientific publication, we would appreciate citation to the following reference:
```bibtex
@article{de2023dynamics,
  title={Dynamics of the Generation of Independent Orbital-Angular-Momentum Modes in a Photonic Chip},
  author={de Oliveira, JM and Rocha, JCA and Santos, LMS and Moura, JVS and Jesus-Silva, AJ and Fonseca, EJS},
  journal={Physical Review Applied},
  volume={20},
  number={2},
  pages={024004},
  year={2023},
  publisher={APS}
}
```

## What it does

- **Propagation** marches the paraxial equation `2ik₀n₀ ∂U/∂z = ∇⊥²U + k₀²(n² − n₀²)U` along `z` using an **ADI (alternating-direction-implicit)** scheme: two tridiagonal half-steps per `z` slice, solved with the Thomas algorithm. Most of the runtime is here (the per-slice solve and building `n(x,y,z)`).
- **PML boundaries** complex-coordinate-stretched **perfectly matched layers** absorb the field at the window edges so there are minimal reflections.
- **Waveguides as functions of z** each guide is a handle `@(z) ...` returning its centre, cross-section and `Δn` at that `z`. Build straight runs ([`straightWaveguide`](functions/straightWaveguide.m)) or smooth S-bends ([`curvedWaveguide`](functions/curvedWaveguide.m)), and stack as many as you like. The cross-section is any function you pass (step-index box, Gaussian, ring, …).
- **Eigenmodes** solve for the guided modes of a given index profile via a sparse eigenproblem ([`findModes`](functions/findModes.m)).
- **Live view & video** updates the intensity/phase image every `plotStep` slices and can record an MP4 of the run.

## How the code looks

A complete run, straight single-mode guide with Gaussian launch (from [`examples/simple.m`](examples/simple.m)):

```matlab
addpath('functions')

% --- grid & optics (micrometres) ---
N = 350;  L = 80;                 % N×N window, L µm wide
n0 = 1.5078;  lambda = 640e-3;    % background index, wavelength
dz = 20;  zFinal = 5000;          % step and total propagation length
zList = linspace(0, zFinal, round(zFinal/dz));

x = -L/2 : L/(N-1) : L/2;  [X,Y] = meshgrid(x, x);
k = 2*pi/lambda * n0;

% --- input field: Gaussian, waist w0 ---
G = @(x,w0) exp(-x.^2./w0.^2);
U = createInitialField(G, X, Y, 10, k, 0, 0, [0 0]);

% --- one straight step-index waveguide, Δn = 1.2e-3 ---
rect = @(x,w) double(abs(x/w) < 1/2);
WG   = @(a,b) rect(X-a, 9.6).*rect(Y-b, 9.6);
waveguides = { @(z) curvedWaveguide(z, [0,0], [0,0], 0, zFinal, WG, 1.2e-3) };

% --- live intensity figure ---
img1 = imagesc(x, x, zeros(N)); axis image; colormap(inferno(1024));

% --- propagate ---
simParams = struct('N',N, 'L',L, 'n0',n0, 'lambda',lambda, 'plotStep',10);
U = FDpropagate(U, simParams, zList, waveguides, {img1});
```

A waveguide handle returns a struct at each `z`; chaining segments with different `zlim`
ranges builds couplers, S-bends, and index/shape transitions (see the examples below).

## The solver

`FDpropagate` is the core. Per `z`-slice it:

1. builds the index `n(x,y)` for that slice via [`applyRefIndex`](functions/applyRefIndex.m) (sum of all active waveguides),
2. assembles the two tridiagonal systems (x-sweep, then y-sweep) including the PML stretch factors,
3. solves each with [`thomas.c`](functions/thomas.c), a batched Thomas solver. Each independent line of the system is a separate tridiagonal solve, so the batch is parallelised across rows with OpenMP. Everything runs in **single-precision complex** (`float complex`) for throughput.

## Repository layout

| Path | What |
|---|---|
| [`FDsimulation.m`](FDsimulation.m) | Top-level driver script. |
| [`functions/FDpropagate.m`](functions/FDpropagate.m) | ADI/PML propagation loo. |
| [`functions/thomas.c`](functions/thomas.c) | Batched tridiagonal solver (MEX, OpenMP). |
| [`functions/applyRefIndex.m`](functions/applyRefIndex.m) | Builds `n(x,y)` at a given `z`. |
| [`functions/straightWaveguide.m`](functions/straightWaveguide.m) · [`curvedWaveguide.m`](functions/curvedWaveguide.m) | Waveguide segment generators. |
| [`functions/createInitialField.m`](functions/createInitialField.m) · [`createInitialFieldLG.m`](functions/createInitialFieldLG.m) | Input-field builders. |
| [`functions/findModes.m`](functions/findModes.m) | Eigenmode solver. |
| [`functions/visualize.m`](functions/visualize.m) | 3D geometry rendering. |
| [`functions/ASpropagate.m`](functions/ASpropagate.m) | Angular-spectrum free-space propagator. |
| [`examples/`](examples/) | Ready-to-run scenarios (see below). |

## Examples

| Script | Scenario |
|---|---|
| [`examples/simple.m`](examples/simple.m) | Single straight waveguide, Gaussian launch. |
| [`examples/manyGuides.m`](examples/manyGuides.m) | Multiple guides in one window. |
| [`examples/ycombiner.m`](examples/ycombiner.m) · [`ycombiner_ydivider.m`](examples/ycombiner_ydivider.m) | Y-junction combiner / divider. |
| [`examples/changingIndexAlong.m`](examples/changingIndexAlong.m) | `Δn` varying along `z`. |
| [`examples/changeIndexAndShape.m`](examples/changeIndexAndShape.m) | Index **and** cross-section transition. |
| [`examples/simpleLG_with_phase.m`](examples/simpleLG_with_phase.m) | Laguerre-Gaussian launch with phase plot. |

## Build

The propagation loop calls a compiled MEX. A prebuilt `thomas.mexw64` (Windows x64) is included;
rebuild it for your platform / to pick up changes. From the `functions` folder in MATLAB:

```matlab
mex -R2018a CFLAGS='$CFLAGS -fopenmp -ffast-math' LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS='-O3 -DNDEBUG' thomas.c
```

- `-R2018a` enables the interleaved-complex API (`mxGetComplexSingles`).
- `-fopenmp` must be in **both** `CFLAGS` and `LDFLAGS` so the OpenMP runtime links.
- The solver uses `omp_get_max_threads()`; control it with the `OMP_NUM_THREADS` environment variable.

## Requirements

- **MATLAB R2018a+** (interleaved-complex MEX API).
- A C compiler with **OpenMP** configured for `mex` (MinGW-w64, MSVC, or gcc/clang). Run `mex -setup C` once.
- No toolboxes required for the core solver; `findModes` uses sparse `eigs`.


## References

The solver follows the FD-BPM and PML formulations in:

- Pedrola, Ginés Lifante. *Beam Propagation Method for Design of Optical Waveguide Devices.* John Wiley & Sons, 2015.
- Huang, W. P., et al. "The perfectly matched layer (PML) boundary condition for the beam propagation method." *IEEE Photonics Technology Letters* 8.5 (1996): 649-651.

Since I couldn't find a reference which contained the finite difference equations for (2+1)D FD-ADI BPM with PML, I derived them by hand in 2020 and used in the code.

## Units

All lengths are in **micrometres (µm)**; angles in `createInitialField` are in degrees.
