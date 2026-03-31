# Synthetic Aperture Radar (SAR) Simulation

> MATLAB implementation of a stripmap SAR imaging system with full Range-Doppler Algorithm (RDA) processing chain — raw echo generation, range compression, RCMC, and azimuth compression — producing focused 2D and 3D target images.

![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange?style=flat-square&logo=mathworks)
![Domain](https://img.shields.io/badge/Domain-Radar%20%7C%20Signal%20Processing-blue?style=flat-square)
![Algorithm](https://img.shields.io/badge/Algorithm-Range--Doppler%20(RDA)-purple?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-lightgrey?style=flat-square)
![Stars](https://img.shields.io/github/stars/Abwahab55/Synthetic-Aperature-radar-SAR-simulation?style=flat-square)

---

## Overview

This project simulates a **stripmap Synthetic Aperture Radar (SAR)** system in MATLAB, implementing the complete **Range-Doppler Algorithm (RDA)** image formation pipeline. Three point targets are placed in the scene and the system generates raw radar echoes, then applies a full two-dimensional focusing chain to produce a sharp SAR image.

The simulation models a squint-capable airborne platform flying at constant altitude and velocity, transmitting a **Linear Frequency Modulated (LFM)** waveform and receiving the scattered echoes. The output is a focused 2D SAR image and 3D waterfall visualization of the target scene.

---

## SAR System Parameters

| Parameter | Symbol | Value |
|---|---|---|
| Platform height | H | 4000 m |
| Platform velocity | v | 100 m/s |
| Squint angle | θT | 0° (broadside) |
| Carrier frequency | fc | 1.5 GHz |
| Wavelength | λ | 0.2 m |
| Antenna aperture | D | 4 m |
| Scene center range | Y | 10000 m |
| Swath half-width | R0 | 500 m |
| Azimuth scene extent | X0 | ±400 m |

### LFM Pulse Parameters

| Parameter | Symbol | Value |
|---|---|---|
| Pulse duration | Tr | 5 μs |
| Bandwidth | Br | 50 MHz |
| Chirp rate | Kr | 10¹³ Hz/s |

### Processing Grid

| Domain | Samples | Description |
|---|---|---|
| Azimuth (u) | Na = 1024 | Cross-range samples |
| Range (t) | Nr = 512 | Fast-time samples |

---

## Target Configuration

Three point targets placed in the scene:

| Target | Range (m) | Azimuth (m) | Notes |
|---|---|---|---|
| T1 | 10000 | 0 | Scene center |
| T2 | 9950 | 0 | 50 m closer in range |
| T3 | 10000 | 50 | 50 m offset in azimuth |

The targets are designed to test range and azimuth resolution independently.

---

## Processing Chain

```
Raw Echo Generation
  │  (LFM phase model, 3 targets, windowed to synthetic aperture)
  ▼
Hamming Windowing
  │  (range: over pulse Tr) + (azimuth: over Na samples)
  ▼
Range Compression
  │  (matched filter in frequency domain via FFT/IFFT)
  ▼
Azimuth FFT → Range-Doppler Domain
  │
  ▼
Range Cell Migration Correction (RCMC)
  │  (Hrcc filter: compensates range walk vs. Doppler frequency)
  ▼
Azimuth Compression
  │  (quadratic phase matched filter in Doppler frequency domain)
  ▼
2D IFFT → Focused SAR Image
```

### Key Signal Processing Equations

```
Doppler centroid:    fdc  =  2·v·sin(θT) / λ
Doppler chirp rate:  fdr  = −2·(v·cos(θT))² / (λ·Rc)
Synthetic aperture:  Lsar =  λ·Rc / D
Slant range:         Rb   =  √(H² + Y²)
RCMC filter:         Hrcc = exp(−jπ·fu²·f / (fc·fdr))
Azimuth filter:      p0   = exp( jπ·(fu − fdc)² / fdr)
```

---

## Output Figures

| Figure | Description |
|---|---|
| Figure 1 | Range profile — normalized amplitude (dB) vs. ground range at center azimuth after range compression |
| Figure 2 | Focused 2D SAR image — grayscale intensity map of all three targets in (Azimuth × Range) |
| Figure 3(a) | Range-Doppler domain — signal before RCMC |
| Figure 3(b) | Corrected RD domain — signal after RCMC showing migration correction |
| Figure 4 | 3D waterfall plot — `waterfall()` visualization of target amplitudes across (Range × Azimuth) |

---

## Repository Structure

```
Synthetic-Aperature-radar-SAR-simulation/
│
├── stripmapSAR.m    # Complete SAR simulation — echo generation + full RDA pipeline
└── README.md
```

---

## Requirements

- MATLAB R2020a or later
- Signal Processing Toolbox (for `hamming()`)
- No additional toolboxes required

---

## Getting Started

**1. Clone the repository**

```bash
git clone https://github.com/Abwahab55/Synthetic-Aperature-radar-SAR-simulation.git
cd Synthetic-Aperature-radar-SAR-simulation
```

**2. Open MATLAB and run**

```matlab
run('stripmapSAR.m')
```

All four figures generate automatically. No configuration needed.

---

## Customization Guide

| What to change | Variable | Notes |
|---|---|---|
| Squint angle | `thetaT` | 0° = broadside, increase for squint mode |
| Carrier frequency | `fc` | Affects λ and Doppler parameters |
| Platform altitude | `H` | Changes slant range geometry |
| Platform velocity | `v` | Affects synthetic aperture length |
| Pulse bandwidth | `Br` | Controls range resolution |
| Pulse duration | `Tr` | Affects range ambiguity |
| Target positions | `Ptar` | Add rows for more targets |
| Azimuth samples | `Na` | Increase for finer azimuth resolution |
| Range samples | `Nr` | Increase for wider swath |

**To disable windowing** (observe sidelobe effects):

```matlab
% s_ut = s_ut .* (wr * ones(1, Na));   % disable range window
% s_ut = s_ut .* (ones(Nr,1) * wa');   % disable azimuth window
```

**To test squint mode**, change:

```matlab
thetaT = 15;  % 15-degree squint angle
```

---

## Background

Synthetic Aperture Radar is an active microwave imaging technique that synthesizes a large virtual antenna aperture by coherently combining echoes received along a flight path, achieving fine azimuth resolution independent of range. The **Range-Doppler Algorithm** is the classical two-step image formation approach:

1. **Range compression** — matched filtering in fast-time to achieve range resolution
2. **Azimuth compression** — matched filtering in slow-time (Doppler domain) to achieve cross-range resolution, after correcting for the range cell migration caused by the curved flight path geometry

The RDA is foundational in SAR signal processing and is the basis for more advanced algorithms including the Chirp Scaling Algorithm (CSA) and ω-k (wavenumber domain) algorithm.

---

## Relation to Published Work

This simulation was developed as part of research in radar signal processing and bistatic SAR systems:

> Abdul Wahab, Kashif Hussain, Haroon Ahmed, Tahir Bashir, Kamil Jalal & Zeeshan Hameed, **"A Novel Bistatic-SAR Simulation-Based On Fixed Receiver"**, *2020 IEEE 23rd International Multitopic Conference (INMIC)*, 2020.

---

## Reference

> I. G. Cumming and F. H. Wong, *Digital Processing of Synthetic Aperture Radar Data: Algorithms and Implementation*. Norwood, MA: Artech House, 2005.

---

## Author

**Abdul Wahab**
AE @ Lumissil Microsystems | SiC power systems → Cloud computing 

- GitHub: [@Abwahab55](https://github.com/Abwahab55)
- Email: wahab.engr55@yahoo.com

---
## License

This project is licensed under the MIT License - see the LICENSE file for details.
Copyright (c) 2026 Abdul Wahab

