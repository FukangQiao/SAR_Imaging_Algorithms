# SAR Image Formation — Range-Doppler Algorithm (RDA)

> MATLAB implementation of a stripmap Synthetic Aperture Radar (SAR) simulator with full Range-Doppler Algorithm processing chain: echo generation → range compression → RCMC → azimuth compression.

![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange?style=flat-square&logo=mathworks)
![Signal Processing](https://img.shields.io/badge/Domain-SAR%20%7C%20Signal%20Processing-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-lightgrey?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)

---

## Overview

This project simulates a **stripmap SAR system** and implements the full **Range-Doppler Algorithm (RDA)** image formation pipeline in MATLAB. Three point targets are placed in the scene and the system simulates raw echo generation, followed by a complete two-dimensional focusing chain.

The simulation models a squint-capable platform flying at constant velocity and altitude, transmitting a **Linear Frequency Modulated (LFM)** waveform. The output is a focused 2D SAR image of the target scene.

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

| Domain | Points | Description |
|---|---|---|
| Azimuth (u) | Na = 1024 | Cross-range samples |
| Range (t) | Nr = 512 | Fast-time samples |

---

## Target Configuration

Three point targets are placed in the scene:

| Target | Range (m) | Azimuth (m) |
|---|---|---|
| T1 | 10000 | 0 |
| T2 | 9950 | 0 |
| T3 | 10000 | 50 |

---

## Processing Chain

```
Raw Echo Generation
        │
        ▼
  Hamming Windowing
  (Range + Azimuth)
        │
        ▼
  Range Compression
  (Matched Filter in f domain)
        │
        ▼
  Azimuth FFT
  → Range-Doppler Domain
        │
        ▼
  Range Cell Migration
  Correction (RCMC)
        │
        ▼
  Azimuth Compression
  (Phase Multiply in fu domain)
        │
        ▼
  2D IFFT → Focused SAR Image
```

### Step-by-step description

**1. Echo generation** — For each target at `(rn, xn)`, the instantaneous slant range `R(u)` is computed as a function of azimuth slow-time `u`. The raw echo is modeled as a 2D LFM signal with quadratic phase in both range and azimuth, windowed to the synthetic aperture length `Lsar`.

**2. Windowing** — Hamming windows are applied in both range (over the pulse duration `Tr`) and azimuth (over `Na` samples) to suppress sidelobes.

**3. Range compression** — A reference LFM chirp centered at `Rc` is generated, transformed to the frequency domain, and used as a matched filter. The range-compressed signal `src_ut` is obtained via element-wise conjugate multiplication followed by IFFT.

**4. Range-Doppler domain** — An azimuth FFT transforms the range-compressed signal to the Range-Doppler (2D frequency) domain `src_fut`.

**5. Range Cell Migration Correction (RCMC)** — The RCMC filter `Hrcc` corrects the range walk of each target as a function of Doppler frequency `fu` and range frequency `f`, using the Doppler chirp rate `fdr`.

**6. Azimuth compression** — A matched filter in the Doppler frequency domain `fu` is applied using the quadratic phase `exp(jπ/fdr · (fu − fdc)²)` to focus targets in azimuth.

**7. Final image** — An azimuth IFFT produces the focused complex SAR image `s2rcac_ut`.

---

## Key Signal Processing Parameters

```matlab
fdc  = 2·v·sin(θT)/λ          % Doppler centroid frequency
fdr  = −2·(v·cos(θT))²/λ/Rc  % Doppler chirp rate
Lsar = λ·Rc/D                 % Synthetic aperture length
Rb   = sqrt(H² + Y²)          % Slant range to scene center
```

---

## Output Figures

| Figure | Title | Description |
|---|---|---|
| Figure 1 | Range Profile | Normalized amplitude (dB) vs. ground range after range compression — center azimuth slice |
| Figure 2 | Focused SAR Image | 2D grayscale image of focused targets in (Azimuth × Range) |
| Figure 3(a) | Range-Doppler Domain | Signal before RCMC |
| Figure 3(b) | Corrected RD Domain | Signal after RCMC |
| Figure 4 | 3D Waterfall Plot | `waterfall` visualization of target amplitudes in (Range × Azimuth) |

---

## Repository Structure

```
SAR-RDA/
│
├── sar_rda.m        # Main MATLAB script — full RDA pipeline
└── README.md
```

---

## Requirements

- MATLAB R2020a or later
- Signal Processing Toolbox (for `hamming`)
- No additional toolboxes required

---

## Getting Started

**1. Clone the repository**

```bash
git clone https://github.com/Abwahab55/SAR-RDA.git
cd SAR-RDA
```

**2. Open MATLAB and run**

```matlab
sar_rda.m
```

All four figures will be generated automatically.

---

## Customization

| What to change | Variable | Location |
|---|---|---|
| Squint angle | `thetaT` | Line 1 |
| Carrier frequency | `fc` | Line 4 |
| Platform altitude | `H` | Line 10 |
| Platform velocity | `v` | Line 18 |
| Pulse bandwidth | `Br` | Line 32 |
| Target positions | `Ptar` | Line 40 |
| Azimuth / range samples | `Na`, `Nr` | Lines 22, 31 |

To **disable windowing**, comment out the respective lines:

```matlab
% s_ut=s_ut.*(wr*ones(1,Na));   % disable range window
% s_ut=s_ut.*(ones(Nr,1)*wa');  % disable azimuth window
```

---

## Algorithm Reference

The Range-Doppler Algorithm (RDA) is a classical SAR focusing algorithm. The key processing steps follow the formulation in:

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

---

