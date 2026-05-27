# SAV Solvers for Nonlinear Stiff String Models
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20413144.svg)](https://doi.org/10.5281/zenodo.20413144)

MATLAB implementation of **Scalar Auxiliary Variable (SAV)** schemes for three geometrically nonlinear stiff string models, using Finite-Difference Time-Domain (FDTD) spatial discretisation. All solvers are fully explicit (Sherman–Morrison rank-1 update) and exactly conserve a discrete energy at every time step.

> R. Russo, M. Ducceschi, S. Bilbao — *"Numerical convergence of the Scalar Auxiliary Variable method applied to nonlinear stiff string models"*, Journal of Nonlinear Science, 2025.

---

## Quick Start

No toolboxes required. Open MATLAB, navigate to the relevant folder and run the script directly:

```matlab
cd GEstring
run('GEStringSAV_FB_split.m')
```

Edit the flags and parameters at the top of each file (see [Flags](#flags) below).

---

## Repository Structure

```
.
├── GEstring/               % Geometrically exact string (u and v)
│   ├── GEStringSAV_FA_noSplit.m
│   ├── GEStringSAV_FA_split.m
│   └── GEStringSAV_FB_split.m
├── CUBstring/              % Cubic string (u only)
│   ├── CubicStringSAV_noSplit.m
│   └── CubicStringSAV_split.m
└── KCstring/               % Kirchhoff–Carrier string (u only)
    ├── KCStringSAV_noSplit.m
    └── KCStringSAV_split.m
```

---

## Models and Files

Each model is available in a **non-split (NS)** and one or more **split (S)** variants.

- **Non-split** (`SAV = true`): the full potential energy `V` is absorbed into the SAV auxiliary variable `psi`. The scheme is unconditionally stable.
- **Split** (`SAV = true`): the potential is decomposed as `V = Q + R`. The quadratic part `Q` is treated with the classical explicit FDTD scheme; only the nonlinear residual `R` enters `psi`. A CFL-type stability condition applies.

In both cases, the update is explicit and requires only a rank-1 linear solve at each time step via the Sherman–Morrison formula.

Each script also contains a reference integrator (Störmer–Verlet for the GE string; an energy-preserving explicit scheme for the cubic and KC strings), activated by setting `SAV = false`.

### Geometrically Exact String — `GEstring/`

Simulates both transverse displacement `u(x,t)` and longitudinal displacement `v(x,t)`. Two decompositions of the potential are provided.

| File | Variant | Form | Key property |
|---|---|---|---|
| `GEStringSAV_FA_noSplit.m` | NS | A | Full `V_G` in `psi`; unconditionally stable; requires `h ≥ h_MAX` for convergence |
| `GEStringSAV_FA_split.m` | S | A | Residual `R_GA` in `psi`; incomplete linear/nonlinear split; accuracy sensitive to `shiftV` |
| `GEStringSAV_FB_split.m` | S | B | Residual `R_GB` in `psi`; complete linear/nonlinear split; Jacobian bounded; insensitive to `shiftV` |

**Form A** — natural Q/R split:

$$Q_{GA} = \frac{T_0}{2}(\zeta^2+\xi^2) + \frac{EI}{2}\chi^2, \qquad R_{GA} = \frac{EA-T_0}{2}\!\left(\sqrt{(1+\xi)^2+\zeta^2}-1\right)^{\!2}$$

**Form B** — quadratic part aligned with the linearised wave equation (most robust convergence):

$$Q_{GB} = \frac{T_0}{2}\zeta^2 + \frac{EA}{2}\xi^2 + \frac{EI}{2}\chi^2, \qquad R_{GB} = (EA-T_0)\!\left(\frac{\zeta^2}{2}+\frac{3}{2}+\xi-\sqrt{(1+\xi)^2+\zeta^2}\right)$$

> ⚠️ Form B requires `v(x,t) > -1` (no string folding) for the residual to remain non-negative.

The split schemes require `h ≥ h_MAX = max(h_L, h_T)`, where:

$$h_L = \sqrt{E/\rho}\; k \qquad h_T = \sqrt{\frac{EAk^2+\sqrt{(EAk^2)^2+16\rho AEIk^2}}{2\rho A}}$$

The `hTypeLong` flag selects which condition drives the grid spacing (see [Flags](#flags)).

---

### Cubic String — `CUBstring/`

Transverse displacement `u(x,t)` only. Obtained by Taylor expansion of the geometrically exact potential and neglecting longitudinal motion.

$$Q_C = \frac{T_0}{2}\zeta^2 + \frac{EI}{2}\chi^2, \qquad R_C = \frac{EA-T_0}{8}\zeta^4$$

| File | Variant | Key property |
|---|---|---|
| `CubicStringSAV_noSplit.m` | NS | Full `V_C` in `psi`; unconditionally stable; requires `h ≥ h_T` for convergence |
| `CubicStringSAV_split.m` | S | Residual `R_C` in `psi`; complete linear/nonlinear split; Jacobian bounded at equilibrium |

---

### Kirchhoff–Carrier String — `KCstring/`

Transverse displacement `u(x,t)` only. Obtained by averaging the longitudinal strain over the string length; the nonlinear term is spatially non-local.

$$Q_K = \frac{T_0}{2}\zeta^2 + \frac{EI}{2}\chi^2, \qquad R_K = \frac{EA}{8L}\,\zeta^2\int_0^L\zeta^2\,\mathrm{d}x$$

| File | Variant | Key property |
|---|---|---|
| `KCStringSAV_noSplit.m` | NS | Full `V_K` in `psi`; unconditionally stable; requires `h ≥ h_T` for convergence |
| `KCStringSAV_split.m` | S | Jacobian bounded independently of `shiftV`; no shift constant needed |

---

## Flags

The following variables are set near the top of every script under the `% Flags` and `% Custom Parameters` sections.

### Solver selection

| Variable | Type | Description |
|---|---|---|
| `SAV` | `true`/`false` | `true` runs the SAV scheme; `false` runs the built-in reference integrator (Störmer–Verlet for GE string; energy-preserving explicit scheme for cubic and KC) |
| `nLinOn` | `true`/`false` | *(cubic and KC only)* Enables/disables the nonlinear term. Useful for sanity checks against the linear stiff string |
| `hTypeLong` | `true`/`false` | *(GE string only)* Selects which stability condition drives the grid: `false` (default) uses the transverse condition `h_T`; `true` uses the longitudinal condition `h_L` |

### Simulation parameters

| Variable | Description |
|---|---|
| `OSfac` | Oversampling factor relative to 44100 Hz. Effective sample rate is `SR = OSfac * 44100`. |
| `durSec` | Simulation duration in seconds (default: `0.06`) |
| `amplitude` | Peak initial displacement in metres (default: `2e-3`) |
| `initType` | Initial condition: `1` = first vibration mode (sinusoid); `2` = raised cosine centred at the string midpoint, width = 25% of string length |
| `small_shift` | Selects the SAV shift constant `shiftV`: `false` (default) sets `shiftV = 1e3`; `true` sets `shiftV = eps` (machine epsilon, ≈ 2.2×10⁻¹⁶). See note below |
| `dampOn` | `true`/`false` — enables frequency-dependent damping (loss parameters `sigma0`, `sigma1`) |

> **Note on `shiftV`:** This is the constant `ε` in `ψ = √(2Φ + ε)`. A large value (`1e3`) improves the accuracy of the auxiliary variable `psi` in schemes where the SAV potential vanishes near equilibrium (NS and Form A split). It has no effect on Form B and KC-split schemes, where the Jacobian of `g` is bounded independently of `ε`.

### Output and display

| Variable | Description |
|---|---|
| `plotPsi` | Plots the auxiliary variable `psi` over time |
| `computeEnergy` | Computes and plots the discrete energy variation `ΔH` at each time step (should be at machine-epsilon level) |
| `plotSpect` | Plots the spectrogram of the output displacement |
| `realTimeDraw` | Animates the string shape during the simulation (slows down execution significantly) |
| `play` | Plays the output signal through the audio device at the end of the simulation |
| `gridOn` | Overlays horizontal lines at multiples of `eps` on the energy plot |

---

## Physical Parameters

Default values match the steel string used in the paper:

| Variable | Quantity | Value |
|---|---|---|
| `rho` | Density | 8050 kg/m³ |
| `radius` | String radius (0.016 gauge) | 3.556×10⁻⁴ m |
| `T0` | Rest tension | 75 N |
| `E` | Young's modulus | 174 GPa |
| `L` | String length | 1 m |

`Area`, `I`, `rA`, `K`, `c` are derived automatically from the above.

---

## Citation

```bibtex
@article{russo2025sav,
  title   = {Numerical convergence of the {S}calar {A}uxiliary {V}ariable method
             applied to nonlinear stiff string models},
  author  = {Russo, Riccardo and Ducceschi, Michele and Bilbao, Stefan},
  journal = {Nonlinear Dynamics},
  year    = {2026}
}
```

---

## Funding

European Union Horizon 2020 — grant **NEMUS-StG-950084**.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
