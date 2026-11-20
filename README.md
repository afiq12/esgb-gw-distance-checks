# esgb-gw-distance-checks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24403951.svg)](https://doi.org/10.5281/zenodo.24403951)

Symbolic and numerical checks of the gravitational-wave luminosity distance on the luminal branch of Einstein-scalar-Gauss-Bonnet (EsGB) cosmology.

This repository holds verification scripts and derivation notes. It is not a research article: it documents that a set of known and derived relations has been checked independently, so that later work can build on them.

## What is checked

Conventions: `M_pl = 1`, action `S = ∫ d^4x sqrt(-g) [ R/2 - (1/2)(∂φ)^2 - V(φ) + f(φ) G ] + S_m`, flat FLRW background.

| Statement | Script | Method |
| --- | --- | --- |
| Tensor quadratic action: `G_T = 1 + 8 H f'`, `F_T = 1 + 8 f''` (dots as primes) | `scripts/check_tensor_action.py` | Direct curvature expansion to second order, Euler-Lagrange equation |
| The same `G_T`, `F_T` from the Horndeski map | `scripts/check_background_and_horndeski_map.py` | Kobayashi, Yamaguchi, Yokoyama (2011) functions |
| Luminal branch: `c_T = 1` if and only if `df/dt = c a`, giving `G_T = 1 + 8 c a H` | `scripts/check_tensor_action.py` | Symbolic substitution |
| Hamiltonian constraint `3 H^2 (1 + 8 H df/dt) = ρ_tot` | `scripts/check_background_and_horndeski_map.py` | Minisuperspace action with lapse |
| Distance ratio `d_GW/d_EM = sqrt(G_T(0)/G_T(z))` | `scripts/check_distance_ratio_numerical.py` | Numerical integration of the wave equation, compared with the WKB result |
| How well the `(Xi_0, n)` and `c_M` templates reproduce that ratio | `scripts/fit_standard_templates.py` | Least-squares fits on a LambdaCDM background |

The distance ratio on the luminal branch is the result of Linder (arXiv:2108.11526); the tensor sector of EsGB with `c_T = 1` is also analysed by Fier et al. (arXiv:2503.01975). The scripts here reproduce those relations independently.

## Results in brief

- The numerical integration agrees with the analytic ratio to a relative difference of at most `1.3e-4` for `beta = 8 c H_0` between 0.1 and 0.6 and source redshifts between 0.3 and 3.
- With `beta > 0` the ratio is above unity for `z < 2.09` and peaks at the deceleration to acceleration transition (`z ≈ 0.67` for `Omega_m = 0.3`), where the anomaly is about `0.066 beta`.
- For `z < 1` the `(Xi_0, n)` template reproduces the ratio to better than 0.3 percent, and the `c_M` template to better than 0.75 percent. When the fit is extended to `z = 2` or `3` the residuals become comparable to the signal.

## Limitations

- All numerical work uses a LambdaCDM background as a test approximation, not a self-consistent EsGB background.
- The stability conditions of the scalar sector are not derived here.
- The order-of-magnitude bounds in the notes use assumed tolerances and are not the result of a data analysis.

## How to run

```
pip install -r requirements.txt
bash run_all.sh
```

The output should match `expected_output.txt`. The numerical check takes about twenty seconds; the others take a few seconds each. Tested with Python 3.13, SymPy 1.14, SciPy 1.18, NumPy 2.5 and Matplotlib 3.11.

## Contents

```
scripts/     five Python scripts (four checks and one figure)
notes/       derivation notes (PDF and LaTeX source, in Indonesian)
figures/     ratio_gb_vs_nmc.png
expected_output.txt
```

## Citation

Agung, A. (2026). *esgb-gw-distance-checks: symbolic and numerical checks of the gravitational-wave luminosity distance on the luminal branch of Einstein-scalar-Gauss-Bonnet cosmology* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24403952

To cite all versions, use the concept DOI https://doi.org/10.5281/zenodo.24403951. See also `CITATION.cff`.

## Licence

Code: MIT. Notes and figure: CC BY 4.0.

## Author

Afiq Agung, Physics, Program Studi Pendidikan Fisika, Universitas Sulawesi Barat, Majene, West Sulawesi, Indonesia. ORCID 0000-0002-8698-9459.
