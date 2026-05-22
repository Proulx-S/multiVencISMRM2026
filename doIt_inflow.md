# Inflow model — section 13 documentation

**Branch:** `dev-inflow`  
**Entry point:** `doIt_inflow.m`  
**Branched from:** `dev-postISMRM`

---

## TL;DR

Section 13 fits a physics-based inflow model to phantom PC-MRI data. Instead of treating the magnitude–velocity relationship $m(v)$ as a free polynomial, it derives $m(v)$ from steady-state spoiled GRE theory (Bianciardi et al. 2016) via `getMz_ss` / `getMxy_ss`. The vessel cross-section is an elliptical parabolic flow profile; the vessel tilt is parametrized by a Cartesian unit-vector projection $(t_x, t_y)$ that simultaneously drives both the ellipse aspect ratio and the effective slab thickness. A three-compartment model (lumen / wall / tissue) is fit jointly to magnitude, velocity, and optionally noFlow images using `lsqnonlin`. Pre-realized complex noise replaces the usual Rician approximation, allowing wall Mxy to be fixed at zero while still predicting nonzero wall signal. Two sequential fits are run: fit 13a with the noFlow constraint, fit 13b warm-started from 13a without it.

---

## Environment and data

The script auto-detects the host (`takoyaki` or laptop) and mounts SSHFS if needed. All tool repos (`util`, `pcMRAsim`, `multiVencSim`) are cloned to `dev-inflow` branches. `pcMRAsim` is re-added to the MATLAB path **after** `multiVencSim` because `multiVencSim` contains a legacy `runSim.m` that shadows `pcMRAsim`'s:

```matlab
addpath(genpath(fullfile(toolDir, 'pcMRAsim')));
```

Phantom data (session `20251010_multiVENCphantom03`) is loaded via `loadPhantom03`. Outputs: `data` (complex, all VENCs), `dataNoFlow`, `PEspacing`, `FEspacing`. Phase sign is flipped: `data = conj(data)`.

A center-of-mass origin is computed from the inf-VENC magnitude image. All subsequent coordinates (`FEgrid`, `PEgrid`, `rGrid`, `pGrid`) are relative to this origin. The angular convention is:

```matlab
pGrid = -atan2(FEgrid, PEgrid);  % theta=0 → +PE direction
```

so `r*cos(p) = PE`, `-r*sin(p) = FE`.

**Physical tube dimensions** (hard-coded):

| Symbol | Value | Meaning |
|---|---|---|
| `ID` | 6.35 mm | inner diameter (lumen boundary) |
| `OD` | 11.11 mm | outer diameter (tissue boundary) |
| `bestVenc` | 10 cm/s | venc for velocity map |

---

## Forward model

### Coordinate system and velocity profile

The velocity field uses polar coordinates `(r, p)` defined relative to the global origin. Internally `velocity_func_ellipse` converts these to Cartesian displacements from the velocity peak:

```matlab
dPE = r.*cos(p)  - PEoffset;
dFE = -r.*sin(p) - FEoffset;
```

The parabolic profile with elliptical cross-section is:

$$v_i = V_{\max} \cdot \max\!\left(0,\ 1 - \left(\frac{r_{v,i}}{R_{\text{eff},i}}\right)^2\right)$$

where $r_{v,i} = \sqrt{dPE_i^2 + dFE_i^2}$ is the radial distance from the velocity peak and the direction-dependent lumen radius is:

$$R_{\text{eff},i} = \frac{R}{\sqrt{A_i^2 + AR^2\, B_i^2}}$$

with $A_i$ and $B_i$ the unit displacement vector $(dPE/r_v, dFE/r_v)$ projected onto the ellipse axes:

$$A_i = \hat{u}_{PE}\cos\alpha + \hat{u}_{FE}\sin\alpha, \qquad B_i = -\hat{u}_{PE}\sin\alpha + \hat{u}_{FE}\cos\alpha$$

At $\alpha = 0$: $A_i = \hat{u}_{PE}$, $B_i = \hat{u}_{FE}$. In the PE direction $R_{\text{eff}} = R$ (semi-major), in the FE direction $R_{\text{eff}} = R/AR$ (semi-minor, $AR \ge 1$). The lumen is the set $\{v_i > 0\}$, equivalently $\{r_{v,i} < R_{\text{eff},i}\}$.

### Vessel angle parametrization

The vessel tilt is represented by a 2D Cartesian projection of the vessel unit vector onto the slice plane:

$$\vec{t} = (t_x, t_y), \qquad |\vec{t}| \le 1$$

Derived quantities (computed by `vessel_angle_params13`):

$$\cos\theta = \sqrt{1 - t_x^2 - t_y^2}, \qquad AR = \frac{1}{\cos\theta}, \qquad \alpha = \text{atan2}(t_y, t_x)$$

The effective slice thickness along the vessel axis is:

$$d_{\text{eff}} = \frac{d_{\text{slice}}}{\cos\theta} = d_{\text{slice}} \cdot AR$$

This single pair $(t_x, t_y)$ simultaneously encodes: (1) the ellipse aspect ratio and rotation of the projected vessel cross-section, and (2) the longer effective slab that tilted vessels traverse, which governs inflow saturation.

For the phantom (perpendicular vessel): $t_x = t_y = 0 \Rightarrow AR = 1$, $d_{\text{eff}} = d_{\text{slice}}$, circular lumen.

### Steady-state magnetization: `getMz_ss`

For a spin with velocity $v$ (cm/s), the number of RF pulses it experiences while crossing the slab is:

$$n(v) = \left\lceil \frac{v_c}{v} \right\rceil, \qquad v_c = 0.1 \cdot \frac{d_{\text{eff}}}{TR}$$

Three regimes:

| Regime | Condition | $M_z^{\text{ss}}(v)$ |
|---|---|---|
| Stationary | $v = 0$ | $\displaystyle M_0 \frac{1 - E_1}{1 - Q_1}$ |
| Partial inflow | $0 < v < v_c$ | $\displaystyle M_{z0} + (M_0 - M_{z0}) \cdot \frac{1 - Q_1^{n}}{n(1-Q_1)}$ |
| Full inflow | $v \ge v_c$ | $M_0$ |

where $E_1 = e^{-TR/T_1}$, $Q_1 = E_1\cos\alpha_{\text{FA}}$, $M_{z0} = M_0(1-E_1)/(1-Q_1)$ is the stationary steady state, and $n = \lceil v_c / v \rceil$.

The partial-inflow formula is a **staircase** in velocity: within each interval $v \in (v_c/(n+1),\, v_c/n]$, the integer $n = \lceil v_c/v \rceil$ is constant and $M_z^{\text{ss}}$ takes a constant value. Steps occur at $v = v_c/n$ for $n = 2, 3, \ldots$ This is the correct per-isochromat physics — the staircase should not be smoothed at this level. At $v \to 0^+$: $n \to \infty$, the factor $(1-Q_1^n)/(n(1-Q_1)) \to 0$, so $M_z^{\text{ss}} \to M_{z0}$. At $v$ just below $v_c$: $n = 2$, giving a value below $M_0$; regime 3 then jumps directly to $M_0$, leaving a small discontinuity at $v_c$. The apparent smoothing seen in voxel data arises from the distribution of spin velocities within each voxel — this is handled in the fitting by the 7×7 sub-voxel spin grid, not by smoothing the staircase itself.

**Implementation note:** the effective pMri struct with scaled `sliceThickness` is passed to `getMz_ss`, so $d_{\text{eff}}$ is automatically used inside the function without any special handling.

### Transverse magnetization: `getMxy_ss`

$$M_{xy}^{\text{ss}}(v) = M_z^{\text{ss}}(v) \cdot Q_2, \qquad Q_2 = e^{-TE/T_2^*} \cdot \sin\alpha_{\text{FA}}$$

This is the expected signal magnitude at echo time TE, before the amplitude scaling $A$.

### Inflow magnitude function: `inflowMag13`

$$m_{\text{inflow}}(v;\, A, \theta_{\text{MRI,eff}}) = A \cdot M_{xy}^{\text{ss}}\!\left(M_z^{\text{ss}}(v;\, \theta_{\text{MRI,eff}}),\, \theta_{\text{MRI,eff}}\right)$$

`inflowMag13` forces double precision to avoid type errors from `getMz_ss` returning single when `pRelax.T1` is single:

```matlab
function m = inflowMag13(v, A, pMri_eff, pRelax)
[Mz, pMri_u] = getMz_ss(pMri_eff, pRelax, double(v(:)'));
m = double(A) * double(getMxy_ss(double(Mz), pMri_u, pRelax));
m = reshape(m, size(v));
end
```

### Three-compartment signal model

The predicted voxel signal is computed by averaging over a 7×7 sub-voxel spin grid. Each sub-spin is assigned to a compartment based on its own position, and its Mxy is computed from the per-spin (staircase) inflow model. Averaging across the 49 sub-spins gives the predicted voxel signal — the same natural smoothing that `simVesselSpins` implements, without interpolating the staircase.

| Compartment (per sub-spin) | Condition | $\hat{M}_{xy}$ |
|---|---|---|
| Lumen | $v_{\text{spin}} > 0$ | $A \cdot M_{xy}^{\text{ss}}(v_{\text{spin}})$ — staircase |
| Wall | $v_{\text{spin}} = 0$ and $r_{v,\text{spin}} < R_{\text{eff},\text{spin}} + WT$ | $0$ (fixed) |
| Tissue | everything else | $S_t$ (free scalar) |

```matlab
% 7×7 sub-voxel spin grid
[dfe_sg, dpe_sg] = ndgrid(linspace(-0.5,0.5,7)*FEspacing, linspace(-0.5,0.5,7)*PEspacing);
v_sg     = velocity_func_ellipse(...);          % velocity per sub-spin
lumen_sg = v_sg > 0;
tissue_sg = ~lumen_sg & (rv_sg >= Reff_sg + WT);
mxy_sg(lumen_sg)  = inflowMag13(v_sg(lumen_sg), A, pMri_eff, pRelax);
mxy_sg(tissue_sg) = S_tissue;
Mxy = reshape(mean(mxy_sg, 2), size(FEgrid));   % voxel-averaged
```

Pixel-center velocities (`v_pred`, `lumen`) are still computed for the velocity residuals only — they are not used for the magnitude prediction.

The outer wall boundary per sub-spin uses a constant radial offset $WT$ from the inner ellipse $R_{\text{eff},\text{spin}}$, consistent with the original fit definition.

**Caveat — figure overlay vs. fit definition:** the outer wall ellipse drawn on the maps uses semi-axes $(R + WT,\; R/AR + WT)$, which is a concentric ellipse. This does **not** match the radial offset used in the fit residuals. For the phantom ($AR \approx 1$) both are nearly identical circles, but for in-vivo data with $AR > 1$ the displayed outer ellipse will not match the fitted wall boundary.

### Pre-realized noise model

Thermal noise is treated as complex Gaussian with standard deviation $\sigma_n$. To avoid re-drawing noise each optimizer iteration (which would make the gradient stochastic), noise is drawn once before optimization with `rng(0)`:

```matlab
rng(0);
noise_grid   = randn(size(M_ph)) + 1i * randn(size(M_ph));  % flow data
noise_noflow = randn(size(M_ph)) + 1i * randn(size(M_ph));  % noFlow (independent draw)
```

The predicted magnitude is then:

$$m_i^{\text{pred}} = \sqrt{\left(\hat{M}_{xy,i} + \sigma_n \operatorname{Re}(\eta_i)\right)^2 + \left(\sigma_n \operatorname{Im}(\eta_i)\right)^2}$$

where $\eta_i \sim \mathcal{CN}(0,1)$ are the pre-realized noise samples. This is the exact signal magnitude under complex Gaussian noise when Mxy is real-valued. $\sigma_n$ is a free parameter — it scales the fixed noise realization each iteration.

**Why wall signal is nonzero:** wall voxels have $\hat{M}_{xy} = 0$, so $m_i^{\text{pred}} = \sigma_n |\eta_i|$. The optimizer can fit wall signal by adjusting $\sigma_n$ alone. The expected wall signal is $\sigma_n \sqrt{\pi/2}$ (Rayleigh mean), which is why the initial estimate is:

```matlab
sigma_n_init = double(mean(M_ph(mask_wall_init))) * sqrt(2/pi);
```

### MRI acquisition parameters (phantom 03)

```matlab
pMri_ph.fieldStrength  = 3;            % T
pMri_ph.sliceThickness = 2.2;          % mm
pMri_ph.TR             = 75.90/(5+1)/1000;  % s  = 12.65 ms
pMri_ph.TE             = 9.8/1000;     % s
pMri_ph.FA             = 50;           % deg
```

Relaxation parameters (`pRelax_ph`) are fetched from `runSim` defaults for `species='phantom'`. Note: `runSim` is called with no arguments first to get default parameters, then again with `pMri_ph` substituted to populate the struct correctly.

---

## Cost function

Minimized by `lsqnonlin` (trust-region-reflective, `MaxFunctionEvaluations = 3×10^4`, `FunctionTolerance = 10^{-9}`).

The residual vector has three stacked components (implemented in `residuals_inflow13_full`):

$$\mathbf{r}(\boldsymbol{\theta}) = \begin{bmatrix} \mathbf{r}_{\text{mag}} \\ \mathbf{r}_{\text{vel}} \\ \mathbf{r}_{\text{noFlow}} \end{bmatrix}$$

**Magnitude residuals** (all $N_{\text{pix}}$ pixels):

$$r_{\text{mag},i} = \frac{m_i^{\text{meas}} - m_i^{\text{pred}}(\boldsymbol{\theta})}{\sigma_m}, \qquad \sigma_m = \text{std}(M^{\text{meas}}_{\text{all}})$$

**Velocity residuals** (blood-masked lumen pixels only):

$$r_{\text{vel},i} = \frac{v_i^{\text{meas}} - v_i^{\text{pred}}(\boldsymbol{\theta})}{\sigma_v}, \qquad \sigma_v = \text{std}(v^{\text{meas}}_{\text{blood}})$$

The velocity mask is `mask_vel & lumen`, where `mask_vel = maskBlood_ph` is a fixed threshold mask (`M > 0.3 * max(M)`) and `lumen` changes each iteration as the vessel parameters change. The data count in the fit report slightly overestimates this as `sum(maskBlood_ph)`.

**noFlow magnitude residuals** (all $N_{\text{pix}}$ pixels, fit 13a only):

$$r_{\text{noFlow},i} = \frac{m_{\text{noFlow},i}^{\text{meas}} - m_{\text{noFlow},i}^{\text{pred}}(\boldsymbol{\theta})}{\sigma_m}$$

For the noFlow prediction, lumen pixels use $v=0$ (fully saturated signal $A \cdot M_{xy}^{\text{ss}}(0)$), wall pixels use $\sigma_n |\eta_i^{\text{noFlow}}|$, and tissue pixels use $S_t$. The noise realization `noise_noflow` is drawn independently from `noise_grid`.

---

## Free parameters

Parameter vector: `theta = [Vmax, R, tx, ty, A, FEoffset, PEoffset, WT, S_tissue, sigma_n]`

| # | Name | Units | Lower bound | Upper bound | Initial value (fit 13a) | Source |
|---|---|---|---|---|---|---|
| 1 | $V_{\max}$ | cm/s | 0 | $\infty$ | `max(abs(v_blood))` | peak measured velocity in blood mask |
| 2 | $R$ | mm | $10^{-6}$ | `ID` = 6.35 | `ID/2` = 3.175 | inner lumen radius (ID/2) |
| 3 | $t_x$ | — | −0.7 | 0.7 | 0 | perpendicular vessel |
| 4 | $t_y$ | — | −0.7 | 0.7 | 0 | perpendicular vessel |
| 5 | $A$ | a.u. | 0 | $\infty$ | `mean(m_noflow_blood) / Mxy_ss(v=0)` | noFlow signal divided by fully-saturated model |
| 6 | FEoffset | mm | $-$`FEspacing` | `FEspacing` | 0 | center of FOV |
| 7 | PEoffset | mm | $-$`PEspacing` | `PEspacing` | 0 | center of FOV |
| 8 | $WT$ | mm | 0 | `OD` = 11.11 | `OD/2 − ID/2` = 2.38 | physical tube wall thickness |
| 9 | $S_t$ | a.u. | 0 | $\infty$ | `mean(M_ph(rGrid > OD/2))` | mean signal outside outer tube |
| 10 | $\sigma_n$ | a.u. | 0 | $\infty$ | `mean(M_ph(mask_wall)) * sqrt(2/pi)` | Rayleigh inversion of mean wall signal |

**Bounds note:** $|t_x|, |t_y| \le 0.7$ enforces $\cos\theta \ge \sqrt{1-2\times0.49} = \sqrt{0.02} \approx 0.14$, i.e., vessel tilt $\theta \le 82°$. The $R$ upper bound is `ID` (full lumen radius, not `ID/2`) which allows the fit to explore a factor-of-2 range.

Fit 13b initial values are the fit 13a final values. Same bounds.

**Derived parameters** (reported but not free):

$$\theta_{\text{vessel}} = \arccos\!\left(\sqrt{1-t_x^2-t_y^2}\right), \quad AR = 1/\cos\theta, \quad \alpha = \text{atan2}(t_y, t_x)$$

---

## Fit sequence

```
theta0_13  →  lsqnonlin(f_res_a)  →  theta_13a
                                           |
                                    lsqnonlin(f_res_b)  →  theta_13b
```

- **Fit 13a** (`f_res_a`): residuals include magnitude + velocity + noFlow. The noFlow term anchors $A$ and $\sigma_n$ independently of the inflow shape.
- **Fit 13b** (`f_res_b`): residuals include magnitude + velocity only. Warm-started from `theta_13a`. The inflow shape must now be inferred purely from the signal variation across the lumen.

Both fits use the same `residuals_inflow13_full` function; fit 13b simply passes `m_noflow = []` which causes `res_nf = []` to be omitted.

The figure and parameter reports all use fit 13b results.

---

## Figure panels (A–H)

The figure is a 2×4 tiled layout (`tiledlayout(f_13, 2, 4, ...)`). All panels have dark (`'Color','k'`) axes backgrounds. Markers: circle (`'o'`), size 4, white face, black edge, line width 0.5.

**Compartment colors** (from 7×7 sub-grid partial volumes):

| Color | Compartment |
|---|---|
| Cyan `[0.15 0.85 1.00]` | lumen (blood) |
| Orange `[1.00 0.45 0.05]` | wall |
| Green `[0.45 0.90 0.45]` | tissue |

The sub-grid assigns each voxel to the dominant compartment by evaluating `velocity_func_ellipse` on a 7×7 grid within that voxel (vectorized as 625×49 matrices) and taking the fraction of sub-points with $v>0$ (lumen), $v=0$ and $r_v < R_{\text{eff}} + WT$ (wall), or the remainder (tissue).

---

### Panel A — magnitude map (tile 1)

`imagesc(PEpos, FEpos, M_ph)` of the inf-VENC averaged magnitude. Overlaid:
- Red solid ellipse: inner lumen boundary ($R_b$, fit 13b)
- Red dashed ellipse: outer wall boundary ($R_b + WT_b$, fit 13b)
- Red `+`: vessel center (`theta_13b(7), theta_13b(6)` = PEoffset, FEoffset)

Ellipse generation:
```matlab
PE_in  = R_b .* cos(t_c);   FE_in  = (R_b/AR_b) .* sin(t_c);
cx_in  = PE_in.*cos(alpha_b_rad) - FE_in.*sin(alpha_b_rad) + theta_13b(7);
cy_in  = PE_in.*sin(alpha_b_rad) + FE_in.*cos(alpha_b_rad) + theta_13b(6);
```
The PE component is plotted as x and FE as y (matching `imagesc` axes). See caveat above about outer ellipse shape.

---

### Panel B — magnitude radial profile (tile 2)

Scatter of `M_ph` vs. `rGridOff_ph` (radial distance from vessel center = `sqrt((FE-FEoffset)² + (PE-PEoffset)²)`). All voxels plotted, colored by dominant compartment.

Three 1D fit curves along $p=0$ (PE direction):
- **Lumen** ($r < R_b$): `inflowMag13(v1D_b(r), theta_13b(5), pMri_eff_b, pRelax_ph)` where `v1D_b(r) = velocity_func_ellipse(r, 0, ...)` (zero angle = PE direction, no offset)
- **Wall** ($R_b \le r < R_b + WT_b$): $\sigma_n \sqrt{\pi/2}$ (Rayleigh expected value of noise magnitude, constant)
- **Tissue** ($r \ge R_b + WT_b$): $S_t$ (constant)

```matlab
m_r_wall_b = theta_13b(10) * sqrt(pi/2) * ones(size(r_wall_plt));
```

Vertical lines: `xline(R_b, 'r-')` and `xline(R_b + WT_b, 'r--')`.

---

### Panel C — predicted vs. measured magnitude scatter (tile 3)

`m_pred_full` vs. `M_ph` for all pixels. `m_pred_full` is computed from the full forward model:

```matlab
v_pred_full = velocity_func_ellipse(rGrid, pGrid, ...);
Mxy_f(lumen_f) = inflowMag13(...);
Mxy_f(tissue)  = theta_13b(9);
% wall: Mxy_f = 0
m_pred_full = sqrt((Mxy_f + theta_13b(10)*real(noise_grid)).^2 + ...
                   (theta_13b(10)*imag(noise_grid)).^2);
```

The identity diagonal `[0,max]→[0,max]` is plotted first (dark gray, `LineWidth=0.8`) then markers on top. Both axes share the same range `xy_lim_m`.

---

### Panel D — inflow model $m(v)$ (tile 4)

Scatter of measured blood-mask magnitudes vs. measured velocities, colored by dominant compartment within the blood mask:
```matlab
cm_b = comp_map(maskBlood_ph);
```

Overlaid: fit 13b inflow curve `inflowMag13(v_plt, theta_13b(5), pMri_eff_b, pRelax_ph)` for $v \in [0,\; 1.1 V_{\max}]$. Horizontal yellow dotted line: `mean(m_noflow_blood)` (noFlow reference level).

---

### Panel E — velocity map (tile 5)

`imagesc(PEpos, FEpos, vFlow_ph, [-bestVenc, bestVenc])` using the `blueBlackRed` colormap. Same ellipse and center overlays as panel A.

`vFlow_ph = phase2vel(angle(cBest_ph), vencToM1(bestVenc))` converts phase to cm/s. The phase-to-velocity relationship is:

$$v = \frac{\phi \cdot 100}{\gamma_{\text{Hz}} \cdot M_1} \qquad \text{(cm/s)}, \quad M_1 = \frac{\pi}{\gamma_{\text{Hz}} \cdot v_{\text{enc}}/100}$$

---

### Panel F — velocity radial profile (tile 6)

Scatter of `vFlow_ph` vs. `rGridOff_ph`, all voxels, compartment-colored. Fit curves:
- **Lumen:** `v1D_b(r) = velocity_func_ellipse(r, 0, Vmax_b, R_b, AR_b, alpha_b_rad, 0, 0)` — parabola evaluated in PE direction
- **Wall and tissue:** zero velocity

Vertical boundary lines same as panel B.

---

### Panel G — complex domain (tile 7)

Data trajectory: `plotComplexDomain(ax_cd, trj_ph_n, all_vencs_ph, 'full', 'markers')`.

`trj_ph_n` is the normalized average complex blood signal at each venc, sorted descending by venc (inf first):

```matlab
trj_ph(1) = mean(cFlow_ph(maskBlood_ph));          % inf-venc reference
for kk = 1:numel(finiteVencs_ph)
    cVenc_kk = squeeze(mean(data(:,:,dataVenc==finiteVencs_ph(kk)),3));
    trj_ph(1+kk) = mean(cVenc_kk(maskBlood_ph));
end
trj_ph_n = trj_ph / abs(trj_ph(1));
```

The predicted trajectory is computed **analytically** from blood pixels only (to avoid dilution from static tissue), then phase-aligned to the data:

```matlab
gamma_hz  = 2.6752218708e8 / (2*pi);          % Hz/T
v_b_ms    = double(v_pred_full(maskBlood_ph)) / 100;   % cm/s → m/s
Mxy_b_cd  = double(Mxy_f(maskBlood_ph));
m1_traj   = linspace(0, vencToM1(min(finiteVencs_ph)), 600)';
phases_cd = gamma_hz .* m1_traj .* v_b_ms(:)';   % 600 × N_blood
I_pred_cd = mean(Mxy_b_cd(:)' .* exp(1j .* phases_cd), 2);  % 600×1
I_pred_cd_n = I_pred_cd / abs(I_pred_cd(1));
I_pred_cd_n = I_pred_cd_n * exp(1j * angle(trj_ph_n(1)));
```

The phase encoding formula is $\phi = \gamma_{\text{Hz}} \cdot M_1 \cdot v\;[\text{m/s}]$. The sweep covers $M_1 = 0$ (inf-venc, no encoding) to $M_1^{\text{max}}$ (most aggressive finite venc). The phase alignment ensures that the predicted trajectory starts at the same complex angle as the data at $M_1=0$.

The predicted curve (`h_sl`, cyan line) is sent to the bottom of the stack with `uistack(h_sl, 'bottom')`.

---

### Panel H — predicted vs. measured phase scatter (tile 8)

```matlab
phase_pred_full = double(pi * v_pred_full / bestVenc);
phase_meas_full = double(angle(cBest_ph));
```

The predicted phase uses $\phi = \pi v / v_{\text{enc}}$, derived from $\phi = \gamma_{\text{Hz}} M_1 v[\text{m/s}]$ with $M_1 = \pi / (\gamma_{\text{Hz}} \cdot v_{\text{enc}}[\text{m/s}])$.

`cBest_ph = squeeze(mean(data(:,:, dataVenc==bestVenc), 3))` is the average complex signal at the best VENC. The phase of this is the background-corrected velocity phase (the `conj` applied at load time has already adjusted the sign convention).

Both axes fixed to $[-\pi, \pi]$ with $\pi/2$ ticks.

---

### Panel labels A–H

Applied after all panels are populated:

```matlab
set(findall(f_13,'Type','axes'),'FontSize',12);
set(findall(f_13,'Type','text'),'FontSize',8);
for k_ = 1:8
    text(nexttile(tl_13,k_), 0.02, 0.97, char('A'+k_-1), ...
        'Units','normalized', 'FontSize',13, 'FontWeight','bold', ...
        'Color','w', 'VerticalAlignment','top', 'HorizontalAlignment','left');
end
```

The `set(... 'FontSize', 8)` call **resets all text objects (including labels) to size 8 first**, then the panel-label loop writes size-13 labels. Order matters.

---

### Figure export

Two exports: PNG (standard) and SVG. For the SVG, circular markers are swapped to dots to produce vector-friendly output:

```matlab
hMkr_ = findobj(f_13,'Marker','o');
origMEC_ = get(hMkr_,{'MarkerEdgeColor'});
arrayfun(@(h) set(h,'MarkerEdgeColor',h.MarkerFaceColor), hMkr_);
set(hMkr_,'Marker','.');
exportgraphics(f_13, fullfile(sec13fig,'phantom_inflow.svg'));
set(hMkr_,'Marker','o');
set(hMkr_,{'MarkerEdgeColor'}, origMEC_);
```

---

## Known limitations and potential issues

1. **Outer wall overlay inconsistency for $AR \ne 1$:** the plotted outer wall ellipse uses semi-axes $(R+WT,\; R/AR+WT)$, which is not a constant radial offset from the inner ellipse. The fit uses `r_v < R_eff_in + WT` (constant radial offset). For the phantom ($AR \approx 1$) both are circles. For tilted in-vivo vessels this will create a visible mismatch between the overlay and the actual fit boundary.

2. **Velocity residual count in fit report:** `dsum_13b` counts `sum(maskBlood_ph)` velocity residuals but the actual count is `sum(maskBlood_ph & lumen_final)`, which can be smaller if the optimizer shrinks the lumen. Minor effect.

3. **Discontinuity in `getMz_ss` at $v = v_c$:** regime 2 with $v \to v_c^-$ gives $n=2$, which does not exactly equal regime 3 ($M_0$). The discrepancy is small but produces a step at $v_c$ in the per-isochromat staircase. This is physically correct — spins crossing the slab in exactly one TR have $n=1$ but that case only occurs at $v \ge v_c$. The sub-voxel averaging in the fitting further softens this edge.

4. **Phase scatter (panel H) background phase:** `angle(cBest_ph)` contains any residual background phase not removed by the reference scan. For the phantom this is negligible, but for in vivo data with B0 inhomogeneity this offset would shift the measured phase cluster relative to the prediction.

5. **`sigma_n_init` can be NaN:** if `mask_wall_init = rGrid > R_ph & rGrid <= OD/2` yields no pixels (e.g., very small FOV), `mean(M_ph(mask_wall_init))` is NaN. There is no guard (unlike `S_tissue_init`). Not an issue for the current phantom geometry.

6. **Legacy `residuals_inflow13` stub** at end of file is unused (superseded by `residuals_inflow13_full`). Can be removed once the branch is stable.

7. **`writeFitParamsMd` docstring mismatch:** the function signature says `dataSummary` should be a struct, but the calling code passes a plain string (e.g. `dsum_13b`). The function body uses `fprintf(..., '%s\n\n', dataSummary)`, so it works correctly with a string. The docstring is misleading.

---

## File and branch map

| File | Role |
|---|---|
| `doIt_inflow.m` | entry point for this branch |
| `doIt_inflow.md` | this file |
| `velocity_func_ellipse.m` | elliptical parabolic flow profile |
| `writeFitParamsMd.m` | writes `phantom_fitParams_13a/b.md` |
| `plotComplexDomain.m` | complex-plane plot (panel G) |
| `pcMRAsim/getMz_ss.m` | steady-state Mz (inflow model) |
| `pcMRAsim/getMxy_ss.m` | transverse magnetization |
| `util/vencToM1.m` | venc [cm/s] → M1 [T·s²/m] |
| `util/phase2vel.m` | phase [rad] → velocity [cm/s] |
| `figures/13-inflow-model/phantom_inflow.png` | output figure |
| `figures/13-inflow-model/phantom_fitParams_13a.md` | fit 13a parameter report |
| `figures/13-inflow-model/phantom_fitParams_13b.md` | fit 13b parameter report |

All repos on `dev-inflow` branch: `multiVencISMRM2026`, `util`, `pcMRAsim`, `multiVencSim`, `db/*`, `data`, `sim`.

---

## Related files

- `complexFitMath.md` — Fit B residual derivation and parameter evolution
- `asymmetricProfile.md` — ellipse wall geometry background
- `doIt.m` sec 11–12 — current Fit A / Fit B implementations
