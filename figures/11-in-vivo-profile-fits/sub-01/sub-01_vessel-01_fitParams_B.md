# Fit B — runSim forward model (all voxels)

## Parameters

`*` marks parameters fixed at their listed value (not optimised).

```
Parameter         Units       Fixed?             LB            UB       Initial         Final
---------------------------------------------------------------------------------------------
Vmax              cm/s        no                  0           Inf         9.049         7.248
R                 mm          no              1e-06         2.667         2.325         2.646
e1                -           no                 -1             1       -0.2693       -0.5732
e2                -           no                 -1             1       -0.0311      0.003403
A_n               -           no                  0           Inf         4.033         5.712
FEoffset *        mm          yes                 0             0             0             0
PEoffset *        mm          yes                 0             0             0             0
```

Derived:

```
Parameter         Units            Initial         Final
--------------------------------------------------------
AR                -                  1.271         1.573
alpha             rad               -1.513         1.568
```

### Initial value sources

| Parameter | Source |
|---|---|
| Vmax | Fit A final value |
| R | Fit A final value |
| e1 | Fit A final e1 |
| e2 | Fit A final e2 |
| A_n | 1/Mxy_blood(Vmax/2) ≈ 4.033 (so A_n*Mxy ≈ 1 at init) |
| FEoffset * | Fixed at 0 |
| PEoffset * | Fixed at 0 |

---

## Model

### Velocity profile

Parabolic profile with elliptical wall; peak at $(\text{FEoffset}, \text{PEoffset}) = (0, 0)$:

$$v_i = V_{\max} \cdot \max\!\left(0,\ 1 - \left(\frac{r_{v,i}}{R_{\text{eff},i}}\right)^2\right)$$

where $r_{v,i}$ is the distance from the velocity peak, and the direction-dependent wall radius is

$$R_{\text{eff},i} = \frac{R}{\sqrt{A_i^2 + AR^2\, B_i^2}}$$

with $(A_i, B_i)$ the unit direction projected onto the ellipse axes (semi-major $R$ at angle $\alpha$ from PE axis, semi-minor $R/AR$, $AR \ge 1$).

### Magnitude–velocity relationship

$$m(v) = B + C_1 v + C_2 v^2$$

---

## Cost function

Minimised by `lsqnonlin` (trust-region reflective).

$$\mathbf{r} = \frac{1}{\sigma_s}\begin{bmatrix}\operatorname{Re}(s_{ik}^{\mathrm{meas}} - \hat{s}_{ik}) \\\operatorname{Im}(s_{ik}^{\mathrm{meas}} - \hat{s}_{ik})\end{bmatrix}$$

where the predicted signal is

$$\hat{s}_{ik}(\mathrm{voxel}\,i, \mathrm{VENC}\,k) = A_n \sum_{j \in \mathrm{voxel}\,i} M_{xy}(v_j) \cdot e^{\,i\pi v_j/\mathrm{venc}_k}$$

$M_{xy}(v)$ from blood T1 physics (`getMxy\_ss`). Partial volume and within-voxel phase dispersion are handled by the spin summation in `runSim` (`allVoxels` gridMode). All $66$ ROI voxels enter the cost — no masking. $\sigma_s = \mathrm{std}(|s^{\mathrm{meas}}|)$ over all voxels and VENCs. Spins used: 264 (≈ 4.0 per voxel).

**Data**: 66 ROI voxels (all), K = 9 finite VENCs (4  5  6  7  8 10 13 20 40 cm/s), 264 spins. Total residual elements: 1188.

