# Fit 13a — Inflow 3-compartment (phantom, with noFlow)

## Parameters

`*` marks parameters fixed at their listed value (not optimised).

```
Parameter         Units       Fixed?             LB            UB       Initial         Final
---------------------------------------------------------------------------------------------
Vmax              cm/s        no                  0           Inf         9.027         9.295
R                 mm          no              1e-06          6.35         3.175         3.208
nx                            no               -0.1           0.1             0             0
ny                            no               -0.1           0.1             0             0
cx_FE             mm          no               -0.5           0.5             0             0
cx_PE             mm          no            -0.8929        0.8929             0             0
cx_SLC *          mm          yes                 0             0             0             0
A                 a.u.        no                  0           Inf     4.663e-06     1.655e-06
WT                mm          no                  0         11.11          2.38          2.38
S_tissue          a.u.        no                  0           Inf     5.888e-08     5.924e-08
sigma_n           a.u.        no                  0           Inf     1.127e-08     4.946e-09
```

Derived:

```
Parameter         Units            Initial         Final
--------------------------------------------------------
theta_vessel_deg  deg                    0             0
AR                                       1             1
alpha_deg         deg                    0             0
```

### Initial value sources

| Parameter | Source |
|---|---|
| Vmax | max|v_blood| from best-VENC image |
| R | ID/2 (inner lumen radius) |
| nx | 0 (perpendicular vessel) |
| ny | 0 (perpendicular vessel) |
| cx_FE | 0 |
| cx_PE | 0 |
| cx_SLC * | 0 |
| A | mean(noFlow_blood)/getMxy_ss(v=0) |
| WT | OD/2-ID/2 (physical tube wall) |
| S_tissue | mean(pixels outside OD) |
| sigma_n | mean(wall pixels)*sqrt(2/pi) |

---

## Model

### Velocity profile

Parabolic profile with elliptical wall; peak at $(\text{FEoffset}, \text{PEoffset}) = (0, 0)$:

$$v_i = V_{\max} \cdot \max\!\left(0,\ 1 - \left(\frac{r_{v,i}}{R_{\text{eff},i}}\right)^2\right)$$

where $r_{v,i}$ is the distance from the velocity peak, and the direction-dependent wall radius is

$$R_{\text{eff},i} = \frac{R}{\sqrt{A_i^2 + AR^2\, B_i^2}}$$

with $(A_i, B_i)$ the unit direction projected onto the ellipse axes (semi-major $R$ at angle $\alpha$ from PE axis, semi-minor $R/AR$, $AR \ge 1$).

### Magnitude–velocity relationship

$$m_{\mathrm{phys}}(v) = |\hat{M}_{xy}(v) + \sigma_n \eta|$$\n\nLumen: $\hat{M}_{xy} = A \cdot M_{xy}^{\mathrm{ss}}(\tau(v))$. Wall: $\hat{M}_{xy}=0$ (fixed). Tissue: $\hat{M}_{xy}=S_t$.

---

## Cost function

Minimised by `lsqnonlin` (trust-region reflective).

$$\mathbf{r} = \begin{bmatrix}(m_i^{\mathrm{meas}} - |\hat{M}_{xy,i} + \sigma_n \eta_i|)/\sigma_m \\(v_i^{\mathrm{meas}} - v_i^{\mathrm{pred}})/\sigma_v\end{bmatrix}$$

where the predicted signal is

$$|\hat{M}_{xy,i} + \sigma_n \eta_i|,\quad \eta_i \sim \mathcal{CN}(0,1),\text{ pre-realized}$$

$\sigma_m = \mathrm{std}(M^{\mathrm{meas}}_{\mathrm{all}})$, $\sigma_v = \mathrm{std}(v^{\mathrm{meas}}_{\mathrm{blood}})$. Three compartments: lumen (inflow model), wall ($\hat{M}_{xy}=0$, WT free), tissue ($S_t$ free). Noise $\eta_i$ pre-realized once; $\sigma_n$ scales it each iteration.

**Data**: 486 px (all) + 59 vel + 486 noFlow. Total: 1031.

