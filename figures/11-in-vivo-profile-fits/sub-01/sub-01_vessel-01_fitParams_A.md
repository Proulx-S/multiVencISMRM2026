# Fit A — Joint real residual (single-VENC)

## Parameters

`*` marks parameters fixed at their listed value (not optimised).

```
Parameter         Units       Fixed?             LB            UB       Initial         Final
---------------------------------------------------------------------------------------------
Vmax              cm/s        no                  0           Inf         8.848         9.049
R                 mm          no              1e-06         4.444         2.222         2.325
e1                -           no                 -1             1             0       -0.2693
e2                -           no                 -1             1             0       -0.0311
B                 a.u.        no                  0           Inf     2.516e-08     2.408e-08
C1                a.u.        no                  0           Inf             0     7.503e-08
C2                a.u.        no               -Inf             0             0    -1.192e-09
FEoffset *        mm          yes                 0             0             0             0
PEoffset *        mm          yes                 0             0             0             0
```

Derived:

```
Parameter         Units            Initial         Final
--------------------------------------------------------
AR                -                      1         1.271
alpha             rad                    0        -1.513
```

### Initial value sources

| Parameter | Source |
|---|---|
| Vmax | fitVelProfile: 1D parabolic fit to velocity map at best VENC |
| R | fitVelProfile: same 1D fit |
| e1 | 0 (circular: AR=1 → e1=e2=0) |
| e2 | 0 (circular: AR=1 → e1=e2=0) |
| B | Physics: Mxy(0)/Mxy(10 cm/s), human blood T1=1.66s, TR=8ms, FA=25deg |
| C1 | 0 (flat inflow curve) |
| C2 | 0 (flat inflow curve) |
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

$$\mathbf{r} = \begin{bmatrix}(v_i^{\mathrm{meas}} - v_i^{\mathrm{pred}})/\sigma_v \\ (m_i^{\mathrm{meas}} - m(v_i^{\mathrm{pred}}))/\sigma_m\end{bmatrix}$$

$\sigma_v = \mathrm{std}(v^{\mathrm{meas}})$, $\sigma_m = \mathrm{std}(m^{\mathrm{meas}})$. Only blood-masked pixels enter the cost (same single-spin-per-pixel approximation as Fit B).

**Data**: 20 blood pixels (M > 0.3·max), 1 VENC (venc = 10 cm/s). Total residual elements: 40.

