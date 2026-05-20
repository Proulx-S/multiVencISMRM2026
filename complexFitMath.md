# Vessel profile fitting — model evolution and planned change

---

## Stage 1 — Branch point (`multiVencISMRM2026`, phantom)

**Velocity**: symmetric circle, peak fixed at ROI centre

$$v(r) = V_{\max}\left(1 - \frac{r^2}{R^2}\right)$$

**Magnitude**: polynomial with $B$ fixed from a no-flow reference scan

$$m(v) = B_{\text{fixed}} + C_1 v + C_2 v^2$$

**Residual** (single VENC, real-valued):

$$\mathbf{r} = \begin{bmatrix}(v^{\text{meas}} - v^{\text{pred}})/\sigma_v \\ (m^{\text{meas}} - m(v^{\text{pred}}))/\sigma_m\end{bmatrix}$$

**Free parameters** (4): $V_{\max},\ R,\ C_1,\ C_2$

---

## Stage 2 — In vivo additions (`dev-postISMRM`)

Three successive generalisations, all preserving the single-VENC real residual.

### 2a. Free $B$

No no-flow reference scan exists in vivo, so $B$ becomes a free parameter initialised from blood $T_1$ physics:

$$B_0 = \frac{M_z^{\text{ss}}(0)}{M_{xy}^{\text{ss}}(V_{\max}/2)} \quad \text{(ratio from \texttt{getMz\_ss} / \texttt{getMxy\_ss})}$$

**Added parameter**: $B$ (now free) → **5 parameters**

### 2b. Centre offset

The velocity peak is no longer forced to the ROI centre; $(\text{FEoffset},\ \text{PEoffset})$ are fitted freely. The effective radius from the peak becomes

$$r_{v,i} = \sqrt{(\delta\text{PE}_i)^2 + (\delta\text{FE}_i)^2}$$

where $\delta\text{PE}_i = r_i\cos p_i - \text{PEoffset}$ and $\delta\text{FE}_i = -r_i\sin p_i - \text{FEoffset}$.

**Added parameters**: $\text{FEoffset},\ \text{PEoffset}$ → **7 parameters**

### 2c. Elliptical wall (`offset='ellipse'`, current)

The circular wall is replaced by an ellipse with semi-major axis $R$, aspect ratio $AR \ge 1$, and orientation $\alpha$. The velocity peak coincides with the ellipse centre. The direction-dependent wall radius is

$$R_{\text{eff},i} = \frac{R}{\sqrt{A_i^2 + AR^2\, B_i^2}}$$

where $A_i,\ B_i$ are the components of the unit direction $\hat{u}_i$ (from peak to pixel) projected onto the ellipse axes:

$$A_i = \hat{u}_{i,\text{PE}}\cos\alpha + \hat{u}_{i,\text{FE}}\sin\alpha, \qquad B_i = -\hat{u}_{i,\text{PE}}\sin\alpha + \hat{u}_{i,\text{FE}}\cos\alpha$$

The velocity profile becomes

$$v_i = V_{\max}\cdot\max\!\left(0,\ 1 - \left(\frac{r_{v,i}}{R_{\text{eff},i}}\right)^2\right)$$

**Added parameters**: $AR,\ \alpha$ → **9 parameters**

---

## Current model summary

| Symbol | Role | Free? |
|--------|------|-------|
| $V_{\max}$ | Peak velocity | yes |
| $R$ | Semi-major wall radius | yes |
| $AR$ | Wall aspect ratio | yes |
| $\alpha$ | Wall orientation | yes |
| $\text{FEoffset}$ | Peak position (FE) | yes |
| $\text{PEoffset}$ | Peak position (PE) | yes |
| $B$ | Zero-velocity magnitude | yes |
| $C_1$ | Linear inflow coefficient | yes |
| $C_2$ | Quadratic inflow coefficient | yes |

**Residual**: single-VENC, real — $2N$ constraints on 9 parameters.

---

## Planned change — complex multi-VENC residual

Keep all 9 parameters and the polynomial magnitude model unchanged. The only change is to the residual.

**Predicted complex signal** for pixel $i$ at VENC $k$:

$$\hat{s}_{ik} = m(v_i)\cdot e^{\,i\pi v_i/\mathrm{venc}_k}, \qquad m(v) = B + C_1 v + C_2 v^2$$

No phase unwrapping needed — wrapping is handled naturally by the complex exponential.

**New residual** (all finite VENCs simultaneously):

$$\mathbf{r} = \frac{1}{\sigma_s}\begin{bmatrix}\operatorname{Re}\!\left(\mathbf{s}^{\text{meas}} - \hat{\mathbf{s}}\right)\\\operatorname{Im}\!\left(\mathbf{s}^{\text{meas}} - \hat{\mathbf{s}}\right)\end{bmatrix} \in \mathbb{R}^{2NK}$$

where $\sigma_s = \mathrm{mean}(|s_{ik}^{\text{meas}}|)$.

**Constraint count** (with $N\approx100$ pixels, $K=9$ VENCs):

$$2NK \approx 1800 \quad\text{vs.}\quad 2N \approx 200 \text{ (current)}$$

Same 9 parameters, roughly $9\times$ more constraints.
