# Asymmetric radial velocity profile — development notes

**Context:** extension of `fitMagVelProfile.m` joint fit to allow an asymmetric parabolic profile
where the velocity peak centre and the vessel wall centre are distinct points.

---

## Motivation

The current model assumes a symmetric parabolic profile:

```
v(r_v) = Vmax * (1 - (r_v / R)^2)
```

where `r_v` is distance from the velocity-peak centre `(FEoffset, PEoffset)` and the v=0
boundary is the circle of radius R centred at the same point.

In practice, several effects can break this symmetry:
- **Oblique slice:** if the vessel axis is not perpendicular to the imaging plane, the
  cross-section is circular but the parabolic peak shifts toward one wall.
- **Dean flow / vessel curvature:** centrifugal effects shift the peak toward the outer wall.
- **Partial-volume asymmetry:** asymmetric wall overlap with voxels biases the apparent peak.

In all these cases the wall remains a circle, but its centre is displaced relative to the
velocity peak. The proposed model separates the two.

---

## Geometry

| Symbol | Meaning |
|--------|---------|
| `(FEoffset, PEoffset)` | velocity peak centre (current parameters) |
| `(eFE, ePE)` | wall-centre offset **relative to** the velocity peak centre |
| `R` | wall-circle radius |
| `r_v` | distance of a pixel from the velocity peak centre |
| `theta_v` | angle of a pixel from the velocity peak centre (`atan2(PE-PEoffset, FE-FEoffset)`) |

Wall circle centre: `(FEoffset + eFE, PEoffset + ePE)`.

---

## Derivation of R_eff(theta_v)

Cast a ray from the velocity peak in direction `theta_v`:

```
(FE(t), PE(t)) = (FEoffset + t·cos(theta_v),  PEoffset + t·sin(theta_v))
```

Squared distance from the wall centre to a point at parameter t:

```
d²(t) = (t·cos(theta_v) - eFE)² + (t·sin(theta_v) - ePE)²
       = t² - 2t·p + |e|²
```

where:
```
p  = eFE·cos(theta_v) + ePE·sin(theta_v)   % projection of wall offset onto ray direction
|e|² = eFE² + ePE²                          % squared wall-offset magnitude
```

Setting `d²(t) = R²` and solving for t (taking the positive root, i.e. forward intersection):

```
R_eff(theta_v) = p + sqrt(R² - |e|² + p²)
```

**Validity:** requires `R² - |e|² + p² ≥ 0`, i.e. `|e| ≤ R` (velocity peak must lie inside
the lumen). For `|e| << R` this simplifies to:

```
R_eff ≈ R + p  =  R + eFE·cos(theta_v) + ePE·sin(theta_v)
```

a first-order sinusoidal modulation of the effective radius in the direction of the wall offset.

---

## Full asymmetric velocity model

```
r_v      = sqrt((FE - FEoffset)² + (PE - PEoffset)²)
theta_v  = atan2(PE - PEoffset, FE - FEoffset)
p        = eFE·cos(theta_v) + ePE·sin(theta_v)
R_eff    = p + sqrt(R² - eFE² - ePE² + p²)
v        = Vmax · max(0, 1 - (r_v / R_eff)²)
```

Free parameters: **Vmax, R, FEoffset, PEoffset, eFE, ePE** (6 total).

---

## Implementation strategies

### Strategy A — sfit with (r, theta) inputs (preferred, consistent with current interface)

The fittype already uses polar inputs `(r, theta)` measured from the velocity peak.  
Extend the fittype body to compute `R_eff(theta)` internally:

```matlab
ft_vel_asym = fittype( ...
    @(Vmax, R, FEoffset, PEoffset, eFE, ePE, r, theta) ...
        velocity_func_radius_asym(r, theta, Vmax, R, eFE, ePE), ...
    'independent', {'r', 'theta'}, ...
    'coefficients', {'Vmax', 'R', 'FEoffset', 'PEoffset', 'eFE', 'PEoffset'});
```

where `velocity_func_radius_asym(r, theta, Vmax, R, eFE, ePE)` encapsulates the
`R_eff` computation above.

**Note:** `FEoffset`/`PEoffset` do NOT appear in the function body (they are encoded in
the `r` and `theta` inputs, which the caller computes from the raw grid). They are stored
as inert sfit coefficients so they remain accessible via `velFit.FEoffset` / `velFit.PEoffset`.

Calling convention remains `velFit(r, theta)` with `r = sqrt((FE-FEoffset)² + (PE-PEoffset)²)`
and `theta = atan2(PE-PEoffset, FE-FEoffset)` — same as the current offset fit.

**lsqnonlin residuals:** `residuals_joint_asym` takes the full FEPE grid and computes
r_v and theta_v internally, then calls the model. This avoids passing FEoffset/PEoffset
twice.

Bounds:
```
lb_eFE = -R0/4,  ub_eFE = R0/4
lb_ePE = -R0/4,  ub_ePE = R0/4
```

The constraint `|e| ≤ R` (velocity peak inside lumen) is automatically satisfied
given these bounds and `R ≥ R0/2`.

### Strategy B — separate `velocity_func_radius_asym.m` helper (simpler to isolate and test)

Write a standalone function:

```matlab
function v = velocity_func_radius_asym(r, theta, Vmax, R, eFE, ePE)
    p     = eFE.*cos(theta) + ePE.*sin(theta);
    R_eff = p + sqrt(max(0, R^2 - eFE^2 - ePE^2 + p.^2));
    v     = max(0, Vmax .* (1 - (r ./ R_eff).^2));
end
```

This is the function body for any of the fittypes above and also makes unit-testing straightforward.

### Strategy C — linear-perturbation fittype (quick approximation, fewer parameters)

For `|e| << R`, use:

```
v ≈ Vmax · (1 - r_v² / (R + eFE·cos(theta_v) + ePE·sin(theta_v))²)
```

Cheaper to evaluate and differentiate; adequate when the asymmetry is small.

---

## Open questions

1. **Identifiability:** with the current blood-only mask (symmetric ROI), `eFE`/`ePE` may
   be weakly identified because the mask is approximately circular and symmetric about the
   velocity peak. Including wall pixels (where the asymmetry is most visible) in the fit
   might help.

2. **Joint m(v) fit with asymmetry:** the magnitude profile `m(v(r,theta))` will also be
   asymmetric if the PSF blurs differently in each direction. The joint fit should remain
   valid since `magFit(v)` is purely a function of velocity, not position.

3. **Regularisation:** adding a small penalty on `|e|²` (ridge on eFE, ePE) would improve
   conditioning if the data don't strongly constrain the wall offset.

4. **Calling convention change:** if `FEoffset`/`PEoffset` are "inert" in the sfit (Strategy A),
   the caller must still supply them to compute `r` and `theta`. A helper function
   `eval_velFit(velFit, FEgrid, PEgrid)` that extracts the offsets and computes
   `r, theta` internally would make the call sites cleaner.

---

## Related files

- `fitMagVelProfile.m` — current symmetric implementation; asymmetric extension goes here
- `velocity_func_radius.m` — current symmetric model; `velocity_func_radius_asym.m` analogous
- `fitVelProfile.m` — symmetric-only velocity fit (separate entry point)
- `doIt.m` — call sites: "Radial profiles fits" section, "Plot matched simulation summary" section
