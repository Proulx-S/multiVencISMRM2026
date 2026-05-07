function [maskBloodOnly, maskWallOnly, maskNonBloodOnly, maskTissueOnly, maskWallLowMag] = makeVesselMasks(FEgrid, PEgrid, FEspacing, PEspacing, ID, OD, m)
% Pixel membership masks for a circular cross-section vessel.
%   FEgrid, PEgrid       — 2D coordinate grids (mm), centred on vessel axis
%   FEspacing, PEspacing — pixel spacings (mm)
%   ID, OD               — inner and outer diameters (mm)
%   m                    — magnitude image; required only for maskWallLowMag
%
%   maskBloodOnly    — pixel entirely inside inner circle  (d_far  < ID/2)
%   maskWallOnly     — pixel entirely within wall annulus  (d_near > ID/2 & d_far < OD/2)
%   maskNonBloodOnly — nearest edge outside inner circle   (d_near > ID/2)
%   maskTissueOnly   — pixel entirely outside outer circle (d_near > OD/2)
%   maskWallLowMag   — m < min tissue magnitude; empty if m not supplied

dFE    = abs(FEgrid);
dPE    = abs(PEgrid);
d_far  = sqrt((dFE + FEspacing/2).^2 + (dPE + PEspacing/2).^2);
d_near = sqrt(max(0, dFE - FEspacing/2).^2 + max(0, dPE - PEspacing/2).^2);

maskBloodOnly    = d_far  < ID/2;
maskWallOnly     = d_near > ID/2 & d_far < OD/2;
maskNonBloodOnly = d_near > ID/2;
maskTissueOnly   = d_near > OD/2;

if nargin >= 7 && ~isempty(m)
    maskWallLowMag = m < min(m(maskTissueOnly));
else
    maskWallLowMag = [];
end
