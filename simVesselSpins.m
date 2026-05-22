function [magMap,vMap,pVessel,pSim,pMri] = simVesselSpins(pVessel, pSim, pMri)

% Handle 'cylinder3D' profile: fully 3D geometry with per-isochromat inflow signal.
% Returns early so the 2D-only machinery below is untouched.
if ischar(pVessel.profile) && strcmp(pVessel.profile, 'cylinder3D')
    [magMap, vMap, pVessel, pSim, pMri] = simCylinder3D(pVessel, pSim, pMri);
    return;
end

% ---- original simVesselSpins code (unchanged below this line) ----

% Define radial coordinates (relative to vessel center)
[gridFE, gridPE] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);  % dim1=FE, dim2=PE
rGrid = sqrt((gridFE - pVessel.posFE).^2 + (gridPE - pVessel.posPE).^2);

% Define compartment masks
% Elliptical vessel if pVessel.AR and pVessel.alpha are provided; circular otherwise.
R_lumen = pVessel.ID / 2;  % semi-major axis [mm]
if isfield(pVessel,'AR') && ~isempty(pVessel.AR) && pVessel.AR ~= 1
    alpha_v = pVessel.alpha;
    AR_v    = pVessel.AR;
    dFE = gridFE - pVessel.posFE;
    dPE = gridPE - pVessel.posPE;
    u   =  dPE.*cos(alpha_v) + dFE.*sin(alpha_v);   % along major axis
    w   = -dPE.*sin(alpha_v) + dFE.*cos(alpha_v);   % along minor axis
    ellipseDist = sqrt(u.^2 + (AR_v .* w).^2);      % = R_lumen at ellipse boundary
    pVessel.mask.lumen    = ellipseDist <= R_lumen;
    pVessel.mask.wall     = ellipseDist >  R_lumen & ellipseDist <= R_lumen + pVessel.WT;
    pVessel.mask.surround = ~pVessel.mask.lumen & ~pVessel.mask.wall;
else
    pVessel.mask.lumen    = rGrid <= R_lumen;
    pVessel.mask.wall     = rGrid >  R_lumen & rGrid <= R_lumen + pVessel.WT;
    pVessel.mask.surround = rGrid >  R_lumen + pVessel.WT;
end
assert(all(pVessel.mask.lumen(:) + pVessel.mask.wall(:) + pVessel.mask.surround(:) == 1), ...
    'pVessel.mask: lumen, wall, surround must be non-overlapping and cover all spins');

% Define spin velocity map
if ischar(pVessel.profile)
    vMap = getVelMap(rGrid, pVessel.ID, pVessel.profile, pVessel.PD); % [cm/s]
    if ~isempty(pVessel.vMax) && isempty(pVessel.vMean)
        vMap = scale2maxVel(vMap, pVessel.vMax); % to the desired maximum velocity
    elseif ~isempty(pVessel.vMean) && isempty(pVessel.vMax)
        vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
    elseif strcmp(pVessel.profile,'parabolic1') && ~isempty(pVessel.vMax) && ~isempty(pVessel.vMean) && pVessel.vMax/2==pVessel.vMean
        vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
    else
        error('Either pVessel.vMax or pVessel.vMean must be specified');
    end
elseif isnumeric(pVessel.profile)
    vMap = zeros(size(rGrid));
    vMap(:) = pVessel.profile;
else
    dbstack; error('Invalid vessel profile');
end


% MR signal magnitude
% vessel lumen signal (flowing)
if isempty(pVessel.S.lumen)
    switch pVessel.profile
        case 'plug'
            [Mz_vMean,pMri] = getMz_ss(          pMri,pMri.relax.blood,pVessel.vMean);
            [Mxy_vMax,pMri] = getMxy_ss(Mz_vMean,pMri,pMri.relax.blood              );
            pVessel.S.lumen = Mxy_vMax;
        case {'parabolic','parabolic1'}
            [Mz ,pMri] = getMz_ss(    pMri,pMri.relax.blood,vMap(pVessel.mask.lumen));
            [Mxy,pMri] = getMxy_ss(Mz,pMri,pMri.relax.blood                         );
            pVessel.S.lumen = Mxy;
        otherwise
            dbstack; error('Invalid vessel profile');
    end
end
% vessel surround (static)
if isempty(pVessel.S.surround)
    pVessel.S.surround = getMxy_ss(getMz_ss(pMri,pMri.relax.GM),pMri,pMri.relax.GM);
end
% map signal magnitude (divide by nSpinPerVox so summing spins within voxels gives S, the hypothetical measured signal if the voxel was single-compartment)
magMap = zeros(size(rGrid));
magMap(pVessel.mask.lumen)    = pVessel.S.lumen    ./pSim.nSpinPerVox;
magMap(pVessel.mask.wall)     = pVessel.S.wall     ./pSim.nSpinPerVox;
magMap(pVessel.mask.surround) = pVessel.S.surround ./pSim.nSpinPerVox;


% Precompute montecarlo tessalation
if pSim.monteCarloN > 0 && (~isfield(pSim,'monteCarloShiftFE') || ~isfield(pSim,'monteCarloShiftPE') || isempty(pSim.monteCarloShiftFE) || isempty(pSim.monteCarloShiftPE))
    nSpinFE = pSim.spinGrid.matFE / pSim.voxGrid.matFE;  % spins per voxel in FE
    shiftFE = (1:nSpinFE)-nSpinFE/2-0.5;
    nSpinPE = pSim.spinGrid.matPE / pSim.voxGrid.matPE;  % spins per voxel in PE
    shiftPE = (1:nSpinPE)-nSpinPE/2-0.5;
    % find all possible combination of FE and PE shifts
    [idx1, idx2] = ndgrid(1:length(shiftFE), 1:length(shiftPE));
    idx = [idx1(:), idx2(:)];
    % remove the no-shift combination since it is always done before
    idx(all(idx==[find(shiftFE==0) find(shiftPE==0)],2),:) = [];
    % shuffle
    if pSim.monteCarloN~=inf
        idx = idx(randperm(size(idx,1),pSim.monteCarloN),:);
    end
    pSim.monteCarloShiftFE = shiftFE(idx(:,1));
    pSim.monteCarloShiftPE = shiftPE(idx(:,2));
end


% =========================================================================
% Local functions
% =========================================================================

function [magMap, vMap, pVessel, pSim, pMri] = simCylinder3D(pVessel, pSim, pMri)
% Full 3D cylinder geometry for the 'cylinder3D' profile.
%
% The vessel is an infinite cylinder defined by:
%   n_hat  = pVessel.n_hat   [nx; ny; nz]  unit axis vector
%   posFE  = pVessel.posFE   [mm]           FE center offset
%   posPE  = pVessel.posPE   [mm]           PE center offset
%   posSLC = pVessel.posSLC  [mm]           slice-direction offset (default 0)
%   ID     = pVessel.ID      [mm]           lumen inner diameter
%   WT     = pVessel.WT      [mm]           wall thickness (OD/2 = ID/2 + WT)
%   Vmax   = pVessel.Vmax    [cm/s]         peak parabolic velocity
%   A      = pVessel.A       [a.u.]         lumen signal amplitude scale
%
% Inflow model: instantaneous Mz after n_k pulses (per-isochromat z-position),
% computed via getMz_n + getMxy_ss. No effective-slice-thickness approximation.

n_hat  = pVessel.n_hat;   % [nx; ny; nz]
nx = n_hat(1);   ny = n_hat(2);   nz = n_hat(3);
if isfield(pVessel,'posSLC');  cx_SLC = pVessel.posSLC;  else;  cx_SLC = 0;  end

% --- Z coordinate grid for the slice direction ---
% Sample Z at the same count as sub-spins per voxel in FE (e.g. 7 for a 7×7 grid).
% This is coarser than the FE/PE spin spacing but sufficient to resolve n_k(z).
% Slab thickness comes from pMri.sliceThickness; matSlice=1 assumed.
d_slab_zg  = pMri.sliceThickness;   % [mm]
nSpSLC     = round(numel(pSim.spinGrid.coorFE) / numel(pSim.voxGrid.coorFE));
nSpSLC     = max(1, nSpSLC);
if mod(nSpSLC, 2) == 0; nSpSLC = nSpSLC + 1; end   % make odd
dSLC       = d_slab_zg / nSpSLC;
coorZ      = ((-(nSpSLC-1)/2) : ((nSpSLC-1)/2)) .* dSLC;   % cell-centered [mm]

% Effective spin count per voxel including Z dimension
nSpinPerVox_3D = pSim.nSpinPerVox * nSpSLC;

% --- 3D spin coordinate grids ---
% Sizes: [nTotalFE, nTotalPE, nSpSLC]
[gridFE, gridPE, gridZ] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE, coorZ);

% --- 3D cylinder r_perp ---
fe_rel = gridFE - pVessel.posFE;
pe_rel = gridPE - pVessel.posPE;
z_rel  = gridZ  - cx_SLC;
dot_n  = fe_rel.*nx + pe_rel.*ny + z_rel.*nz;
r_perp = sqrt(max(0, fe_rel.^2 + pe_rel.^2 + z_rel.^2 - dot_n.^2));

% --- Compartment masks (3D) ---
R_blood = pVessel.ID / 2;
R_wall  = R_blood + pVessel.WT;
lumen3  = r_perp <= R_blood;
wall3   = r_perp >  R_blood & r_perp <= R_wall;
surr3   = ~lumen3 & ~wall3;

% --- Parabolic velocity in cylinder coordinates ---
Vmax  = pVessel.Vmax;
A_scl = pVessel.A;
v3    = zeros(size(r_perp));
v3(lumen3) = Vmax .* (1 - (r_perp(lumen3) ./ R_blood).^2);

% --- Per-isochromat inflow signal (instantaneous Mz from Z position) ---
d_slab = pMri.sliceThickness;   % [mm]  true slab thickness

% Path along cylinder axis from slab entry (Z = -d_slab/2 + cx_SLC) to current Z position.
% Assumes nz > 0 (blood flows in +z direction through the slab).
z_lumen    = z_rel(lumen3);
d_trav     = (z_lumen + d_slab/2) ./ max(nz, eps);   % [mm] path from slab entry
v_lumen    = v3(lumen3);                              % [cm/s]
% Convert d_trav mm → cm before dividing by v [cm/s] and TR [s].
n_k        = ceil(max(0, d_trav/10) ./ max(v_lumen, eps) ./ pMri.TR);
n_k(v_lumen < eps) = Inf;   % stationary spins → steady-state saturation

% getMz_n: instantaneous Mz after exactly n_k pulses from fully-relaxed M0.
[Mz_k, pMri] = getMz_n(pMri, pMri.relax.blood, n_k);
% getMxy_ss: apply sin(FA) flip and T2* decay to get Mxy.
[Mxy_k, pMri] = getMxy_ss(Mz_k, pMri, pMri.relax.blood);

% --- Surround signal ---
if ~isempty(pVessel.S.surround)
    S_surr = pVessel.S.surround;   % user-specified (e.g. S_tissue from fit)
else
    % Auto-compute from GM relaxation (default for non-fitting use)
    [Mz_gm, pMri]  = getMz_ss(pMri, pMri.relax.GM, 0);
    [S_surr, pMri] = getMxy_ss(Mz_gm, pMri, pMri.relax.GM);
end

% --- Build 3D signal array and Z-sum to 2D ---
% S3: [nTotalFE, nTotalPE, nSpSLC]
S3 = zeros(size(r_perp));
S3(lumen3) = A_scl .* Mxy_k;   % per-lumen-spin inflow signal
S3(surr3)  = S_surr;           % scalar, broadcast over surround spins
% wall3 contribution = 0 (already zero)

% Sum over Z dimension → [nTotalFE, nTotalPE]
S2 = sum(S3, 3);

% --- 2D masks (projection over Z: spin belongs to lumen if ANY Z slice is in lumen) ---
lumen2 = any(lumen3, 3);
wall2  = any(wall3,  3) & ~lumen2;
surr2  = ~lumen2 & ~wall2;
assert(all(lumen2(:) + wall2(:) + surr2(:) == 1), ...
    'simCylinder3D: 2D masks must be non-overlapping and exhaustive');

% --- Assign 2D pVessel fields ---
pVessel.mask.lumen    = lumen2;
pVessel.mask.wall     = wall2;
pVessel.mask.surround = surr2;
pVessel.S.lumen       = S2(lumen2);    % Z-sum per 2D lumen spin
pVessel.S.wall        = 0;
pVessel.S.surround    = S_surr;        % keep as scalar (magMap uses nSpinPerVox_3D below)

% --- 2D magMap (per-spin, divided by 3D nSpinPerVox) ---
% By using nSpinPerVox_3D in the denominator, the sum over the 2D spin grid
% per voxel correctly recovers the 3D mean signal:
%   sum_{2D spins} magMap = sum_{2D} S2 / (nSpFE*nSpPE*nSpSLC)
%                         = sum_{2D} sum_Z(S3) / (nSpFE*nSpPE*nSpSLC)
%                         = mean_{2D,Z}(S3)   (3D mean per voxel)
magMap = zeros(size(gridFE, 1), size(gridFE, 2));   % [nTotalFE, nTotalPE]
magMap(lumen2) = S2(lumen2)         ./ nSpinPerVox_3D;
% wall2 contribution is 0
magMap(surr2)  = (S_surr .* nSpSLC) ./ nSpinPerVox_3D;   % nSpSLC factor cancels in aggregation

% --- 2D velocity map (at Z=0, slab midplane) ---
fe_2d = gridFE(:,:,1) - pVessel.posFE;
pe_2d = gridPE(:,:,1) - pVessel.posPE;
z_2d  = 0 - cx_SLC;   % scalar
dot_n_2d  = fe_2d.*nx + pe_2d.*ny + z_2d.*nz;
r_perp_2d = sqrt(max(0, fe_2d.^2 + pe_2d.^2 + z_2d.^2 - dot_n_2d.^2));
vMap  = Vmax .* max(0, 1 - (r_perp_2d ./ R_blood).^2);
vMap(r_perp_2d > R_blood) = 0;
