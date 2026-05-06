function [velFit, magFit] = fitMagVelProfile(r, v, m, mNoFlow, R, fitType, polyOrder, offset, FEPE)
% Fit v(r) and m(v) jointly or sequentially to flow and no-flow radial profile data.
%   r       — Nx1 radial positions from nominal centre [mm]; rGrid(mask)
%   v, m    — Nx1 velocity [cm/s] and magnitude of flow pixels
%   mNoFlow — Mx1 magnitudes of no-flow pixels; fixes B = mean(mNoFlow)
%   R       — vessel radius start point [mm] (required when offset=true)
%   fitType — 'joint' (default) or 'sequential'
%   polyOrder — polynomial order for m(v) (default 2; supports 1–4)
%   offset  — []   : no centre offset (default)
%             true : fit (FEoffset,PEoffset) as free parameters; FEPE (arg 9) required
%   FEPE    — Nx2 [FEgrid(mask), PEgrid(mask)]: col1 = FE (row axis, vertical
%             in imagesc display), col2 = PE (col axis, horizontal); required
%             when offset=true
%
%   Models:  v(r)    = Vmax*(1-(r/R)^2)           velocity_func_radius
%            m(v)    = B + C1*v + ... + Cn*v^n    magnitude_func_velocity
%
%   velFit calling convention:
%     offset=[]:   velFit(r)         — 1D cfit; r = rGrid(mask)
%     offset=true: velFit(r, theta)  — 2D sfit
%                  theta = atan2(PE, FE) [rad]; theta=0 along +FE axis
%                  (vertical in imagesc display, i.e. downward for typical FEpos)
%                  Use pGrid = atan2(PEgrid, FEgrid) at the call site.
%   FEoffset = offset along FE axis (vertical in display); velFit.FEoffset
%   PEoffset = offset along PE axis (horizontal in display); velFit.PEoffset
if ~exist('offset'   ,'var') || isempty(offset);    offset    = []; end
if ~exist('fitType'  ,'var') || isempty(fitType);   fitType   = 'joint'; end
if ~exist('polyOrder','var') || isempty(polyOrder); polyOrder = 2; end
if exist('R','var') && ~isempty(R); R = double(R); end

r = double(r(:)); v = double(v(:)); m = double(m(:)); mNoFlow = double(mNoFlow(:));
B = mean(mNoFlow);

fitOff = isequal(offset, true);
if fitOff
    if nargin < 9 || isempty(FEPE)
        error('fitMagVelProfile: offset=true requires FEPE as 9th argument [FEgrid(mask), PEgrid(mask)].');
    end
    FEPE = double(FEPE);
end

% v(r) fittype — 1D, no offset
ft_vel = fittype(@(Vmax, R, r) velocity_func_radius(r, Vmax, R), ...
    'independent', 'r', 'coefficients', {'Vmax', 'R'});

% v(r,theta) fittype — 2D sfit, polar coords; offset correction applied internally
ft_vel_off = fittype(@(Vmax, R, FEoffset, PEoffset, r, theta) ...
    velocity_func_radius(sqrt((r.*cos(theta)-FEoffset).^2+(r.*sin(theta)-PEoffset).^2), Vmax, R), ...
    'independent', {'r', 'theta'}, 'coefficients', {'Vmax', 'R', 'FEoffset', 'PEoffset'});

% m(v) fittype — built dynamically for the requested polynomial order
switch polyOrder
    case 1
        ft_mag = fittype(@(B, C1, v) magnitude_func_velocity(v, B, C1), ...
            'independent', 'v', 'coefficients', {'B', 'C1'});
    case 2
        ft_mag = fittype(@(B, C1, C2, v) magnitude_func_velocity(v, B, C1, C2), ...
            'independent', 'v', 'coefficients', {'B', 'C1', 'C2'});
    case 3
        ft_mag = fittype(@(B, C1, C2, C3, v) magnitude_func_velocity(v, B, C1, C2, C3), ...
            'independent', 'v', 'coefficients', {'B', 'C1', 'C2', 'C3'});
    case 4
        ft_mag = fittype(@(B, C1, C2, C3, C4, v) magnitude_func_velocity(v, B, C1, C2, C3, C4), ...
            'independent', 'v', 'coefficients', {'B', 'C1', 'C2', 'C3', 'C4'});
    otherwise
        error('fitMagVelProfile: polyOrder must be 1–4.');
end
sp_mag = [B, zeros(1, polyOrder)];

switch fitType

    case 'sequential'
        [velFit_, ~]   = fitVelProfile(r, v, R);
        velFit         = cfit(ft_vel, velFit_.Vmax, velFit_.R);
        v_pred         = velocity_func_radius(r, velFit_.Vmax, velFit_.R);
        v_all          = [v_pred(:);  zeros(numel(mNoFlow), 1)];
        m_all          = [m(:);       mNoFlow(:)              ];
        magFit         = fit(v_all, m_all, ft_mag, 'StartPoint', sp_mag);

    case 'joint'
        sv = max(std(v),            eps);
        sm = max(std([m; mNoFlow]), eps);
        [velFit0, ~] = fitVelProfile(r, v, R);

        if fitOff
            % theta = [Vmax, R, B, C1..Cn, FEoffset, PEoffset]
            R0     = double(R);  % use input radius (velFit0.R can be unreliable)
            theta0 = [velFit0.Vmax, R0, B, zeros(1, polyOrder), 0, 0];
            lb     = [0,    1e-6,  -inf(1, polyOrder+1), -R0/4, -R0/4];
            ub     = [inf,  2*R0,   inf(1, polyOrder+1),  R0/4,  R0/4];
            opts   = optimoptions('lsqnonlin', 'Display', 'off');
            theta  = lsqnonlin(@(th) residuals_joint_off(th, FEPE, v, m, mNoFlow, sv, sm), ...
                               theta0, lb, ub, opts);
            % sfit uses definition order: Vmax, R, FEoffset, PEoffset
            velFit = sfit(ft_vel_off, theta(1), theta(2), theta(end-1), theta(end));
        else
            theta0 = [velFit0.Vmax, velFit0.R, B, zeros(1, polyOrder)];
            lb     = [0, 1e-6, -inf(1, polyOrder+1)];
            opts   = optimoptions('lsqnonlin', 'Display', 'off');
            theta  = lsqnonlin(@(th) residuals_joint(th, r, v, m, mNoFlow, sv, sm), ...
                               theta0, lb, [], opts);
            velFit = cfit(ft_vel, theta(1), theta(2));
        end
        C_vals = num2cell(theta(3 : 3 + polyOrder));
        magFit = cfit(ft_mag, C_vals{:});

    otherwise
        error('fitMagVelProfile: unknown fitType ''%s''. Choose: sequential | joint.', fitType);

end
end


function res = residuals_joint(theta, r, v, m, mNoflow, sv, sm)
Vmax = theta(1);  R = theta(2);
C    = num2cell(theta(3:end));
v_pred  = velocity_func_radius(r, Vmax, R);
res_v   = (v       - v_pred)                                 / sv;
res_m   = (m       - magnitude_func_velocity(v_pred, C{:})) / sm;
res_mNF = (mNoflow - magnitude_func_velocity(0,      C{:})) / sm;
res = [res_v(:); res_m(:); res_mNF(:)];
end


function res = residuals_joint_off(theta, FEPE, v, m, mNoflow, sv, sm)
Vmax = theta(1);  R = theta(2);
C    = num2cell(theta(3 : end-2));
FEoffset = theta(end-1);  PEoffset = theta(end);
r    = sqrt((FEPE(:,1) - FEoffset).^2 + (FEPE(:,2) - PEoffset).^2);
v_pred  = velocity_func_radius(r, Vmax, R);
res_v   = (v       - v_pred)                                 / sv;
res_m   = (m       - magnitude_func_velocity(v_pred, C{:})) / sm;
res_mNF = (mNoflow - magnitude_func_velocity(0,      C{:})) / sm;
res = [res_v(:); res_m(:); res_mNF(:)];
end


%% -------------------------------------------------------------------------
% DEVELOPMENT NOTE: asymmetric radial profile
% Full description: asymmetricProfile.md (same directory)
%
% Idea: separate the velocity-peak centre from the wall-circle centre.
%   Current:  peak at (FEoffset, PEoffset), wall circle also centred there → symmetric.
%   Extended: wall circle centred at (FEoffset+eFE, PEoffset+ePE) → asymmetric.
%
% Key result — angle-dependent effective radius:
%   Given wall offset (eFE,ePE) relative to the peak, a ray cast from the peak
%   in direction theta_v hits the wall circle (radius R) at distance:
%
%     p        = eFE*cos(theta_v) + ePE*sin(theta_v)    % projection onto ray
%     R_eff    = p + sqrt(R^2 - eFE^2 - ePE^2 + p^2)   % exact
%     R_eff   ~= R + p                                   % first-order (|e|<<R)
%
%   Asymmetric model (drop-in replacement for velocity_func_radius):
%     v(r_v, theta_v) = Vmax * max(0, 1 - (r_v/R_eff(theta_v))^2)
%
%   New parameters: eFE, ePE  (bounds: ±R0/4, same logic as FEoffset/PEoffset)
%   Total free velocity params: Vmax, R, FEoffset, PEoffset, eFE, ePE  (6)
%
% Preferred implementation strategy (Strategy A in asymmetricProfile.md):
%   - New helper: velocity_func_radius_asym(r, theta, Vmax, R, eFE, ePE)
%   - New fittype ft_vel_asym with independent {'r','theta'} and 6 coefficients;
%     FEoffset/PEoffset stored as inert coefficients (same trick as current sfit).
%   - New residuals_joint_asym taking FEPE, computing r_v/theta_v from FEoffset/PEoffset.
%   - Calling convention unchanged: velFit(r, theta) evaluated at offset-corrected coords.
%
% Caution — identifiability: eFE/ePE are weakly constrained by a symmetric blood-only
%   mask; wall pixels or a mask penalty may be needed. See asymmetricProfile.md §Open questions.
% -------------------------------------------------------------------------
