function [velFit, magFit, velFit1d, fitInfo] = fitMagVelProfile(r, p, v, m, mNoFlow, R, fitType, polyOrder, offset, B_init_norm, theta0_phys, maskBlood)
% Fit v(r,p) and m(v) jointly or sequentially to flow and no-flow radial profile data.
%   r           — Nx1 radial distance from nominal centre [mm]; rGrid(mask)
%   p           — Nx1 polar angle [rad]; pGrid(mask)
%                 p = -atan2(FEgrid, PEgrid); p=0 → +PE axis (right in imagesc display)
%   v, m        — Nx1 velocity [cm/s] and magnitude of flow pixels
%   mNoFlow     — Mx1 magnitudes of no-flow pixels; fixes B = mean(mNoFlow)
%                 Pass [] to let B be a free parameter (in-vivo mode).
%   R           — vessel radius start point [mm] (required when offset=true)
%   fitType     — 'joint' (default) or 'sequential'
%   polyOrder   — polynomial order for m(v) (default 2; supports 1–4)
%   offset      — []   : no centre offset (default)
%                 true : fit (FEoffset,PEoffset) as free parameters
%   B_init_norm — Initial B value normalised to mean(m); used only when mNoFlow=[].
%                 Default 0.7. Derive from getMz_ss/getMxy_ss for physics-based prior.
%
%   Models:  v(r)    = Vmax*(1-(r/R)^2)           velocity_func_radius
%            m(v)    = B + C1*v + ... + Cn*v^n    magnitude_func_velocity
%
%   velFit calling convention:
%     offset=[]:   velFit(r, p)      — 2D sfit; centre offset fixed at (0,0)
%     offset=true: velFit(r, p)      — 2D sfit; FEoffset/PEoffset fitted
%   r*cos(p) = PE component, r*sin(p) = -FE component  (p=0 → right, CCW positive)
%   FEoffset = offset along FE axis (vertical); velFit.FEoffset
%   PEoffset = offset along PE axis (horizontal); velFit.PEoffset
fitInfo = struct();   % populated for 'joint'+'ellipse' and 'complex' when nargout>=4
if ~exist('offset'      ,'var') || isempty(offset);       offset      = []; end
if ~exist('maskBlood'   ,'var');                           maskBlood   = []; end
if ~exist('fitType'     ,'var') || isempty(fitType);      fitType     = 'joint'; end
if ~exist('polyOrder'   ,'var') || isempty(polyOrder);    polyOrder   = 2; end
if ~exist('B_init_norm' ,'var') || isempty(B_init_norm);  B_init_norm = 0.7; end
if exist('R','var') && ~isempty(R); R = double(R); end

r = double(r(:)); p = double(p(:));

complexMode = isequal(fitType, 'complex');

if complexMode
    % v arg = s_all (N×K complex, phase-corrected; N = blood pixels, or all ROI pixels)
    % m arg = venc_list (1×K cm/s)
    % maskBlood: logical N×1 — true=blood, false=tissue. [] means all-blood.
    s_all     = double(v);
    venc_list = double(m(:)');
    [~, kHigh] = max(venc_list);
    hasTissue = ~isempty(maskBlood) && any(~maskBlood(:));
    mb        = hasTissue && true;   % flag shorthand
    if hasTissue
        bl = maskBlood(:);
    else
        bl = true(size(s_all, 1), 1);
    end
    % mScale: mean blood signal at highest VENC — matches Fit A scaling
    mScale  = mean(abs(s_all(bl, kHigh)));
    s_scale = mScale;
    s_n     = s_all / s_scale;
    % s_std normalised across blood pixels only (tissue contrast would inflate it)
    s_std   = max(std(abs(s_n(bl, :)), 0, 'all'), eps);
    v_init  = phase2vel(angle(s_all(bl, kHigh)), vencToM1(venc_list(kHigh)));
    B       = B_init_norm;
    noFlowFree = true;
else
    v = double(v(:)); m = double(m(:));
    % Normalise magnitudes internally so polynomial coefficients are O(1).
    % magFit is un-scaled on output.
    % When mNoFlow=[] (in-vivo mode), B is a free parameter initialized at B_init_norm.
    noFlowFree = isempty(mNoFlow);
    if noFlowFree
        mScale = mean(m);
        m      = m / mScale;
        B      = B_init_norm;
    else
        mNoFlow = double(mNoFlow(:));
        mScale  = mean(mNoFlow);
        m       = m       / mScale;
        mNoFlow = mNoFlow / mScale;
        B = 1;  % = mean(mNoFlow) after normalisation
    end
end

fitOff     = isequal(offset, true);
fitEllipse = isequal(offset, 'ellipse');
fitAsym    = isequal(offset, 'asym');
if (fitOff || fitEllipse || fitAsym) && (nargin < 6 || isempty(R))
    error('fitMagVelProfile: offset=true/''asym'' requires R as 6th argument.');
end

% v(r) fittype — 1D, no offset
ft_vel = fittype(@(Vmax, R, r) velocity_func_radius(r, Vmax, R), ...
    'independent', 'r', 'coefficients', {'Vmax', 'R'});

% v(r,p) fittype — 2D sfit, symmetric circle with centre offset
% p=0 → +PE (right in display): r*cos(p)=PE component, r*sin(p)=-FE component
ft_vel_off = fittype(@(Vmax, R, FEoffset, PEoffset, r, p) ...
    velocity_func_radius(sqrt((r.*cos(p)-PEoffset).^2+(-r.*sin(p)-FEoffset).^2), Vmax, R), ...
    'independent', {'r', 'p'}, 'coefficients', {'Vmax', 'R', 'FEoffset', 'PEoffset'});

% v(r,p) fittype — 2D sfit, elliptical wall centred at velocity peak
% Semi-major axis R along angle alpha from PE axis; semi-minor axis R/AR (AR>=1)
ft_vel_ellipse = fittype( ...
    @(Vmax, R, AR, alpha, FEoffset, PEoffset, r, p) ...
        velocity_func_ellipse(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset), ...
    'independent', {'r','p'}, ...
    'coefficients', {'Vmax','R','AR','alpha','FEoffset','PEoffset'});

% v(r,p) fittype — 2D sfit, asymmetric elliptical wall (peak offset from wall centre)
% Velocity peak at (FEoffset,PEoffset); wall ellipse centred at (FEoffset+eFE,PEoffset+ePE)
% Semi-major axis R along angle alpha from PE axis; semi-minor axis R/AR (AR>=1)
ft_vel_asym = fittype( ...
    @(Vmax, R, AR, alpha, FEoffset, PEoffset, eFE, ePE, r, p) ...
        velocity_func_ellipse_asym(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset, eFE, ePE), ...
    'independent', {'r','p'}, ...
    'coefficients', {'Vmax','R','AR','alpha','FEoffset','PEoffset','eFE','ePE'});

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
        velFit         = sfit(ft_vel_off, velFit_.Vmax, velFit_.R, 0, 0);
        velFit1d       = cfit(ft_vel,     velFit_.Vmax, velFit_.R);
        v_pred         = velocity_func_radius(r, velFit_.Vmax, velFit_.R);
        v_all          = [v_pred(:);  zeros(numel(mNoFlow), 1)];
        m_all          = [m(:);       mNoFlow(:)              ];
        magFit_n     = fit(v_all, m_all, ft_mag, 'StartPoint', sp_mag);
        coeffs_scaled = num2cell(coeffvalues(magFit_n) * mScale);
        magFit       = cfit(ft_mag, coeffs_scaled{:});

    case 'joint'
        sv = max(std(v),            eps);
        if noFlowFree
            sm = max(std(m), eps);
        else
            sm = max(std([m; mNoFlow]), eps);
        end
        [velFit0, ~] = fitVelProfile(r, v, R);
        mNoFlow_pass = mNoFlow;  % [] when noFlowFree, Mx1 otherwise

        opts = optimoptions('lsqnonlin', 'Display', 'off');
        R0   = double(R);

        if fitEllipse
            % Internal theta = [Vmax, R, e1, e2, B, C1..Cn, FEoffset, PEoffset]
            % (e1,e2) = (AR-1)*[cos(2*alpha), sin(2*alpha)] — Cartesian ellipse coords
            % FEoffset and PEoffset fixed at 0 (lb==ub enforces equality).
            theta0 = [velFit0.Vmax, R0, 0, 0, B, zeros(1,polyOrder), 0, 0];
            [poly_lb, poly_ub] = polyBounds(polyOrder);
            lb     = [0,    1e-6,  -1, -1, 0,   poly_lb, 0, 0];
            ub     = [inf,  2*R0,   1,  1, inf, poly_ub, 0, 0];
            theta  = lsqnonlin(@(th) residuals_joint_ellipse(th, r, p, v, m, mNoFlow_pass, sv, sm), ...
                               theta0, lb, ub, opts);
            [AR_f, alpha_f] = e_to_ar_alpha(theta(3), theta(4));
            velFit = sfit(ft_vel_ellipse, theta(1),theta(2),AR_f,alpha_f,theta(end-1),theta(end));
            c_start = 5;
            if nargout >= 4
                sc = [1,1,1,1, repmat(mScale,1,polyOrder+1), 1,1];
                mag_names = [{'B'}, arrayfun(@(i)sprintf('C%d',i),1:polyOrder,'UniformOutput',false)];
                fitInfo.names  = [{'Vmax','R','e1','e2'}, mag_names, {'FEoffset','PEoffset'}];
                fitInfo.units  = [{'cm/s','mm','-','-'}, repmat({'a.u.'},1,polyOrder+1), {'mm','mm'}];
                fitInfo.fixed  = [false(1, 5+polyOrder), true, true];
                fitInfo.theta0 = theta0 .* sc;
                fitInfo.lb     = [0, 1e-6, -1, -1, 0, poly_lb, 0, 0];
                fitInfo.ub     = [inf, 2*R0, 1, 1, inf, poly_ub, 0, 0];
                fitInfo.theta  = theta  .* sc;
                fitInfo.derived.names  = {'AR', 'alpha'};
                fitInfo.derived.units  = {'-', 'rad'};
                fitInfo.derived.theta0 = [1, 0];
                fitInfo.derived.theta  = [AR_f, alpha_f];
            end
        elseif fitAsym
            % Internal theta = [Vmax, R, e1, e2, B, C1..Cn, FEoffset, PEoffset, eFE, ePE]
            theta0 = [velFit0.Vmax, R0, 0, 0, B, zeros(1,polyOrder), 0, 0, 0, 0];
            [poly_lb, poly_ub] = polyBounds(polyOrder);
            lb     = [0,    1e-6, -1, -1, 0,   poly_lb, -R0/2,-R0/2, -R0/4,-R0/4];
            ub     = [inf,  2*R0,  1,  1, inf, poly_ub,  R0/2, R0/2,  R0/4, R0/4];
            theta  = lsqnonlin(@(th) residuals_joint_asym(th, r, p, v, m, mNoFlow_pass, sv, sm), ...
                               theta0, lb, ub, opts);
            [AR_f, alpha_f] = e_to_ar_alpha(theta(3), theta(4));
            velFit = sfit(ft_vel_asym, theta(1),theta(2),AR_f,alpha_f, ...
                          theta(end-3),theta(end-2),theta(end-1),theta(end));
            c_start = 5;
        elseif fitOff
            % theta = [Vmax, R, B, C1..Cn, FEoffset, PEoffset]
            theta0 = [velFit0.Vmax, R0, B, zeros(1,polyOrder), 0, 0];
            lb     = [0,    1e-6,  0,  -inf(1,polyOrder), -R0/4, -R0/4];
            ub     = [inf,  2*R0,  inf(1,polyOrder+1),     R0/4,  R0/4];
            theta  = lsqnonlin(@(th) residuals_joint_off(th, r, p, v, m, mNoFlow_pass, sv, sm), ...
                               theta0, lb, ub, opts);
            velFit = sfit(ft_vel_off, theta(1), theta(2), theta(end-1), theta(end));
            c_start = 3;
        else
            theta0 = [velFit0.Vmax, velFit0.R, B, zeros(1,polyOrder)];
            lb     = [0, 1e-6, 0, -inf(1,polyOrder)];
            theta  = lsqnonlin(@(th) residuals_joint(th, r, v, m, mNoFlow_pass, sv, sm), ...
                               theta0, lb, [], opts);
            velFit = sfit(ft_vel_off, theta(1), theta(2), 0, 0);
            c_start = 3;
        end
        velFit1d = cfit(ft_vel, theta(1), theta(2));  % 1D: v(r) along major axis, no offset
        C_vals   = num2cell(theta(c_start : c_start+polyOrder) * mScale);
        magFit   = cfit(ft_mag, C_vals{:});

    case 'complex'
        % Complex multi-VENC residual.
        % Blood pixels: theta = [Vmax,R,AR,alpha, B,C1..Cn (,Bt if tissue present)]
        % Tissue pixels (maskBlood=false): predicted signal = Bt (constant across VENCs)
        % FEoffset/PEoffset fixed at 0.
        if ~isequal(offset, 'ellipse')
            error('fitMagVelProfile: fitType=''complex'' currently requires offset=''ellipse''.');
        end
        opts = optimoptions('lsqnonlin', 'Display', 'off');
        R0       = double(R);
        venc_row = reshape(venc_list, 1, []);

        % Internal theta = [Vmax, R, e1, e2, B, C1..Cn (, Bt)]
        % (e1,e2) = (AR-1)*[cos(2*alpha), sin(2*alpha)]
        [poly_lb, poly_ub] = polyBounds(polyOrder);
        if hasTissue
            Bt_init = mean(abs(s_n(~bl, kHigh)));
            lb  = [0, 1e-6, -1, -1, 0, poly_lb, 0  ];
            ub  = [inf, 2*R0, 1, 1, inf, poly_ub, inf];
        else
            Bt_init = [];
            lb  = [0, 1e-6, -1, -1, 0, poly_lb];
            ub  = [inf, 2*R0, 1, 1, inf, poly_ub];
        end

        if nargin >= 11 && ~isempty(theta0_phys)
            % theta0_phys in physical units: [Vmax,R,AR,alpha, B,C1..Cn (,Bt)]
            th = theta0_phys(:)';
            [e1_0, e2_0] = ar_alpha_to_e(th(3), th(4));
            theta0 = [th(1), th(2), e1_0, e2_0, th(5:end) / mScale];
            theta0 = min(max(theta0, lb), ub);
        else
            [velFit0, ~] = fitVelProfile(r(bl), v_init, R);
            theta0 = [velFit0.Vmax, R0, 0, 0, B, zeros(1,polyOrder)];
            if hasTissue; theta0 = [theta0, Bt_init]; end
        end

        theta = lsqnonlin( ...
            @(th) residuals_complex_ellipse(th, r, p, s_n, venc_row, s_std, bl, hasTissue), ...
            theta0, lb, ub, opts);

        [AR_f, alpha_f] = e_to_ar_alpha(theta(3), theta(4));
        velFit   = sfit(ft_vel_ellipse, theta(1),theta(2),AR_f,alpha_f, 0, 0);
        velFit1d = cfit(ft_vel, theta(1), theta(2));
        c_start  = 5;
        C_vals   = num2cell(theta(c_start : c_start+polyOrder) * mScale);
        magFit   = cfit(ft_mag, C_vals{:});

        if nargout >= 4
            mag_names = [{'B'}, arrayfun(@(i)sprintf('C%d',i),1:polyOrder,'UniformOutput',false)];
            [AR_0, alpha_0] = e_to_ar_alpha(theta0(3), theta0(4));
            if hasTissue
                sc = [1,1,1,1, repmat(mScale,1,polyOrder+2)];
                fitInfo.names  = [{'Vmax','R','e1','e2'}, mag_names, {'Bt','FEoffset','PEoffset'}];
                fitInfo.units  = [{'cm/s','mm','-','-'}, repmat({'a.u.'},1,polyOrder+2), {'mm','mm'}];
                fitInfo.fixed  = [false(1, 6+polyOrder), true, true];
                fitInfo.lb     = [0, 1e-6, -1, -1, 0, poly_lb, 0,   0, 0];
                fitInfo.ub     = [inf, 2*R0, 1, 1, inf, poly_ub, inf, 0, 0];
            else
                sc = [1,1,1,1, repmat(mScale,1,polyOrder+1)];
                fitInfo.names  = [{'Vmax','R','e1','e2'}, mag_names, {'FEoffset','PEoffset'}];
                fitInfo.units  = [{'cm/s','mm','-','-'}, repmat({'a.u.'},1,polyOrder+1), {'mm','mm'}];
                fitInfo.fixed  = [false(1, 5+polyOrder), true, true];
                fitInfo.lb     = [0, 1e-6, -1, -1, 0, poly_lb, 0, 0];
                fitInfo.ub     = [inf, 2*R0, 1, 1, inf, poly_ub, 0, 0];
            end
            fitInfo.theta0 = [theta0 .* sc, 0, 0];
            fitInfo.theta  = [theta  .* sc, 0, 0];
            fitInfo.derived.names  = {'AR', 'alpha'};
            fitInfo.derived.units  = {'-', 'rad'};
            fitInfo.derived.theta0 = [AR_0, alpha_0];
            fitInfo.derived.theta  = [AR_f, alpha_f];
            if hasTissue; fitInfo.Bt = theta(end) * mScale; end
        end

    otherwise
        error('fitMagVelProfile: unknown fitType ''%s''. Choose: sequential | joint | complex.', fitType);

end
end


function res = residuals_joint(theta, r, v, m, mNoflow, sv, sm)
Vmax = theta(1);  R = theta(2);
C    = num2cell(theta(3:end));
v_pred = velocity_func_radius(r, Vmax, R);
res_v  = (v - v_pred)                                / sv;
res_m  = (m - magnitude_func_velocity(v_pred, C{:})) / sm;
if isempty(mNoflow)
    res_mNF = [];
else
    res_mNF = (mNoflow - magnitude_func_velocity(0, C{:})) / sm;
end
res = [res_v(:); res_m(:); res_mNF(:)];
end


function res = residuals_joint_off(theta, r, p, v, m, mNoflow, sv, sm)
Vmax = theta(1);  R = theta(2);
C    = num2cell(theta(3 : end-2));
FEoffset = theta(end-1);  PEoffset = theta(end);
% p=0 → +PE: r*cos(p)=PE component, r*sin(p)=-FE component
r_off  = sqrt((r.*cos(p) - PEoffset).^2 + (-r.*sin(p) - FEoffset).^2);
v_pred = velocity_func_radius(r_off, Vmax, R);
res_v  = (v - v_pred)                                / sv;
res_m  = (m - magnitude_func_velocity(v_pred, C{:})) / sm;
if isempty(mNoflow)
    res_mNF = [];
else
    res_mNF = (mNoflow - magnitude_func_velocity(0, C{:})) / sm;
end
res = [res_v(:); res_m(:); res_mNF(:)];
end


function res = residuals_joint_ellipse(theta, r, p, v, m, mNoflow, sv, sm)
% theta = [Vmax, R, e1, e2, B, C1..Cn, FEoffset, PEoffset]  (e1,e2 = Cartesian ellipse)
Vmax=theta(1); R=theta(2);
[AR, alpha] = e_to_ar_alpha(theta(3), theta(4));
C = num2cell(theta(5:end-2));
FEoffset=theta(end-1); PEoffset=theta(end);
v_pred = velocity_func_ellipse(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset);
res_v  = (v - v_pred)                                / sv;
res_m  = (m - magnitude_func_velocity(v_pred, C{:})) / sm;
if isempty(mNoflow)
    res_mNF = [];
else
    res_mNF = (mNoflow - magnitude_func_velocity(0, C{:})) / sm;
end
% Penalty: m'(Vmax) = C1 + 2*C2*Vmax >= 0 (nonlinear constraint; box bounds enforce C1>=0, C2<=0)
if numel(C) >= 3
    res_mono = max(0, -(C{2} + 2*C{3}*Vmax)) / sm;
else
    res_mono = [];
end
res = [res_v(:); res_m(:); res_mNF(:); res_mono];
end


% velocity_func_ellipse is a standalone file: velocity_func_ellipse.m


function res = residuals_joint_asym(theta, r, p, v, m, mNoflow, sv, sm)
% theta = [Vmax, R, e1, e2, B, C1..Cn, FEoffset, PEoffset, eFE, ePE]
Vmax = theta(1); R = theta(2);
[AR, alpha] = e_to_ar_alpha(theta(3), theta(4));
C    = num2cell(theta(5 : end-4));
FEoffset = theta(end-3); PEoffset = theta(end-2);
eFE      = theta(end-1); ePE      = theta(end);
v_pred = velocity_func_ellipse_asym(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset, eFE, ePE);
res_v  = (v - v_pred)                                / sv;
res_m  = (m - magnitude_func_velocity(v_pred, C{:})) / sm;
if isempty(mNoflow)
    res_mNF = [];
else
    res_mNF = (mNoflow - magnitude_func_velocity(0, C{:})) / sm;
end
res = [res_v(:); res_m(:); res_mNF(:)];
end


function v = velocity_func_ellipse_asym(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset, eFE, ePE)
% Parabolic velocity profile with elliptical wall and offset velocity peak.
%   Velocity peak at (FEoffset, PEoffset); wall ellipse centred at (FEoffset+eFE, PEoffset+ePE).
%   Semi-major axis R along angle alpha from PE axis; semi-minor axis R/AR (AR>=1).
%   For each spin, R_eff is the ray distance from the velocity peak to the ellipse wall.
%   v = Vmax * max(0, 1 - (r_v/R_eff)^2)
%
%   Derivation: in the ellipse-aligned frame (rotated by alpha), a ray from the velocity peak
%   in direction (uPE,uFE) intersects the ellipse at t satisfying:
%     (t*A - C1)^2/R^2 + AR^2*(t*B - C2)^2/R^2 = 1
%   where A,B are rotated direction components and C1,C2 are rotated ellipse centre offset.
%   Solving the quadratic gives t = R_eff.

% Spin position relative to velocity peak (in PE,FE coordinates)
dPE = r.*cos(p)  - PEoffset;
dFE = -r.*sin(p) - FEoffset;
r_v = sqrt(dPE.^2 + dFE.^2);

% Unit direction from velocity peak toward spin
uPE = dPE ./ max(r_v, eps);
uFE = dFE ./ max(r_v, eps);

% Rotate into ellipse-aligned frame (major axis along alpha from PE axis)
A  = uPE.*cos(alpha) + uFE.*sin(alpha);   % along major axis
B  = -uPE.*sin(alpha) + uFE.*cos(alpha);  % along minor axis
C1 = ePE.*cos(alpha)  + eFE.*sin(alpha);  % ellipse offset along major axis
C2 = -ePE.*sin(alpha) + eFE.*cos(alpha);  % ellipse offset along minor axis

% Quadratic coefficients: (t*A - C1)^2 + AR^2*(t*B - C2)^2 = R^2
a_c  = A.^2 + AR.^2.*B.^2;
b_c  = -2.*(A.*C1 + AR.^2.*B.*C2);
c_c  = C1.^2 + AR.^2.*C2.^2 - R.^2;
disc = max(b_c.^2 - 4.*a_c.*c_c, 0);
R_eff = (-b_c + sqrt(disc)) ./ (2.*a_c);

v = Vmax .* max(0, 1 - (r_v ./ max(R_eff, eps)).^2);
v(r_v < eps) = Vmax;  % at peak: v = Vmax
end


function res = residuals_complex_ellipse(theta, r, p, s_n, venc_row, s_std, bl, hasTissue)
% theta = [Vmax, R, e1, e2, B, C1..Cn]            (all-blood mode)
%       = [Vmax, R, e1, e2, B, C1..Cn, Bt]        (blood+tissue mode)
% s_n      — (N×K) complex, normalised (all pixels: blood first if hasTissue)
% bl       — logical N×1, true=blood  (pass true(N,1) for all-blood mode)
% hasTissue— scalar logical
if nargin < 8; hasTissue = false; end
if nargin < 7 || isempty(bl); bl = true(size(r)); end

if hasTissue
    Bt    = theta(end);
    theta = theta(1:end-1);
end
Vmax=theta(1); R_v=theta(2);
[AR, alpha] = e_to_ar_alpha(theta(3), theta(4));
C = num2cell(theta(5:end));

% Blood pixels: parabolic velocity + inflow polynomial
v_pred  = velocity_func_ellipse(r(bl), p(bl), Vmax, R_v, AR, alpha, 0, 0);  % Nb×1
m_pred  = magnitude_func_velocity(v_pred, C{:});                               % Nb×1
s_pred_b = m_pred .* exp(1i .* pi .* v_pred ./ venc_row);                     % Nb×K

% Tissue pixels: static (v=0), predicted signal = Bt (constant across VENCs)
if hasTissue
    Nt       = sum(~bl);
    K        = size(venc_row, 2);
    s_pred_t = repmat(Bt, Nt, K);   % phase = 0 → complex signal is real = Bt
    diff = zeros(size(s_n), 'like', 1+1i);
    diff(bl,  :) = s_n(bl,  :) - s_pred_b;
    diff(~bl, :) = s_n(~bl, :) - s_pred_t;
else
    diff = s_n - s_pred_b;
end
res = [real(diff(:)); imag(diff(:))] / s_std;
% Penalty: m'(Vmax) = C1 + 2*C2*Vmax >= 0 (nonlinear constraint; box bounds enforce C1>=0, C2<=0)
if numel(C) >= 3
    res = [res; max(0, -(C{2} + 2*C{3}*Vmax)) / s_std];
end
end


function [e1, e2] = ar_alpha_to_e(AR, alpha)
% (AR, alpha) → Cartesian ellipse coords: e = (AR-1)*[cos(2α), sin(2α)]
d  = AR - 1;
e1 = d * cos(2*alpha);
e2 = d * sin(2*alpha);
end

function [AR, alpha] = e_to_ar_alpha(e1, e2)
% Cartesian ellipse coords → (AR, alpha)
d     = sqrt(e1^2 + e2^2);
AR    = 1 + d;
alpha = atan2(e2, e1) / 2;
end

function [lb, ub] = polyBounds(polyOrder)
% Box bounds for [C1..Cn]: C1 >= 0 (increasing at v=0), Cn <= 0 (decelerating).
if polyOrder >= 2
    lb = [0, -inf(1, polyOrder-1)];
    ub = [inf(1, polyOrder-1), 0];
else
    lb = zeros(1, polyOrder);   % C1 >= 0 only for linear polynomial
    ub = inf(1, polyOrder);
end
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
