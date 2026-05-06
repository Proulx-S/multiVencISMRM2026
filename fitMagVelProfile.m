function [velFit, magFit] = fitMagVelProfile(r, v, m, mNoFlow, R, fitType, polyOrder)
% Fit v(r) and m(v) to flow and no-flow data.
%   fitType  — 'sequential' (default) or 'joint'
%     'sequential' — fit v(r) then m(v) independently; v(r) uses free parabolic fit
%     'joint'      — single lsqnonlin over [Vmax, R, B, C1...Cn] minimising both
%                    velocity error and magnitude error simultaneously
%   polyOrder — polynomial order for m(v) (default 2; supports 1–4)
%   Models: velocity_func_radius   →  v(r) = Vmax*(1-(r/R)^2)
%           magnitude_func_velocity →  m(v) = B + C1*v + ... + Cn*v^n
%   r, v, m   — radial positions, velocities, magnitudes of flow pixels
%   mNoFlow   — magnitudes of no-flow pixels (v=0); B = mean(mNoFlow)
%   R         — cylinder radius [mm]; used as start point in both cases
r = double(r(:)); v = double(v(:)); m = double(m(:)); mNoFlow = double(mNoFlow(:));
if exist('R'        ,'var') && ~isempty(R);        R        = double(R); end
if ~exist('fitType' ,'var') || isempty(fitType);   fitType  = 'sequential'; end
if ~exist('polyOrder','var') || isempty(polyOrder); polyOrder = 2; end

B = mean(mNoFlow);

% v(r) fittype
ft_vel = fittype(@(Vmax, R, r) velocity_func_radius(r, Vmax, R), ...
    'independent', 'r', 'coefficients', {'Vmax', 'R'});

% m(v) fittype and start point — built dynamically for the requested polynomial order
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
sp_mag = [B, zeros(1, polyOrder)];   % start point: [B, C1=0, ..., Cn=0]

switch fitType

    case 'sequential'
        [velFit_, ~] = fitVelProfile(r, v, R);
        velFit       = cfit(ft_vel, velFit_.Vmax, velFit_.R);
        v_pred       = velocity_func_radius(r, velFit_.Vmax, velFit_.R);
        v_all        = [v_pred(:);  zeros(numel(mNoFlow), 1)];
        m_all        = [m(:);       mNoFlow(:)              ];
        magFit       = fit(v_all, m_all, ft_mag, 'StartPoint', sp_mag);

    case 'joint'
        sv = max(std(v),            eps);
        sm = max(std([m; mNoFlow]), eps);

        [velFit0, ~] = fitVelProfile(r, v, R);
        theta0 = [velFit0.Vmax, velFit0.R, B, zeros(1, polyOrder)];
        lb     = [0, 1e-6, -inf(1, polyOrder + 1)];

        opts  = optimoptions('lsqnonlin', 'Display', 'off');
        theta = lsqnonlin(@(th) residuals_joint(th, r, v, m, mNoFlow, sv, sm), ...
                          theta0, lb, [], opts);

        velFit    = cfit(ft_vel, theta(1), theta(2));
        C_vals    = num2cell(theta(3:end));   % [B, C1, ..., Cn]
        magFit    = cfit(ft_mag, C_vals{:});

    otherwise
        error('fitMagVelProfile: unknown fitType ''%s''. Choose: sequential | joint.', fitType);

end
end


function res = residuals_joint(theta, r, v, m, mNoflow, sv, sm)
Vmax = theta(1);  R = theta(2);
C    = num2cell(theta(3:end));   % {B, C1, ..., Cn}
v_pred  = velocity_func_radius(r, Vmax, R);
res_v   = (v       - v_pred)                                       / sv;
res_m   = (m       - magnitude_func_velocity(v_pred, C{:}))       / sm;
res_mNF = (mNoflow - magnitude_func_velocity(0,      C{:}))       / sm;
res = [res_v(:); res_m(:); res_mNF(:)];
end
