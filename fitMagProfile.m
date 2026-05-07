function [magFit, magFitFixR, magFitFixB, magFitFixBR] = fitMagProfile(r, m, R, B, polyOrder)
% Polynomial m(r) magnitude function of radial position in a cylinder, with dm/dr=0 at r=0.
%   polyOrder   — polynomial degree (default 2; supports 2–4; no linear term)
%   magFit      — B free,  R free:  {B, C2[..CN]} fitted
%   magFitFixR  — B free,  R fixed: {B, C2[..CN]} fitted, m=B at r=R enforced
%   magFitFixB  — B fixed, R free:  {C2[..CN], R} fitted, m=B at r=R enforced
%   magFitFixBR — B fixed, R fixed: {C2[..CN]}    fitted, m=B at r=R enforced
%   Ck ≤ 0 for k=2..N (magnitude non-increasing from centre to wall)
%   B           — wall magnitude; start point or fixed value
%   R           — wall radius [mm]; start point or fixed value
r = double(r(:)); m = double(m(:));
if ~exist('polyOrder','var') || isempty(polyOrder); polyOrder = 2; end
if ~ismember(polyOrder, 2:4)
    error('fitMagProfile: polyOrder must be 2, 3, or 4.');
end
if exist('B','var') && ~isempty(B);  B = double(B);  end
if exist('R','var') && ~isempty(R);  R = double(R);  end

haveB = exist('B','var') && ~isempty(B);
haveR = exist('R','var') && ~isempty(R);

B_start  = haveB * B          + (~haveB) * min(m);
R_start  = haveR * R          + (~haveR) * max(r) * 0.9;
C2_start = (B_start - max(m)) / R_start^2;
sp_C     = [C2_start, zeros(1, polyOrder-2)];
ub_C     = [0, inf(1, polyOrder-2)];

switch polyOrder
    case 2
        ft_free  = fittype('B + C2*r^2',                                         'independent','r',                     'coefficients',{'B','C2'      });
        ft_fixR  = fittype('B + C2*(r^2-R^2)',                                   'independent','r', 'problem',{'R'},    'coefficients',{'B','C2'      });
        ft_fixB  = fittype('B + C2*(r^2-R^2)',                                   'independent','r', 'problem',{'B'},    'coefficients',{'C2','R'      });
        ft_fixBR = fittype('B + C2*(r^2-R^2)',                                   'independent','r', 'problem',{'B','R'},'coefficients',{'C2'          });
    case 3
        ft_free  = fittype('B + C2*r^2 + C3*r^3',                               'independent','r',                     'coefficients',{'B','C2','C3'      });
        ft_fixR  = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3)',                   'independent','r', 'problem',{'R'},    'coefficients',{'B','C2','C3'      });
        ft_fixB  = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3)',                   'independent','r', 'problem',{'B'},    'coefficients',{'C2','C3','R'      });
        ft_fixBR = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3)',                   'independent','r', 'problem',{'B','R'},'coefficients',{'C2','C3'          });
    case 4
        ft_free  = fittype('B + C2*r^2 + C3*r^3 + C4*r^4',                     'independent','r',                     'coefficients',{'B','C2','C3','C4'      });
        ft_fixR  = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3) + C4*(r^4-R^4)',   'independent','r', 'problem',{'R'},    'coefficients',{'B','C2','C3','C4'      });
        ft_fixB  = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3) + C4*(r^4-R^4)',   'independent','r', 'problem',{'B'},    'coefficients',{'C2','C3','C4','R'      });
        ft_fixBR = fittype('B + C2*(r^2-R^2) + C3*(r^3-R^3) + C4*(r^4-R^4)',   'independent','r', 'problem',{'B','R'},'coefficients',{'C2','C3','C4'          });
end

% B free, R free: {B, C2[..CN]} fitted
magFit = fit(r, m, ft_free, 'StartPoint', [max(m), sp_C], 'Upper', [inf, ub_C]);

% B free, R fixed: m=B at r=R enforced; {B, C2[..CN]} fitted
if haveR
    magFitFixR = fit(r, m, ft_fixR, 'problem', {R}, 'StartPoint', [B_start, sp_C], 'Upper', [inf, ub_C]);
else
    magFitFixR = [];
end

% B fixed, R free: m=B at r=R enforced; {C2[..CN], R} fitted
if haveB
    magFitFixB = fit(r, m, ft_fixB, 'problem', {B}, 'StartPoint', [sp_C, R_start], 'Upper', [ub_C, inf]);
else
    magFitFixB = [];
end

% B fixed, R fixed: m=B at r=R enforced; {C2[..CN]} fitted
if haveB && haveR
    magFitFixBR = fit(r, m, ft_fixBR, 'problem', {B, R}, 'StartPoint', sp_C, 'Upper', ub_C);
else
    magFitFixBR = [];
end
