function [magFit, magFitFixR, magFitFixB, magFitFixBR] = fitMagProfile(r, m, R, B)
% Order-2 polynomial m(r) magnitude function of radial position in a cylinder, with dm/dr=0 at r=0.
%   magFit      — B free,  R free:  {B, C2} fitted
%   magFitFixR  — B free,  R fixed: {B, C2} fitted, m=B at r=R enforced
%   magFitFixB  — B fixed, R free:  {C2, R} fitted, m=B at r=R enforced
%   magFitFixBR — B fixed, R fixed: {C2}    fitted, m=B at r=R enforced
%   C2 constrained to ≤0 in all fits (magnitude decreases from center to wall)
%   B           — wall magnitude; start point or fixed value
%   R           — wall radius [mm]; start point or fixed value
r = double(r(:)); m = double(m(:));
if exist('B','var') && ~isempty(B);  B = double(B);  end
if exist('R','var') && ~isempty(R);  R = double(R);  end

haveB = exist('B','var') && ~isempty(B);
haveR = exist('R','var') && ~isempty(R);

B_start  = haveB * B          + (~haveB) * min(m);
R_start  = haveR * R          + (~haveR) * max(r) * 0.9;
C2_start = (B_start - max(m)) / R_start^2;

% B free, R free: {B, C2} fitted
magFit = fittype('B + C2*r^2', 'independent', 'r', 'coefficients', {'B', 'C2'});
magFit = fit(r, m, magFit, 'StartPoint', [max(m), C2_start], 'Upper', [inf, 0]);

% B free, R fixed: m=B at r=R enforced; {B, C2} fitted
if haveR
    ftFixR     = fittype('B + C2*(r^2-R^2)', 'independent', 'r', 'problem', {'R'}, 'coefficients', {'B', 'C2'});
    magFitFixR = fit(r, m, ftFixR, 'problem', {R}, 'StartPoint', [B_start, C2_start], 'Upper', [inf, 0]);
else
    magFitFixR = [];
end

% B fixed, R free: m=B at r=R enforced; {C2, R} fitted
if haveB
    ftFixB     = fittype('B + C2*(r^2-R^2)', 'independent', 'r', 'problem', {'B'}, 'coefficients', {'C2', 'R'});
    magFitFixB = fit(r, m, ftFixB, 'problem', {B}, 'StartPoint', [C2_start, R_start], 'Upper', [0, inf]);
else
    magFitFixB = [];
end

% B fixed, R fixed: m=B at r=R enforced; {C2} fitted
if haveB && haveR
    ftFixBR     = fittype('B + C2*(r^2-R^2)', 'independent', 'r', 'problem', {'B', 'R'}, 'coefficients', {'C2'});
    magFitFixBR = fit(r, m, ftFixBR, 'problem', {B, R}, 'StartPoint', C2_start, 'Upper', 0);
else
    magFitFixBR = [];
end
