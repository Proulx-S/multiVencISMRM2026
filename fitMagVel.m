function [magFit, magFitFixR, magFitFixB, magFitFixBR] = fitMagVel(m, v, B)
% Order-2 polynomial m(v) magnitude function of velocity in blood.
%   magFit      — B free,  Vref free:  {B, C1, C2} fitted, flow data only
%   magFitFixR  — B free,  Vref=0 fixed: {B, C1, C2} fitted, flow + (0,B) anchor point
%   magFitFixB  — B fixed, Vref free:  {C1, C2, Vref} fitted, m(Vref)=B enforced
%   magFitFixBR — B fixed, Vref=0 fixed: {C1, C2} fitted, m(0)=B enforced
%   B           — m(0), no-flow magnitude; start point (free fits) or fixed value (constrained fits)
v = double(v(:)); m = double(m(:));
if exist('B','var') && ~isempty(B);  B = double(B);  end

haveB   = exist('B','var') && ~isempty(B);
B_start = haveB * B + (~haveB) * min(m);

% B free, Vref free: {B, C1, C2} fitted
magFit = fittype('B + C1*v + C2*v^2', 'independent', 'v', 'coefficients', {'B', 'C1', 'C2'});
magFit = fit(v, m, magFit, 'StartPoint', [B_start, 0, 0]);

% B free, Vref=0 fixed: {B, C1, C2} fitted, (0,B) anchor point added to flow data
if haveB
    ftFixR     = fittype('B + C1*v + C2*v^2', 'independent', 'v', 'coefficients', {'B', 'C1', 'C2'});
    magFitFixR = fit([v; 0], [m; B], ftFixR, 'StartPoint', [B_start, 0, 0]);
else
    magFitFixR = [];
end

% B fixed, Vref free: m(Vref)=B enforced; {C1, C2, Vref} fitted
if haveB
    ftFixB     = fittype('B + C1*(v-Vref) + C2*(v^2-Vref^2)', ...
        'independent', 'v', 'problem', {'B'}, 'coefficients', {'C1', 'C2', 'Vref'});
    magFitFixB = fit(v, m, ftFixB, 'problem', {B}, 'StartPoint', [0, 0, max(v)]);
else
    magFitFixB = [];
end

% B fixed, Vref=0 fixed: {C1, C2} fitted, m(0)=B enforced
if haveB
    ftFixBR     = fittype('B + C1*v + C2*v^2', ...
        'independent', 'v', 'problem', {'B'}, 'coefficients', {'C1', 'C2'});
    magFitFixBR = fit(v, m, ftFixBR, 'problem', {B}, 'StartPoint', [0, 0]);
else
    magFitFixBR = [];
end
