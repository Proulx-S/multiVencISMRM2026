function [velFit, velFitFixR] = fitVelProfile(r, v, R)
% Parabolic fits of velocity v at radial position r in cylinder of radius R with maximum velocity Vmax.
%   velFit      — free fit: both Vmax and R fitted
%   velFitFixR  — constrained fit: R fixed at input R, Vmax free
%   R           — cylinder radius start point (free fit) or fixed value (constrained fit)

r = double(r(:)); v = double(v(:));

% Free fit
velFit = fittype('Vmax * (1 - (r/R)^2)', 'independent', 'r', 'coefficients', {'Vmax', 'R'});
if exist('R','var') && ~isempty(R)
    % use R as starting point
    velFit = fit(r, v, velFit, 'Lower', [0, 0], 'StartPoint', [max(v), R]);
else
    % R estimated from the size of the ROI
    velFit = fit(r, v, velFit, 'Lower', [0, 0], 'StartPoint', [max(v), sqrt(max(r)^2/2)/2]);
end

% Constrained radius
if exist('R','var') && ~isempty(R)
    velFitFixR  = fittype('Vmax * (1 - (r/R)^2)', 'independent', 'r', 'problem', 'R', 'coefficients', {'Vmax'});
    velFitFixR  = fit(r, v, velFitFixR, 'problem', {R}, 'Lower', 0, 'StartPoint', max(v));
else
    velFitFixR = [];
end

