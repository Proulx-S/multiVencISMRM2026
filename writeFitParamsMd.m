function writeFitParamsMd(filepath, fitLabel, fitInfo, costDesc, dataSummary)
% Write a fit-parameter report as a Markdown file.
%
%   filepath    — output path (e.g. '.../sub-01_vessel-01_fitParams_A.md')
%   fitLabel    — short label for the H1 title (e.g. 'A — Joint real residual')
%   fitInfo     — struct from fitMagVelProfile (names, units, fixed, lb, theta0, theta)
%   costDesc    — struct with fields:
%                   .equation   — LaTeX string for the residual vector
%                   .predicted  — LaTeX string for the predicted quantity (or '')
%                   .normDesc   — plain-text description of the normalisation
%                   .dataDesc   — plain-text description of input data
%   dataSummary — struct with fields:
%                   .N_pixels, .K_vencs, .vencs (vector, cm/s), .bestVenc (cm/s)

fid = fopen(filepath, 'w');

% -------------------------------------------------------------------------
% Title
% -------------------------------------------------------------------------
fprintf(fid, '# Fit %s\n\n', fitLabel);

% -------------------------------------------------------------------------
% Parameter table — in a fenced code block so it is monospace in both
% Markdown preview and plain-text viewers.
% -------------------------------------------------------------------------
fprintf(fid, '## Parameters\n\n');
fprintf(fid, '`*` marks parameters fixed at their listed value (not optimised).\n\n');
fprintf(fid, '```\n');

% Header
hdr = sprintf('%-16s  %-10s  %-7s  %12s  %12s  %12s  %12s', ...
    'Parameter', 'Units', 'Fixed?', 'LB', 'UB', 'Initial', 'Final');
fprintf(fid, '%s\n', hdr);
fprintf(fid, '%s\n', repmat('-', 1, numel(hdr)));

% Rows
for ii = 1:numel(fitInfo.names)
    name  = fitInfo.names{ii};
    unit  = fitInfo.units{ii};
    fixed = fitInfo.fixed(ii);
    if fixed; name = [name ' *']; end  %#ok<AGROW>
    fprintf(fid, '%-16s  %-10s  %-7s  %12s  %12s  %12s  %12s\n', ...
        name, unit, yesno(fixed), ...
        fmtNum(fitInfo.lb(ii)), fmtNum(fitInfo.ub(ii)), ...
        fmtNum(fitInfo.theta0(ii)), fmtNum(fitInfo.theta(ii)));
end
fprintf(fid, '```\n\n');

% Derived parameters (e.g. AR, alpha from Cartesian e1/e2)
if isfield(fitInfo, 'derived') && ~isempty(fitInfo.derived)
    d = fitInfo.derived;
    fprintf(fid, 'Derived:\n\n');
    fprintf(fid, '```\n');
    dhdr = sprintf('%-16s  %-10s  %12s  %12s', 'Parameter', 'Units', 'Initial', 'Final');
    fprintf(fid, '%s\n', dhdr);
    fprintf(fid, '%s\n', repmat('-', 1, numel(dhdr)));
    for ii = 1:numel(d.names)
        fprintf(fid, '%-16s  %-10s  %12s  %12s\n', ...
            d.names{ii}, d.units{ii}, fmtNum(d.theta0(ii)), fmtNum(d.theta(ii)));
    end
    fprintf(fid, '```\n\n');
end

% -------------------------------------------------------------------------
% Initial value sources
% -------------------------------------------------------------------------
if isfield(fitInfo, 'init_notes') && ~isempty(fitInfo.init_notes)
    fprintf(fid, '### Initial value sources\n\n');
    fprintf(fid, '| Parameter | Source |\n');
    fprintf(fid, '|---|---|\n');
    for ii = 1:numel(fitInfo.names)
        name = fitInfo.names{ii};
        if fitInfo.fixed(ii); name = [name ' *']; end  %#ok<AGROW>
        fprintf(fid, '| %s | %s |\n', name, fitInfo.init_notes{ii});
    end
    fprintf(fid, '\n');
end

fprintf(fid, '---\n\n');

% -------------------------------------------------------------------------
% Model
% -------------------------------------------------------------------------
fprintf(fid, '## Model\n\n');
fprintf(fid, '### Velocity profile\n\n');
fprintf(fid, 'Parabolic profile with elliptical wall; ');
fprintf(fid, 'peak at $(\\text{FEoffset}, \\text{PEoffset}) = (0, 0)$:\n\n');
fprintf(fid, '$$v_i = V_{\\max} \\cdot \\max\\!\\left(0,\\ 1 - \\left(\\frac{r_{v,i}}{R_{\\text{eff},i}}\\right)^2\\right)$$\n\n');
fprintf(fid, 'where $r_{v,i}$ is the distance from the velocity peak, and the direction-dependent wall radius is\n\n');
fprintf(fid, '$$R_{\\text{eff},i} = \\frac{R}{\\sqrt{A_i^2 + AR^2\\, B_i^2}}$$\n\n');
fprintf(fid, 'with $(A_i, B_i)$ the unit direction projected onto the ellipse axes ');
fprintf(fid, '(semi-major $R$ at angle $\\alpha$ from PE axis, semi-minor $R/AR$, $AR \\ge 1$).\n\n');
fprintf(fid, '### Magnitude–velocity relationship\n\n');
fprintf(fid, '$$m(v) = B + C_1 v + C_2 v^2$$\n\n');
fprintf(fid, '---\n\n');

% -------------------------------------------------------------------------
% Cost function
% -------------------------------------------------------------------------
fprintf(fid, '## Cost function\n\n');
fprintf(fid, 'Minimised by `lsqnonlin` (trust-region reflective).\n\n');
fprintf(fid, '$$%s$$\n\n', costDesc.equation);
if ~isempty(costDesc.predicted)
    fprintf(fid, 'where the predicted signal is\n\n');
    fprintf(fid, '$$%s$$\n\n', costDesc.predicted);
end
fprintf(fid, '%s\n\n', costDesc.normDesc);
fprintf(fid, '**Data**: %s\n\n', dataSummary);

fclose(fid);
end


function s = fmtNum(x)
if isinf(x) && x > 0
    s = '         Inf';
elseif isinf(x) && x < 0
    s = '        -Inf';
else
    s = sprintf('%12.4g', x);
end
end

function s = yesno(tf)
if tf; s = 'yes'; else; s = 'no'; end
end
