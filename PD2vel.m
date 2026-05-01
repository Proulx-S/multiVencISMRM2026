function vel = PD2vel(theta, venc)
    % theta: phase difference (PD) [rad]
    % venc: velocity encoding [cm/s]
    % vel: velocity [cm/s]

    % Model: theta = pi * vel / venc
    
    % Inverse model:
    vel = theta .* venc./pi;
end