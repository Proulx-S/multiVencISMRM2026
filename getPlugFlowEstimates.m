function [velCD,phi,velPD,If,IfSign,f] = getPlugFlowEstimates(venc,CD,I0,Ss,PD,velCDlamiCorFac,velPDlamiCorFac,unwrapVec)
    % venc: velocity encoding [cm/s]
    %   CD: complex difference signal
    %   Ss: magnitude of signal from individual static spins
    %   I0: measured signal under flow compensation
    %   PD: phase difference signal [rad]
    % velCDlamiCorFac: correction factor for the CD-based velocity and phi estimate to account for laminar flow bias
    if ~exist('velCDlamiCorFac','var') || isempty(velCDlamiCorFac)
        velCDlamiCorFac = 0.75;
    end
    if ~exist('velPDlamiCorFac','var') || isempty(velPDlamiCorFac)
        velPDlamiCorFac = 1;
    end
    if ~exist('unwrapVec','var') || isempty(unwrapVec)
        unwrapVec = zeros(size(venc));
    end

    %  velCD: CD-based estimate of mean velocity [cm/s]
    %    phi: estimate of the phase of the bulk signal from the flowing spins compartment
    %  velPD: PD-based estimate of mean velocity [cm/s]
    %     If: estimate of the magnitude of signal contributed by the flowing spin compartment (before any velocity phase dispersion)
    % IfSign: sign of the estimate of If
    %      f: estimate of the volume fraction of flowing spins, assuming single flowing and static compartments


    % Velocity estimate from CD [equation 4]
    [velCD,phi] = CDphase2vel( angle(CD)+unwrapVec.*2*pi , venc );

    % Bias correction for laminar flow
    velCD = velCD * velCDlamiCorFac;
    phi   = phi   * velCDlamiCorFac;

    % Velocity estimate from PD
    velPD = PD2vel( PD,venc );

    % Bias correction for laminar flow
    velPD = velPD * velPDlamiCorFac;

    % If estimate from CD [equation 5]
    [If,IfSign] = CDmag2fSf( abs(CD),phi );

    % f estimate from I0, If and Ss [equation 6]
    if ~isempty(I0) && ~isempty(If) && ~isempty(Ss)
        f = getffromI0andIf( I0,If,Ss );
    else
        f = [];
    end
end