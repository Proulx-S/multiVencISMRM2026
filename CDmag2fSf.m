function [fSf,fSf_sign] = CDmag2fSf(CDmag, phi)
    % CDmag : magnitude of the complex differente bulk signal (CD) [a.u.]
    % phi   : phase of flowing spin bulk signal [rad]
    % fSf   : magnitude of the flowing spin bulk signal [a.u.]

    % Forward model:
    % CDmag = fSf * 2 * sin(phi/2);

    % Inverse model:
    fSf = CDmag ./ (2 * sin(phi/2));
    fSf_sign = sign(fSf);
    fSf = abs(fSf);