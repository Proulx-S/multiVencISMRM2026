function [vMean,phi] = CDphase2vel(CDphase, venc)
    % CDphase: ∠CD phase of bulk signal complex difference (CD) [rad]
    % venc   : velocity encoding values [cm/s]
    % vMean  : mean velocity [cm/s]
    % phi    : 𝜙 phase of flowing-spin component vector of bulk signal [rad]
    %          a.k.a. corrected phase difference


    % Forward model:
    % ∠CD =  𝜋/2 - 𝜙/2; for ∠CD > 0
    % ∠CD = -𝜋/2 - 𝜙/2; for ∠CD < 0
    % 𝜙   = vMean/venc*𝜋;

    % Inverse model:
    phi = nan(size(CDphase));
    phi(CDphase>=0) = -pi+2*CDphase(CDphase>=0);
    phi(CDphase<0 ) =  pi+2*CDphase(CDphase< 0);
    vMean = phi./pi.*venc;




    % % Forward model:
    % % ∠CD = 𝜋/2 + 𝜙/2;
    % % 𝜙   = vMean/venc*𝜋;

    % % Inverse model:
    % phi = wrapToPi( 2.*CDphase-pi );
    % vMean = phi./pi.*venc;
end