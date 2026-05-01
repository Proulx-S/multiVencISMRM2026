function f = getffromI0andIf(I0, If, Ss)
    % I0: measured signal under flow compensation
    % If: magnitude of signal contributed by the flowing spin compartment
    % Ss: magnitude of signal contributed by the static spin compartment
    %  f: volume fraction of flowing spins, assuming single flowing and static compartments

    % Forward model:
    % I0=(1-f)Ss+If

    f = ( If-abs(I0) )/ Ss + 1; % taking the absolute of I0 here assumes there is some noise-related phase variations around 0. If we take the real value, noisy phase would lead to underestimation of measured magnitude.
end