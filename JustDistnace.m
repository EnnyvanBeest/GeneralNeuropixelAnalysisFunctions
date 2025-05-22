Shank0_Probe0 = [-1650 680]; % Probe 0, Shank 0 position
Shank0_Probe1Prov = [-2480 -2460]; % Probe 1, Shank 0 estimated position

dist = sqrt(sum((Shank0_Probe0 - Shank0_Probe1Prov).^2))
    