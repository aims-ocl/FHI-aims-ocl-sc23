%% Takes energy levels and weights (DOS) and an energy interval and
% returns the total states associated with the given energy interval using
% the Fermi-Dirac distribution
% Input:
% q, k, T - electronic charge, Boltzmann const, temperature
% Output:
% sumStates - Sum of states within interval
function sumOfStates = sumStatesInInterval(q, k , T, energyLevels, DOS, E_interval)
% Reshape interval to two columns
E_interval = reshape(E_interval, 1, 2);
% Reshape energy levels, DOS in columns
numEnergies = numel(energyLevels);
energyLevels = reshape(energyLevels, numEnergies, 1);
DOS = reshape(DOS, numEnergies, 1);

% Fermi-Dirac distributions for interval
occupation_distribution = 1 ./ (exp((energyLevels * ones(1, 2) - ones(numEnergies, 1) * E_interval) .* q / (k * T)) + 1);     % exp((eV*C/e)/J) - Unitless between 0 and 1

occupied_DOS = (DOS * ones(1, 2)) .* occupation_distribution;
sumOfStates = diff(sum(occupied_DOS, 1));

end