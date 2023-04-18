%% Find Fermi level relative to a reference energy
% Input:
% q, k, T - electronic charge, Boltzmann const, temperature
% energyLevels - Energy levels found for system
% DOS - DOS found for system
% rel_occ_states - Total states to populate in system
% Output:
% relEf - Fermi level at which the levels occupied sum to the given,
% relative to reference energy given
function relEf = findDeltaEf(q, k, T, energyLevels, DOS, rel_occ_states, E_reference)
% Threshold
threshold = 1e-12;

% Return reference E if sum of states is less than threshold
if abs(rel_occ_states) <= threshold
    relEf = 0;
    return
end

% Determine interval to seek in
minEnergy = min(energyLevels) - 100 * k * T / q;     % eV + J/K * K / C/e = eV + V*C / C/e = eV
maxEnergy = max(energyLevels) + 100 * k * T / q;     % eV
if rel_occ_states < 0
    E_interval = [E_reference, minEnergy];
else
    E_interval = [E_reference, maxEnergy];
end

sum_occupied = 0;
Ef_guess = 0;
while abs(sum_occupied - rel_occ_states) > threshold
    Ef_guess = mean(E_interval);
    % Find sum of occupied states from reference to guess
    sum_occupied = sumStatesInInterval(q, k , T, energyLevels, DOS, [E_reference, Ef_guess]);
    if abs(sum_occupied) < abs(rel_occ_states)
        E_interval = [Ef_guess, E_interval(2)];
    else
        E_interval = [E_interval(1), Ef_guess];
    end
end

relEf = Ef_guess - E_reference;

end