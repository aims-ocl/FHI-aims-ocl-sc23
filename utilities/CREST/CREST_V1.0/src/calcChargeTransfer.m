%% Receives atomistic data for right-hand slab, bulk data for left-hand
% semi-infinite bulk, and finds the amount of charge transfer required to
% obtain simultaneously Fermi-level equilibrium and charge conservation.

function [deltaQ, dphi_b, Ef_out, z_d, U_d, E_s, E_b] = calcChargeTransfer(U, enLev_DOS, Ef0, Qsheet, ...
    useDOS, extrapolateBulklikeU, UAveragingInterval, UblAveragingInterval, ...
    degeneracyLimit_kT, ignoreNumericFail, cc_Qsheet_minimal, cc_dphi_b, cc_Ef, cc_Qsheet, ...
    T, Nv, Nc, epsilon, Eg, Ndop, Wv_FP, eDelta_phi_v_FP, ...
    vacZ, bulklikeZs, z_d_est, surfaceCellArea, Qintrinsic, outputFile)

% LOAD AND UNPACK physical constants
physConsts = initPhysicalConsts;
q = physConsts.q;     % C/e
k = physConsts.k;     % J/K
e = physConsts.e;     % e
epsilon_0 = physConsts.epsilon_0;     % e/(V-cm)


% FIND BULK EF IN CASE OF NO CHARGE TRANSFER
% Find the local vacuum level at the vacuum z
vacU = findYsFromXs(vacZ, U, UAveragingInterval);     % eV

% Try to work out position to stitch bulk to slab at.
z_d = [];
U_d = [];

if abs(Qsheet) <= cc_Qsheet_minimal
    % It doesn't matter what z_d is, E_b is found from the vacuum and
    % the flatband work function
    U_vac_at_z_d = vacU;
    outputMessage(sprintf('Sheet charge negligible; using U in vacuum = %0.6f eV', ...
        vacU), outputFile);
    
else
    % Calculate e*field (U slope) in vacuum from the sheet charge
    % e * e/cm^2 / e/(V-cm) * cm/Ang = e * V/cm * cm/Ang = eV/Ang
    vacField = -e * -Qsheet / epsilon_0 * 1e-8;     % eV/Ang
    % Extrapolate U in the sheet-slab gap. Result is the coefficients
    Uvac_line = getLineFromPAndSlope([vacZ, vacU], vacField);
    outputMessage(sprintf('Electrostatic energy in vacuum extrapolated to : U(z) = %0.5e z + %0.5e', ...
        Uvac_line.a1, Uvac_line.a0), outputFile, 'silent');
    
    % User specified no extrapolation from bulk
    if ~extrapolateBulklikeU
        z_d = z_d_est;
        outputMessage(sprintf('Using fixed z_d = %0.5f %s', ...
            z_d, char(197)), outputFile);
        
    else
        % Correct vacuum line by virtual surface potential drop
        Uvac_corr = Uvac_line;
        Uvac_corr.a0 = Uvac_line.a0 - eDelta_phi_v_FP;
        
        % Find energies at bulk-like positions
        bulklikeUs = findYsFromXs(bulklikeZs, U, UblAveragingInterval);
        % Zip (z,U) pairs
        zUs = [bulklikeZs bulklikeUs];
        % Extrapolate U from bulk-like Zs. Result is the coefficients
        Ubl_parab = getParabolaFrom3P(zUs(1, :), zUs(2, :), zUs(3, :));
        outputMessage(sprintf('Bulk-like potential extrapolated to : U(z) = %0.5e z^2 + %0.5e z + %0.5e', ...
            Ubl_parab.a2, Ubl_parab.a1, Ubl_parab.a0), outputFile, 'silent');
        
        % Examine the intersections of vacuum and BL extrapolations. This
        % scheme is most reliable where (of course) intersections exist,
        % and the leftmost of them (the only physically meaninful
        % intersection assuming epsilon > 1) obeys:
        % Extremum < left intersection < leftmost bulk-like z.
        % However, the left intersection is still probably the best bet in
        % all cases where the curvature and vacuum slope have the same sign.
        % If they do not have the same sign and the user has not specified
        % an estimate for z_d, a heuristic estimation must be made.
        outputMessage('Searching for intersection of bulk-like and vacuum U curves (z_d,U_d)', outputFile, 'silent');
        [P1, P2] = getParabolaLineIntersection(Ubl_parab, Uvac_corr);
        % Identify intersections
        if isempty(P1)
            left_intersect = Inf;
            right_intersect = Inf;
        elseif isempty(P2)
            left_intersect = P1;
            right_intersect = P1;
        elseif P1(1) < P2(1)
            left_intersect = P1;
            right_intersect = P2;
        else
            left_intersect = P2;
            right_intersect = P1;
        end
        
        % Identify correct candidate z-coordinate in parabola case
        if Ubl_parab.a2 / Uvac_corr.a1 > 0
            parab_z_d_candidate = left_intersect(1);
        else
            parab_z_d_candidate = right_intersect(1);
        end
        
        % Best choice: curvature and vac line have same sign and the
        % candidate is left of the leftmost bulk-z position
        if parab_z_d_candidate < min(bulklikeZs) && Ubl_parab.a2 / Uvac_corr.a1 > 0
            z_d = parab_z_d_candidate;
            outputMessage(sprintf('Found intersection at z_d = %0.5f %s', ...
                z_d, char(197)), outputFile);
        else
            outputMessage(['Couldn''t find a reliable intersection of '...
                'bulk-like and vacuum U curves'], outputFile, 'silent');
        end
        
        % Next best choice: user-estimated z_d
        if isempty(z_d) && ~isempty(z_d_est)
            z_d = z_d_est;
            outputMessage(sprintf('Falling back on user-specified z_d = %0.5f %s', ...
                z_d, char(197)), outputFile);
        end
        
        % So-so choice: a hueristic attempt to determine z_d
        if isempty(z_d)
            outputMessage(['***Warning***: Trying to estimate z_d from bulk U. ' ...
                'This is usually not a reliable estimate!'], outputFile, 'silent');
            % Find the intersection with the line passing through a
            % bulk-like potential point with the "bulk" slope at z_d
            E_d = vacField / epsilon;     % eV/Ang
            % Take the mean of the bulklike points
            blZ = mean(bulklikeZs);
            blU = mean(bulklikeUs);
            % Get line
            Ubulk_line = getLineFromPAndSlope([blZ, blU], E_d);
            % Find intersection with corrected vacuum line
            line_z_d_candidate = (Uvac_corr.a0 - Ubulk_line.a0) / (Ubulk_line.a1 - Uvac_corr.a1);
            
            % Use rightmost of the two candidates that is still left of
            % the leftmost bulk-like z
            z_d_candidates = [parab_z_d_candidate line_z_d_candidate];
            z_d = max(z_d_candidates(z_d_candidates <= min(bulklikeZs)));
            if ~isempty(z_d)
                outputMessage(sprintf('Using rough estimation of z_d = %0.5f %s', ...
                    z_d, char(197)), outputFile);
            end
        end
        
        % Bad: No way to obtain z_d
        if isempty(z_d)
            % Exit with error
            err = MException('MATLAB:CREST_main:BulkExtrapolationFailed', ...
                'bulk extrapolation failed but no estimation for z_d provided');
            throw(err)
        end
    end
    
    % Calculate E_b from z_d and vacuum line
    U_vac_at_z_d = Uvac_line.a1 .* z_d + Uvac_line.a0;
    outputMessage(sprintf('Taking vacuum U at z_d (z_d/%s,U_vac(z_d)/eV) = (%0.5f,%0.6f)', ...
        char(197), z_d, U_vac_at_z_d), outputFile);
end

E_b = U_vac_at_z_d - Wv_FP;     % ==> E_b = outputFP.Ef
outputMessage(sprintf('E_b = %0.6f eV', E_b), outputFile, 'silent');

% Calculate U_d if eDelta_phi_v_FP provided
if ~isempty(eDelta_phi_v_FP)
    U_d = U_vac_at_z_d - eDelta_phi_v_FP;
end


if useDOS
    % GENERAL CASE: USE CALCULATION DOS TO EXPEDITE CONVERGENCE
    % Unpack for convenience
    energyLevels = enLev_DOS(:, 1);
    DOS = enLev_DOS(:, 2);
    
    % FIND SLAB EF IN CASE OF NO CHARGE TRANSFER
    % i.e, Ef calculated from the system's intrinsic non-fixed charge
    % alone (without the reservoir-derived charge), based on the
    % system's DOS
    safelyLowRefEnergy = min(energyLevels) - 100 * k * T / q;     % eV + J/K * K / C/e = eV + V*C / C/e = eV
    E_s = safelyLowRefEnergy + findDeltaEf(q, k, T, ...
        energyLevels, DOS, Qintrinsic / (-e), safelyLowRefEnergy);
    outputMessage(sprintf('E_s = %0.6f eV', E_s), outputFile, 'silent');
    
    
    % FIND NEW CHARGING
    % Feed range into function that provides surface (slab) charging
    [deltaQ, dphi_b, Ef_out] = calcCTBinary(E_s, E_b, ...
        T, Ndop, Nv, Nc, epsilon, Eg, ...
        energyLevels, DOS, surfaceCellArea, ...     % cm^2
        cc_Qsheet_minimal, cc_dphi_b, cc_Ef, cc_Qsheet, ...     % e/cm^2,  e/cm^2, eV
        degeneracyLimit_kT, ignoreNumericFail, ...
        outputFile);
    
else
    % GET CHARGE UNDER ASSUMPTION OF HIGH DENSITY OF SLAB STATES
    % Use output Ef as E_s, do not search for equilibrium value
    
    % Find B.B. potential drop Ef0 corresponds to. Convention is that
    % a negative potential drop bends bands UP.
    dphi_b = -(E_b - Ef0) / e;     % eV / e = V
    
    % Find surface charge Qs that matches the band bending potential
    % drop. Qs is of opposing sign to the band bending; this means that
    % a negative dphi_b, implying an upward bend and an accumulation of
    % negative charge on the surface, corresponds to a POSITIVE Qs.
    % Therefore reverse the sign.
    deltaQ = -findQsc(dphi_b, T, Ndop, Nv, Nc, epsilon, Eg, ...
    cc_Qsheet_minimal, degeneracyLimit_kT, ...
    ignoreNumericFail, outputFile);     % e/cm^2
    if abs(deltaQ) < cc_Qsheet_minimal     % e/cm^2
        deltaQ = 0;
        dphi_b = 0;
    end
    
    % Output Ef value
    Ef_out = Ef0;
    
    % Dummy assignment
    E_s = [];
end

end