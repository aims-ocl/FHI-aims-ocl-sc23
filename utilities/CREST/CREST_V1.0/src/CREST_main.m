%% Works out the correct electronic structure, band bending and Fermi level
% of a system with any reasonable degree of charge transfer to the surface.
% Input:
% inputFile - path to file containing properly formatted input flags and data
% Output:
% outputFile - path to file to which to write output and messages

function CREST_main(inputFile, outputFile)

% Begin output messages
outputMessage(sprintf('CREST ver. 1.0; program begins: %s',...
    datestr(clock, 'dd/mm/yyyy HH:MM:SS')), outputFile, 'new');

% REPORT PHYSICAL CONSTANTS USED
physConsts = initPhysicalConsts;
outputMessage(sprintf(['\nUsing the following physical constants:\n' ...
    'Electron charge     = %0.8e Coulomb/electron\n' ...
    'Boltzmann constant  = %0.8e Joule/Kelvin\n' ...
    'Vacuum permittivity = %0.8e e/V-cm' ...
    ],...
    physConsts.q, physConsts.k, physConsts.epsilon_0), outputFile, 'silent');

% PARSE AND UNPACK DATA from CREST user input
crest_in = read_crest_in(inputFile);
% Report input flags
crest_in_str = stringify_crest_in(crest_in);
outputMessage(sprintf('\nRead input flags:\n%s', crest_in_str), outputFile, 'silent');

% Control
atomisticCommand = crest_in.run_electr_structure_command;
atom2crestDir = crest_in.atomistic_output_dir;
calcETotCorrection = crest_in.calculate_energy_correction;
useDOS = crest_in.use_DOS;
extrapolateBulklikeU = crest_in.extrapolate_bulklike_U;
UAveragingInterval = crest_in.U_averaging_interval;     % A
UblAveragingInterval = crest_in.Ubl_averaging_interval;     % A
degeneracyLimit_kT = crest_in.degeneracy_limit;     % kT
ignoreNumericFail = crest_in.fallback_analytical_PBsol;
% SC
cc_Qsheet = crest_in.convergence_criterion_Qsheet;     % Electronic charge
cc_Ef = crest_in.convergence_criterion_Ef;     % eV
cc_dphi_b = crest_in.convergence_criterion_delta_phi_b;     % V
max_nIterations = crest_in.max_n_iterations;
mixingParameter = crest_in.mixing_parameter;
initial_Qsheet = crest_in.initial_Qsheet;
initial_dphi_b = crest_in.initial_delta_phi_b;
prevOutputFileXML = crest_in.xml_from_previous;
% Bulk
T = crest_in.T;     % K
Nv = crest_in.bulk_Nv;     % 1/cm^3
Nc = crest_in.bulk_Nc;     % 1/cm^3
epsilon = crest_in.bulk_epsilon;
Eg = crest_in.bulk_Eg;     % eV
Ndop = crest_in.doping_density;     % e/cm^3
% Non-charged (flatband) reference data
Wv_FP = crest_in.flatband_WF;     % eV
eDelta_phi_v_FP = crest_in.flatband_eDelta_phi_v;     % eV
% Slab geometry & charge
vacZ = crest_in.z_vacuum;     % A
sheetZ = crest_in.z_sheet;     % A
bulklikeZ1 = crest_in.z_bulklike_1;     % A
bulklikeZ2 = crest_in.z_bulklike_2;     % A
bulklikeZ3 = crest_in.z_bulklike_3;     % A
z_d_est = crest_in.z_bisection_estimated;     % A
z_d_uncertainty = crest_in.z_bisection_uncertainty;     % A
surfaceCellArea = crest_in.surface_cell_area * 1e-16;     % A^2 * cm^2/A^2 = cm^2
Qintrinsic = -crest_in.intrinsic_num_electrons;     % e


% PROCESS INPUT DATA
% Ensure bulk-like z's given are sufficiently separated
if extrapolateBulklikeU
    bulklikeZs = sort([bulklikeZ1; bulklikeZ2; bulklikeZ3]);
    if min(diff(bulklikeZs)) < UAveragingInterval
        % Error
        err = MException('MATLAB:CREST_main:BLPositionsTooTight', ...
            'Positions for bulk potential extrapolation closer together than general averaging interval');
        throw(err)
    end
else
    bulklikeZs = [];
end
% Derive convergence criteria
cc_energy = min(cc_dphi_b, cc_Ef);     % eV (V)
dphi_b_widest_bounds = [floor(-1.2*Eg / cc_energy), ceil(1.2*Eg / cc_energy)] * cc_energy;
cc_Qsheet_minimal = estimateQsheetCriterion(cc_dphi_b, cc_Ef, cc_Qsheet, ...
    dphi_b_widest_bounds, [], T, Ndop, Nv, Nc, epsilon, Eg, ...
    false, [], [], surfaceCellArea);
outputMessage(sprintf('Minimal convergence criterion on sheet charge is: %0.3e e/cm^2', ...
    cc_Qsheet_minimal), outputFile, 'silent');
% Load previous history if available
if ~isempty(prevOutputFileXML)
    [prev_Qsheet_history, prev_output_history, ~] = load_prev_crest_history(prevOutputFileXML);
    prev_output_Qsheet_history = prev_output_history.Qsheet;
    prev_output_Ef_history = prev_output_history.Ef;
    prev_output_dphi_b_history = prev_output_history.dphi_b;
else
    prev_Qsheet_history = [];
    prev_output_Qsheet_history = [];
    prev_output_Ef_history = [];
    prev_output_dphi_b_history = [];
end


% Prepare to save history
Qsheet_history = zeros(max_nIterations, 1);
output_Qsheet_history = zeros(max_nIterations, 1);
output_Ef_history = zeros(max_nIterations, 1);
output_dphi_b_history = zeros(max_nIterations, 1);
atomistic_data_history = cell(max_nIterations, 1);

% Determine initial guess for sheet charge
if ~isempty(initial_Qsheet)
    next_Qsheet = initial_Qsheet;     % e/cm^2
else
    next_Qsheet = findQsc(initial_dphi_b, T, Ndop, Nv, Nc, epsilon, Eg, ...
        cc_Qsheet_minimal, degeneracyLimit_kT, ...
        ignoreNumericFail, outputFile);     % e/cm^2
end

converged = false;
for iter_num = 1 : max_nIterations
    % Use the Qsheet mixed in the previous step as input
    Qsheet = next_Qsheet;
    outputMessage(sprintf(['\nBeginning iteration %d: %s\n' ...
        '      Input Qsheet = %0.10e e/cm^2'], ...
        iter_num, datestr(clock, 'dd/mm HH:MM:SS'), Qsheet), outputFile, 'silent');
    % Save Qsheet history
    Qsheet_history(iter_num) = Qsheet;
    
    
    % SETUP, SOLVE AND SAVE INFO FROM system with given charged sheet
    outputMessage('Running atomistic code', outputFile, 'silent');
    runExternalCommand(atomisticCommand, Qsheet);
    outputMessage('Loading data from atomistic code', outputFile, 'silent');
    [U, enLev_DOS, Ef0] = read_atomistic2crest(useDOS, atom2crestDir);
    
    
    % FIND DELTA Q (EQUILIBRIUM CHARGE TRANSFER) for this atomistic data
    [deltaQ, dphi_b, Ef_out, z_d, U_d, E_s, E_b] = calcChargeTransfer(U, enLev_DOS, Ef0, Qsheet, ...
        useDOS, extrapolateBulklikeU, UAveragingInterval, UblAveragingInterval, ...
        degeneracyLimit_kT, ignoreNumericFail, cc_Qsheet_minimal, cc_dphi_b, cc_Ef, cc_Qsheet, ...
        T, Nv, Nc, epsilon, Eg, Ndop, Wv_FP, eDelta_phi_v_FP, ...
        vacZ, bulklikeZs, z_d_est, surfaceCellArea, Qintrinsic, outputFile);
    
    % Convert slab net charge to sheet charge
    output_Qsheet = -deltaQ;     % e/cm^2
    
    % Determine Qsheet convergence criterion for this step
    dphi_b_bounds = [max(dphi_b_widest_bounds(1), dphi_b - 5*cc_energy), ...
        min(dphi_b_widest_bounds(2), dphi_b + 5*cc_energy)];
    Ef_bounds = [];
    if useDOS
        Ef_bounds = [max(min(E_s, E_b)-0.2*Eg, Ef_out - 5*cc_Ef), ...
            min(max(E_s, E_b)+0.2*Eg, Ef_out + 5*cc_Ef)];
    end
    cc_Qsheet_thisstep = estimateQsheetCriterion(cc_dphi_b, cc_Ef, cc_Qsheet, ...
        dphi_b_bounds, Ef_bounds, T, Ndop, Nv, Nc, epsilon, Eg, ...
        useDOS, enLev_DOS, E_s, surfaceCellArea);
    
    % Report iteration results
    outputMessage(sprintf(['Iteration %d results:\n' ...
        '      Output Qsheet = %0.10e e/cm^2\n' ...
        '      Ef = %0.6f eV\n' ...
        '      delta_phi_b = %0.6f V\n' ...
        '      Conv. criterion on Qsheet in this step = %0.3e e/cm^2'], ...
        iter_num, output_Qsheet, Ef_out, dphi_b, cc_Qsheet_thisstep), outputFile, 'silent');
    
    % Update output history
    output_Qsheet_history(iter_num) = output_Qsheet;
    output_Ef_history(iter_num) = Ef_out;
    output_dphi_b_history(iter_num) = dphi_b;
    atomistic_data_history{iter_num} = struct;
    atomistic_data_history{iter_num}.U = U;
    if ~isempty(enLev_DOS); atomistic_data_history{iter_num}.DOS = enLev_DOS; end
    if ~isempty(Ef0); atomistic_data_history{iter_num}.Ef = Ef0; end
    
    
    % CHECK CONVERGENCE
    % Combine current and previous run histories
    combined_Qsheet_history = [prev_Qsheet_history; Qsheet_history(1 : iter_num)];
    combined_output_Qsheet_history = [prev_output_Qsheet_history; output_Qsheet_history(1 : iter_num)];
    combined_output_Ef_history = [prev_output_Ef_history; output_Ef_history(1 : iter_num)];
    combined_output_dphi_b_history = [prev_output_dphi_b_history; output_dphi_b_history(1 : iter_num)];
    
    % Two conditions for convergence:
    % 1. That the input and output Qsheet values in this step
    % match within the convergence criterion for Qsheet in this step
    QsheetConsistent = abs(combined_Qsheet_history(end) - combined_output_Qsheet_history(end)) < cc_Qsheet_thisstep;
    
    % 2. That the difference between the outputs of the current and
    % previous steps is less than the convergence criteria specified
    % (ignored if this is the first step)
    if numel(combined_output_Qsheet_history) < 2
        outputConverged = true;
    else
        outputConverged = abs(diff(combined_output_Qsheet_history(end-1 : end))) < cc_Qsheet && ...
            abs(diff(combined_output_Ef_history(end-1 : end))) < cc_Ef && ...
            abs(diff(combined_output_dphi_b_history(end-1 : end))) < cc_dphi_b;
    end
    
    if QsheetConsistent && outputConverged
        % Perform series of tests to ensure convergence
        converged = true;
        
        if mixingParameter < 0 && mixingParameter > -1
            % Must first reach a mixing parameter of -1 or >0
            mixingParameter = -1;
            converged = false;
        elseif mixingParameter == -1 || mixingParameter < 1
            % Then must reach a mixing parameter of 1
            mixingParameter = 1;
            converged = false;
        end
        
        if converged
            % Finally converged
            outputMessage(sprintf('\nCONVERGENCE REACHED'), outputFile);
            finalIteration = iter_num;
            break;
        end
    end
    
    % Didn't break from the loop, so did not converge. Mix next Qsheet.
    next_Qsheet = mixQsheet(combined_Qsheet_history, combined_output_Qsheet_history, mixingParameter);
    
    % If we are at the maximum no. of iterations, stop loop
    if iter_num == max_nIterations
        outputMessage(sprintf('\nWARNING: CREST SELF-CONSISTENT CYCLE DID NOT CONVERGE'), outputFile);
        finalIteration = iter_num;
        break;
    end
end


% CALCULATE TOTAL ENERGY CORRECTION
if converged && calcETotCorrection
    outputMessage('Calculating total energy correction', outputFile);
    if ~isempty(z_d)
        [Qsc, z_SCR, V_SCR, E_SCR, ~, trunc_est_from, errmsg] = solvePoisson_numeric(dphi_b, ...
            T, Ndop, Nv, Nc, epsilon, Eg, degeneracyLimit_kT, output_Qsheet, ...
            'outgridnpts', 10001);
        if ~strcmp(errmsg, '')
            outputMessage(sprintf('Numeric Poisson solver may have failed. Error message:\n%s', ...
                errmsg), outputFile);
        elseif abs(Qsc - output_Qsheet) > cc_Qsheet_thisstep
            % Warn. If numerical solver failed, this need not be checked
            % since in the estimation, Qsc == output_Qsheet by definition.
            outputMessage(sprintf(['Warning: Numerical output for SCR charge does not match expected value.\n' ...
                '  delta-phi_b obtained in last iteration (input to solver): %0.6f V\n' ...
                '  Q_SCR result from solver:                                 %0.10e e/cm^2\n' ...
                '  Q_sheet obtained in last iteration:                       %0.10e e/cm^2' ...
                ], dphi_b, Qsc, output_Qsheet), outputFile);
        end
        if ~isempty(trunc_est_from)
            if trunc_est_from == 0
                outputMessage(sprintf(['Estimated all potential development and '...
                    'total energy from Qsheet and band-bending']), outputFile);
            else
                outputMessage(sprintf(['Very wide SCR! potential development more than ' ...
                    '%0.2e %s from z_d approximated.'], ...
                    abs(trunc_est_from / 1e-8), char(197)), outputFile);
            end
        end
        Utot_correction = calcTotalEnergyCorrection(epsilon, z_SCR, E_SCR, Qsheet, z_d, sheetZ);
        U_SCR = [z_SCR*1e8 -V_SCR];
    else
        Utot_correction = 0;
        U_SCR = [0 0];
    end
end


% ESTIMATE ERROR CAUSED BY UNCERTAINTY IN Z_D
% Relevant only if:
% 1. The calculation converged
% 2. The DOS was used to facilitate convergence; otherwise it is impossible
% to estimate the error in values that depend on both the slab and the bulk
% without performing 2 additional atomistic calculations. In this case the
% user must manage manually
if converged && useDOS
    outputMessage('Estimating uncertainty in results', outputFile);
    Qsheet_uncertainty = 0;
    dphi_b_uncertainty = 0;
    Ef_out_uncertainty = 0;
    if isempty(z_d)
        z_d_uncertainty = [];
    elseif z_d_uncertainty > 0
        Qsheet_uncertainty_limits = ones(3, 1) * Qsheet;
        dphi_b_uncertainty_limits = ones(3, 1) * dphi_b;
        Ef_out_uncertainty_limits = ones(3, 1) * Ef_out;
        for li = [1 3]
            z_d_try = z_d + (li-2) * z_d_uncertainty;
            [deltaQ_bound, dphi_b_bound, Ef_out_bound, ~, ~, ~, ~] = calcChargeTransfer(U, enLev_DOS, Ef0, Qsheet, ...
                useDOS, false, UAveragingInterval, UblAveragingInterval, ...
                degeneracyLimit_kT, ignoreNumericFail, cc_Qsheet_minimal, cc_dphi_b, cc_Ef, cc_Qsheet, ...
                T, Nv, Nc, epsilon, Eg, Ndop, Wv_FP, eDelta_phi_v_FP, ...
                vacZ, bulklikeZs, z_d_try, surfaceCellArea, Qintrinsic, outputFile);
            Qsheet_uncertainty_limits(li) = -deltaQ_bound;
            dphi_b_uncertainty_limits(li) = dphi_b_bound;
            Ef_out_uncertainty_limits(li) = Ef_out_bound;
        end
        Qsheet_uncertainty = max(abs(diff(Qsheet_uncertainty_limits)));
        dphi_b_uncertainty = max(abs(diff(dphi_b_uncertainty_limits)));
        Ef_out_uncertainty = max(abs(diff(Ef_out_uncertainty_limits)));
    end
    outputMessage(sprintf(' '), outputFile);
end


% Pack up output data in struct
crest_out = struct;
crest_out.converged = converged;
crest_out.finalIteration = finalIteration;
if converged
    crest_out.Qsheet = Qsheet;
    crest_out.Ef = Ef_out;
    crest_out.dphi_b = dphi_b;
    crest_out.z_d = z_d;
    crest_out.U_d = U_d;
    if useDOS
        crest_out.z_d_uncertainty = z_d_uncertainty;
        crest_out.Qsheet_uncertainty = Qsheet_uncertainty;
        crest_out.dphi_b_uncertainty = dphi_b_uncertainty;
        crest_out.Ef_uncertainty = Ef_out_uncertainty;
    end
    if calcETotCorrection
        crest_out.Utot_correction = Utot_correction;
        crest_out.U_SCR = U_SCR;
    end
else
    crest_out.next_Qsheet = next_Qsheet;
end
crest_out.Qsheet_history = Qsheet_history(1 : finalIteration);
crest_out.output_Qsheet_history = output_Qsheet_history(1 : finalIteration);
crest_out.output_Ef_history = output_Ef_history(1 : finalIteration);
crest_out.output_dphi_b_history = output_dphi_b_history(1 : finalIteration);
crest_out.atomistic_data_history = atomistic_data_history(1 : finalIteration);
% Write to file(s)
write_crest_out(outputFile, crest_out);

% Report final results
if converged
    if ~isempty(surfaceCellArea)
        finalPerCellQsheetText = sprintf(' = %0.4e e/cell', ...
            crest_out.Qsheet * surfaceCellArea);     % e/cm^2 * cm^2 = e
    else
        finalPerCellQsheetText = '';
    end
    outputMessage(sprintf(['Final results:\n' ...
        'Qsheet      = %0.4e e/cm^2%s\n' ...
        'E_F         = %0.4f eV\n' ...
        'Delta_phi_b = %0.4f V' ...
        ], crest_out.Qsheet, finalPerCellQsheetText, ...
        crest_out.Ef, crest_out.dphi_b), outputFile);
    if calcETotCorrection
        if ~isempty(surfaceCellArea)
            finalPerCellTotEnText = sprintf(' = %0.4e eV/cell', ...
                crest_out.Utot_correction * surfaceCellArea);     % eV/cm^2 * cm^2 = eV
        else
            finalPerCellTotEnText = '';
        end
        outputMessage(sprintf( ...
            'Total energy correction = %0.4e eV/cm^2%s', ...
            crest_out.Utot_correction, ...
            finalPerCellTotEnText), outputFile);
    end
    if useDOS && ~isempty(z_d)
        outputMessage(sprintf(['Estimated uncertainty due to z_d uncertainty of %0.3f %s:\n' ...
            '  Qsheet -------> +/- %0.1e e/cm^2\n' ...
            '  E_F ----------> +/- %0.1e eV\n' ...
            '  Delta_phi_b --> +/- %0.1e V' ...
            ], z_d_uncertainty, char(197), crest_out.Qsheet_uncertainty, ...
            crest_out.Ef_uncertainty, crest_out.dphi_b_uncertainty), ...
            outputFile);
    end
else
    outputMessage(sprintf('Next Qsheet guess = %0.12e e/cm^2\n', ...
        crest_out.next_Qsheet), outputFile, 'silent');
end

% Report completion time
outputMessage(sprintf('\nCREST program completed: %s', datestr(clock, 'dd/mm/yyyy HH:MM:SS')), outputFile);

end
