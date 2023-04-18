%% Load CREST input from external file. Arguments:
% infilePath - path to input file
function crest_in = read_crest_in(infilePath)

% Required flags
requiredFlags = {...
    'run_electr_structure_command'; ...
    'atomistic_output_dir'; ...
    'T'; ...
    'bulk_Nv'; ...
    'bulk_Nc'; ...
    'bulk_epsilon'; ...
    'bulk_Eg'; ...
    'doping_density'; ...
    'flatband_WF'; ...
    'z_vacuum'; ...
    };

% Conditionally required flags
conditionalFlags = {...
    'flatband_eDelta_phi_v', 'extrapolate_bulklike_U', true; ...
    'z_sheet', 'calculate_energy_correction', true; ...
    'z_bulklike_1', 'extrapolate_bulklike_U', true; ...
    'z_bulklike_2', 'extrapolate_bulklike_U', true; ...
    'z_bulklike_3', 'extrapolate_bulklike_U', true; ...
    'z_bisection_estimated', 'extrapolate_bulklike_U', false; ...
    'surface_cell_area', 'use_DOS', true; ...
    'intrinsic_num_electrons', 'use_DOS', true; ...
    };

% Flags with default values
defaultFlags = {...
    'calculate_energy_correction', true; ...
    'use_DOS', true; ...
    'extrapolate_bulklike_U', true; ...
    'U_averaging_interval', 0.5; ...
    'Ubl_averaging_interval', []; ...
    'degeneracy_limit', 3.5; ...
    'fallback_analytical_PBsol', true; ...
    'convergence_criterion_Qsheet', Inf; ...
    'convergence_criterion_Ef', 1e-3; ...
    'convergence_criterion_delta_phi_b', 1e-3; ...
    'max_n_iterations', 20; ...
    'mixing_parameter', 1; ...
    'initial_Qsheet', []; ...
    'initial_delta_phi_b', []; ...
    'xml_from_previous', []; ...
    'z_bisection_uncertainty', 2; ...
    };

% Load data from external file
if ~exist(infilePath, 'file')
    % Error
    err = MException('MATLAB:read_crest_in:FileNotFound', ...
        sprintf('Could not find file %s', infilePath));
    throw(err)
end
infid = fopen(infilePath, 'r');
crest_in = struct;
parseErr = false;
while true
    lineText = fgets(infid);
    % End of lines
    if lineText == -1; break; end
    % Remove remarks
    noRemarksText = regexprep(lineText, '\s*#.*', '');
    % Check for empty line
    if isempty(regexprep(noRemarksText, '\s', '')); continue; end
    % Check for flag format
    flag_parts = regexp(noRemarksText, '\s*=\s*', 'split');
    if numel(flag_parts) ~= 2; parseErr = true; break; end
    % Save flag components without leading/trailing spaces
    flag_name = regexprep(flag_parts{1}, '^\s+|\s+$', '');
    flag_value = regexprep(flag_parts{2}, '^\s+|\s+$', '');
    % Handle specific non-number flags
    if strcmp(flag_name, 'run_electr_structure_command') || ...
            strcmp(flag_name, 'xml_from_previous') || ...
            strcmp(flag_name, 'atomistic_output_dir')
        crest_in.(flag_name) = eval(flag_value);
    elseif strcmp(flag_name, 'calculate_energy_correction') || ...
            strcmp(flag_name, 'use_DOS') || ...
            strcmp(flag_name, 'extrapolate_bulklike_U')
        if strcmp(flag_value, 'true') || strcmp(flag_value, '1');
            crest_in.(flag_name) = true;
        elseif strcmp(flag_value, 'false') || strcmp(flag_value, '0');
            crest_in.(flag_name) = false;
        else
            parseErr = true; break
        end
    else
        try
            crest_in.(flag_name) = str2double(flag_value);
        catch s2dErr
            parseErr = true; break
        end
        if isnan(crest_in.(flag_name))
            crest_in.(flag_name) = [];
        end
    end
end
if parseErr
    fclose(infid);
    err = MException('MATLAB:read_crest_in:ParseError', ...
        ['Could not parse string "' noRemarksText '" in CREST input file']);
    throw(err)
end
fclose(infid);

% Set defaults
for di = 1 : size(defaultFlags, 1)
    if ~isfield(crest_in, defaultFlags{di, 1}) || isempty(crest_in.(defaultFlags{di, 1}))
        crest_in.(defaultFlags{di, 1}) = defaultFlags{di, 2};
    end
end
% Set Ubl averaging interval default
if isempty(crest_in.Ubl_averaging_interval);
    crest_in.Ubl_averaging_interval = crest_in.U_averaging_interval;
end
% Initial guesses
if isempty(crest_in.initial_Qsheet) && isempty(crest_in.initial_delta_phi_b)
    % No initial guess specified explicitly.
    default_initial_Qsheet = 0;
    % Update default initial Qsheet by with guess from file, if available
    if ~isempty(crest_in.xml_from_previous)
        [~, ~, prev_next_Qsheet] = load_prev_crest_history(crest_in.xml_from_previous);
        if ~isempty(prev_next_Qsheet)
            default_initial_Qsheet = prev_next_Qsheet;
        end
    end
    % Set default value
    crest_in.initial_Qsheet = default_initial_Qsheet;
elseif ~isempty(crest_in.initial_Qsheet) && ~isempty(crest_in.initial_delta_phi_b)
    % Initial guesses are over-determined. Error
    err = MException('MATLAB:read_crest_in:FlagForbidden', ...
        ['Problem over-determined: both flag "initial_Qsheet"' ...
        '" and flag "initial_delta_phi_b" specified in CREST input file']);
    throw(err)
end

% Add conditionally required flags to required list
for ci = 1 : size(conditionalFlags, 1)
    if crest_in.(conditionalFlags{ci, 2}) == conditionalFlags{ci, 3}
        % Flag required. Add to list
        requiredFlags{end+1, 1} = conditionalFlags{ci, 1};
    elseif ~isfield(crest_in, conditionalFlags{ci, 1})
        % Ensure flag existence
        crest_in.(conditionalFlags{ci, 1}) = [];
    end
end

% Ensure all required flags exist and have non-empty values
for ri = 1 : size(requiredFlags, 1)
    if ~isfield(crest_in, requiredFlags{ri}) || isempty(crest_in.(requiredFlags{ri}))
        err = MException('MATLAB:read_crest_in:FlagRequired', ...
            ['Could not find/resolve value of flag "' requiredFlags{ri} '" in CREST input file']);
        throw(err)
    end
end

end
