%% Runs the command given by the string after substituting the given Qsheet
% in. Input: 
% command - the command to run. Will be parsed.
% Qsheet - will substitute 'xxQsheetxx' in the 'command' parameter.
function runExternalCommand(command, Qsheet)

% Substitute Qsheet in
command2run  = regexprep(command, 'xxQsheetxx', sprintf('%0.20E', Qsheet));

% Run command
status = system(command2run);

% Exit with error if status is not 0
% Load data from external file
if status ~= 0
    err = MException('MATLAB:runExternalCommand:CommandError', ...
        sprintf('Command exited with nonzero status %s', status));
    throw(err)
end

% For TB code running on Matlab
%eval(command2run);

end
