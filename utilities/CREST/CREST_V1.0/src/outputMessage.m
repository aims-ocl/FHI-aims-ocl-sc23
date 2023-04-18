%% Saves message to file and displays it in command window
% Input:
% messageToAppend - text to append
% filePath - file to append text to
% Option:
% 'silent' - do not display in command window.
% 'new' - do not append, write file from scratch
function outputMessage(messageToAppend, filepath, varargin)

% Display in command window?
optionIndex = find(strcmp('silent', varargin), 1);
if ~isempty(optionIndex)
    silent = 1;
else
    silent = 0;
end

% New file?
optionIndex = find(strcmp('new', varargin), 1);
if ~isempty(optionIndex)
    makeNew = 1;
else
    makeNew = 0;
end

% Very simple write to file
if ~exist(filepath, 'file') || makeNew
    formatSpec = '%s';
    fid = fopen(filepath, 'w');
else
    formatSpec = '\n%s';
    fid = fopen(filepath, 'a');
end
fprintf(fid, formatSpec, messageToAppend);
% Close file
fclose(fid);

% Also display
if ~silent; disp(messageToAppend); end

end