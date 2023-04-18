%% Load CREST input from external file. Arguments:
% infilePath - path to input file
function [U, DOS, Ef0] = read_atomistic2crest(useDOS, atom2crestDir)
% Files suffix
suffix = '.dat';

% Add separator to dir if not included
if ~strcmp(atom2crestDir(end), filesep)
    atom2crestDir = [atom2crestDir filesep];
end

fileNotFound = [];
% Load electrostatic potential energy
if exist([atom2crestDir 'U' suffix], 'file')
    fid = fopen([atom2crestDir 'U' suffix], 'r');
    contents = textscan(fid, '%f\t%f', 'CollectOutput', 1, 'commentStyle', '#', 'MultipleDelimsAsOne', 1);
    U = contents{1};
    fclose(fid);
else
    fileNotFound = [atom2crestDir 'U' suffix];
end

% Load DOS
if exist([atom2crestDir 'DOS' suffix], 'file')
    fid = fopen([atom2crestDir 'DOS' suffix], 'r');
    contents = textscan(fid, '%f\t%f', 'CollectOutput', 1, 'commentStyle', '#', 'MultipleDelimsAsOne', 1);
    DOS = contents{1};
    fclose(fid);
else
    if useDOS
        fileNotFound = [atom2crestDir 'DOS' suffix];
    else
        DOS = [];
    end
end

% Load Ef
if exist([atom2crestDir 'Ef0' suffix], 'file')
    fid = fopen([atom2crestDir 'Ef0' suffix], 'r');
    contents = textscan(fid, '%f', 'CollectOutput', 1, 'commentStyle', '#', 'MultipleDelimsAsOne', 1);
    Ef0 = contents{1};
    fclose(fid);
else
    if useDOS
        Ef0 = [];
    else
        fileNotFound = [atom2crestDir 'Ef0' suffix];
    end
end

if ~isempty(fileNotFound)
    % Error
    err = MException('MATLAB:read_atomistic2crest:FileNotFound', ...
        'Could not find file %s', fileNotFound);
    throw(err)
end

end
