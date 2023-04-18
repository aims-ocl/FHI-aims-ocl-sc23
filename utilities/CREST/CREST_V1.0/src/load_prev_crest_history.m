%% Load Qsheet input and output from a previous CREST run via the xml 
% produced in that run.
function [prev_Qsheet_history, prev_output_history, prev_next_Qsheet] = load_prev_crest_history(prevOutputFileXML)

% Check file exists
if ~exist(prevOutputFileXML, 'file')
    % Error
    % Build absolute file path
    absPOFXML = strrep([pwd filesep prevOutputFileXML], '\', '\\');
    err = MException('MATLAB:read_prev_Qsheet_history:FileNotFound', ...
        sprintf('Could not find file %s', absPOFXML));
    throw(err)
end
    
% Load DOM node
crestDOMNode = xmlread(prevOutputFileXML);
% Extract Qsheet input histories
prev_Qsheet_history = cellfun(@str2double, strsplit(char(crestDOMNode.getElementsByTagName('Qsheet_history').item(0).getFirstChild.getData), '\n')');
% Extract output histories
prev_output_history = struct;
prev_output_history.Qsheet = cellfun(@str2double, strsplit(char(crestDOMNode.getElementsByTagName('output_Qsheet_history').item(0).getFirstChild.getData), '\n')');
prev_output_history.Ef = cellfun(@str2double, strsplit(char(crestDOMNode.getElementsByTagName('output_Ef_history').item(0).getFirstChild.getData), '\n')');
prev_output_history.dphi_b = cellfun(@str2double, strsplit(char(crestDOMNode.getElementsByTagName('output_dphi_b_history').item(0).getFirstChild.getData), '\n')');
% Extract next Qsheet guess if available
if ~isempty(crestDOMNode.getElementsByTagName('next_Qsheet').item(0))
    prev_next_Qsheet = str2double(char(crestDOMNode.getElementsByTagName('next_Qsheet').item(0).getFirstChild.getData));
else
    prev_next_Qsheet = [];
end
