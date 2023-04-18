%% Write CREST output to external file. Arguments:
% outfilePath - path to output file
% crest_out - crest output in matlab struct form
function write_crest_out(outfilePath, crest_out)

% Save structure in external matlab file
save([outfilePath '.mat'], 'crest_out');

% Save structure to xml format
c = struct;
c.crest_out = crest_out;
struct2xml( c, [outfilePath '.xml'])

end
