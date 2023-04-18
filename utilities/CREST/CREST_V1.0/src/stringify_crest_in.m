%% Write structure fields in pair format. Only scalars / strings.
function crest_in_str = stringify_crest_in(crest_in)

crest_in_str = '';

allnames = fieldnames(crest_in);
for fi = 1 :numel(allnames)
    if strcmp(allnames{fi}, 'run_electr_structure_command') || ...
      strcmp(allnames{fi}, 'xml_from_previous') || ...
      strcmp(allnames{fi}, 'atomistic_output_dir')
        if ~isempty(crest_in.(allnames{fi}))
            restored_command = strrep(crest_in.(allnames{fi}), '\\', '\\\\');
            restored_command = ['''' strrep(restored_command, '''', '''''') ''''];
        else
            restored_command = [];
        end
        crest_in_str = [crest_in_str sprintf('%s = %s\n', allnames{fi}, restored_command)];
    else
        crest_in_str = [crest_in_str sprintf('%s = %0.7g\n', allnames{fi}, crest_in.(allnames{fi}))];
    end
end

end
