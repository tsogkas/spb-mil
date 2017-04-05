fields_wt = fieldnames(structure);
for k_structure=1:length(fields_wt),
    eval(sprintf('%s = structure.%s;',fields_wt{k_structure},fields_wt{k_structure}));
end