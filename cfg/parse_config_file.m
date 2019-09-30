function cfg = parse_config_file(cfg_type)

    [folder, ~, ~] = fileparts(mfilename('fullpath'));    
    fprintf('Parsing CFG file: %s.json\n', cfg_type);
    
    assert(ismember(cfg_type, {'Sanger'; 'BulkDNA'; 'BulkRNA'; 'BulkRNA_UMI8'; ...
                               'sc10xV2'; 'sc10xV3'; 'sc10xV3_UMI10'; 'scInDropsV2'; 'scInDropsV3'}), ...
           'Unrecognized CFG Type');
    
    cfg_file = sprintf('%s/%s.json', folder, cfg_type);
    assert(exist(cfg_file, 'file') == 2, sprintf('Missing config file: %s', cfg_file));
    cfg = jsondecode(fileread(cfg_file));
    
    if (startsWith(cfg_type, 'Bulk'))
       if (cfg.read_perspective.ShouldComplement == 'N')
           assert(any(strcmp(cfg.trim.Primer5, {'exact'; 'misplaced'})), ...
               'Invalid Primer5 trim setting for idenitfying UMIs in bulk run. Valid settings are [exact, misplaced].');
       else
           assert(any(strcmp(cfg.trim.Primer3, {'exact'; 'misplaced'})), ...
               'Invalid Primer3 trim setting for idenitfying UMIs in bulk run. Valid settings are [exact, misplaced].');
       end
    end
   
end
