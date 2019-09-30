function ref_BCs = get_SC_ref_BCs(barcode_file)
    
    fprintf('Parsing CB reference list %s\n', barcode_file);

    assert(exist(barcode_file, 'file') == 2, 'Missing CB reference file: %s', barcode_file);
    fid = fopen(barcode_file);
    assert(fid ~= -1, 'Failed to open CB reference file: %s', barcode_file);
    
    ref_BCs = textscan(fid, '%s');    
    fclose(fid);
    ref_BCs = ref_BCs{1};
    assert(all(cellfun(@(x) all(ismember(x, 'ACGT')), ref_BCs)), 'Barcodes must only contain {ACGT}');
    assert(length(unique(ref_BCs)) == length(ref_BCs), 'Barcodes should be unique');
    fprintf('%d CBs in reference list\n', size(ref_BCs,1));
    
end