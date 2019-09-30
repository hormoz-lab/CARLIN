function apc = alleles_per_cell(summary)
    assert(isa(summary, 'ExperimentSummary'));
    apc = size(summary.alleles,1)/sum(summary.allele_freqs);
end
