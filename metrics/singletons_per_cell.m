function spc = singletons_per_cell(summary)
    assert(isa(summary, 'ExperimentSummary'));
    spc = sum(summary.allele_freqs==1)/sum(summary.allele_freqs);
end