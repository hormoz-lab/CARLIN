function di = diversity_index(summary, normalize_by_edited)
    assert(isa(summary, 'ExperimentSummary'));
    refseq = CARLIN_def.getInstance.seq.CARLIN;
    is_template = ismember(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), refseq);
    ec = effective_alleles(summary);
    if (normalize_by_edited)
        di = max(0,ec/sum(summary.allele_freqs(~is_template)));
    else
        di = ec/sum(summary.allele_freqs);
    end
end