classdef BulkExperimentReport < ExperimentReport

    methods (Access = public)
        
        function obj = BulkExperimentReport(CARLIN_def, alleles, allele_colony, N, reads)
            obj = obj@ExperimentReport(CARLIN_def, alleles, allele_colony, N, reads);
        end
        
    end
    
    methods (Static)
        
        function obj = create(CARLIN_def, UMI_collection_denoised, UMI_denoise_map, UMI_called_allele, FQ, threshold)
            
            fprintf('Generating report for bulk experiment\n');
            
            [alleles, ~, which_allele] = unique_by_freq(cellfun(@(x) degap(x.get_seq), {UMI_called_allele.allele}', 'un', false));            
            allele_colony = arrayfun(@(i) {UMI_called_allele(which_allele==i).UMI}', [1:size(alleles,1)]', 'un', false);
            template_index = strcmp(alleles, CARLIN_def.seq.CARLIN);
            umi_weights = cellfun(@(x) sum(x), {UMI_collection_denoised.UMIs.SEQ_weight}');
            
            N.uncleaned_tags = length(keys(UMI_denoise_map));
            N.cleaned_tags   = length(unique(values(UMI_denoise_map)));
            N.common_tags    = sum(umi_weights >= threshold.chosen);
            N.called_tags    = size(UMI_called_allele,1);
            N.eventful_tags  = sum(cellfun(@length, allele_colony(~template_index)));
                        
            reads.in_fastq            = FQ.Nreads;
            reads.in_fastq            = FQ.Nreads;
            for fn = fieldnames(FQ.masks)'
                reads.(fn{1}) = length(FQ.masks.(fn{1}));
            end            
            reads.common_tags          = sum(umi_weights(umi_weights >= threshold.chosen));
            
            [is, where] = ismember({UMI_called_allele.UMI}', {UMI_collection_denoised.UMIs.UMI}');
            assert(all(is));
            reads.called_tags_total  = umi_weights(where);
            reads.called_tags_allele = arrayfun(@(i) sum(UMI_collection_denoised.UMIs(where(i)).SEQ_weight(UMI_called_allele(i).constituents)), ...
                                                    [1:length(where)]');
            
            if (any(template_index))
                is = ismember({UMI_called_allele.UMI}', allele_colony{template_index});
                reads.eventful_tags_total = sum(reads.called_tags_total(~is));
                reads.eventful_tags_allele = sum(reads.called_tags_allele(~is));
            else
                reads.eventful_tags_total = sum(reads.called_tags_total);
                reads.eventful_tags_allele = sum(reads.called_tags_allele);
            end
            
            reads.called_tags_total = sum(reads.called_tags_total);
            reads.called_tags_allele = sum(reads.called_tags_allele);
            
            obj = BulkExperimentReport(CARLIN_def, alleles, allele_colony, N, reads);
        end    
        
        function obj = loadobj(s)
            
            if (~isempty(s.CARLIN_def))
                obj = s;
                return;
            end
            
            assert(isa(s, 'BulkExperimentReport'), 'loadobj only defined for existing BulkExperimentReport objects');
            fprintf('Reformatting legacy BulkExperimentReport object\n');
            CARLIN_def = ExperimentSummary.get_legacy_amplicon(s);
            obj = BulkExperimentReport(CARLIN_def, s.alleles, s.allele_colony, s.N, s.reads);
            
        end
        
    end
    
end