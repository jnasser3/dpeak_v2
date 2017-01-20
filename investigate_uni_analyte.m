function investigate_uni_analyte(brew_prefix)
% Takes as input a brew prefix and outputs diagnostics for analytes 11 and
% 499

qnorm = parse_gctx(strtrim(ls(fullfile('/Volumes/cmap_obelix/pod/custom/REP/brew/pc',brew_prefix,'*_QNORM_*.gctx'))));

gene_info = parse_tbl('/cmap/data/vdb/chip/L1000_EPSILON.R2.chip','outfmt','record');
lm_idx = strcmp({gene_info.pr_type},'LM');
gene_info = gene_info(lm_idx);
qnorm = annotate_ds(qnorm,gene_info,...
    'dim','row',...
    'keyfield','pr_id');

qnorm.rhd

uni_qnorm = ds_subset(qnorm,'row','pr_analyte_id',{'Analyte 11','Analyte 499'});

figure
histogram(uni_qnorm.mat(1,:))
figure
histogram(uni_qnorm.mat(2,:))
end

