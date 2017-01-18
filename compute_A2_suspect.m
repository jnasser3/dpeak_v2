%Compute suspect zscores on A2

cells = {'A375','A549','HT29','HEPG2','HCC515','HA1E','MCF7','PC3','VCAP'};

for ii = 1:length(cells)
    query_str = sprintf('{"cell_id":"%s","pert_type":"trt_cp"}',cells{ii});
    sig_info = mongo_info('sig_info',query_str,...
        'fields',{'sig_id','cell_id','distil_id','pert_type','distil_nsample'});
    
    three_rep_idx = ([sig_info.distil_nsample] == 3);
    sig_info = sig_info(three_rep_idx);
    
    sig_ids = {sig_info.sig_id};
    lm = mortar.common.Spaces.probe('lm');
    lm = lm.asCell;
    
    %Hack to get instids in the right format
    inst_ids = {};
    for jj = 1:length(sig_ids)
        temp = sig_info(jj).distil_id;
        temp = cell(temp.toArray);
        sig_info(jj).distil_id = strjoin(temp,'|');
        inst_ids = vertcat(inst_ids,temp);
    end
    
    brew_data = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx',...
        'rid',lm,...
        'cid',sig_ids);
    
    brew_data = annotate_ds(brew_data,sig_info);
    
    %Get inst ids from sig_info
    %temp = cellfun(@(x) cell(x.toArray), {sig_info.distil_id}, 'Uni', false);
    %inst_ids = vertcat(sig_info(:).distil_id);
     
    roast_data = parse_gctx('/cmap/data/build/a2y13q1/zspc.gctx',...
        'rid',lm,...
        'cid',inst_ids);
    
    instinfo = inst_info(inst_ids);
    roast_data = annotate_ds(roast_data,instinfo);
    
    gene_info = parse_tbl('/cmap/data/vdb/chip/L1000_EPSILON.R2.chip','outfmt','record');
    lm_idx = strcmp({gene_info.pr_type},'LM');
    gene_info = gene_info(lm_idx);
    roast_data = annotate_ds(roast_data,gene_info,...
        'dim','row',...
        'keyfield','pr_id');

    susp = compute_misscall_suspects_per_brew_gene(brew_data,roast_data);
    
    %%Save files
    mkgctx(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_SUSP.gctx',cells{ii})),susp); 
    mkgctx(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_ROAST.gctx',cells{ii})),roast_data);   
    mkgctx(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_BREW.gctx',cells{ii})),brew_data);
    
    
    mkgct(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_SUSP.gct',cells{ii})),susp); 
    mkgct(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_ROAST.gct',cells{ii})),roast_data);   
    mkgct(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices',sprintf('%s_BREW.gct',cells{ii})),brew_data);
   
end