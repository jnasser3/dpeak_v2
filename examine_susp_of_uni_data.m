
%% Get UNI data
info = inst_info('{"det_mode":"uni"}');
inst_ids = {info.distil_id};

%% Get the roast level data and then just make a custom brew
lm = mortar.common.Spaces.probe('lm');
lm = lm.asCell;
roast_data = parse_gctx('/cmap/data/build/a2y13q1/zspc.gctx',...
    'rid',lm,...
    'cid',inst_ids);

roast_data = annotate_ds(roast_data,info,'dim','column');

gene_info = parse_tbl('/cmap/data/vdb/chip/L1000_EPSILON.R2.chip','outfmt','record');
lm_idx = strcmp({gene_info.pr_type},'LM');
gene_info = gene_info(lm_idx);
roast_data = annotate_ds(roast_data,gene_info,...
    'dim','row',...
    'keyfield','pr_id');

%% qnorm level
qnorm_data = parse_gctx('/cmap/data/build/a2y13q1/q2norm.gctx',...
    'rid',lm,...
    'cid',inst_ids);

qnorm_data = annotate_ds(qnorm_data,info,'dim','column');
qnorm_data = annotate_ds(qnorm_data,gene_info,...
    'dim','row',...
    'keyfield','pr_id');

%% Brew it
brew_data = distil(roast_data,...
    'group_by','cell_id,pert_id,pert_dose,pert_time',...
    'pert_desc_field','pert_desc');

%% Test 4 replicate subset
four_rep_idx = cell2mat(brew_data.cdesc(:,brew_data.cdict('distil_nsample'))) == 4;
four_rep_brew = ds_slice(brew_data,'cidx',find(four_rep_idx));
susp = compute_misscall_suspects_per_brew_gene(four_rep_brew,roast_data);

run_single_susp(four_rep_brew,susp,10)