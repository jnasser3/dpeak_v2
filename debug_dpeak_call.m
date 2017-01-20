function debug_dpeak_call(basedir,brew_prefix,analyte,well)
%Make diagnostic plots to debug a dpeak call

%Wants
%lxb histogram
%GEX scatter plots
%QNORM platewide plots


mad_factor = 1.4826;

%Find number of replicates we have
brewdata = parse_gctx(strtrim(ls(fullfile(basedir,brew_prefix,'by_rna_well','*MODZ*SCORE*LM*.gctx'))));

this_data = ds_subset(brewdata,'column','rna_well',well);
num_reps = cell2mat(this_data.cdesc(:,this_data.cdict('distil_nsample')));
temp = this_data.cdesc(:,this_data.cdict('det_plate'));
roast_plates = strsplit(char(temp),'|');

%Get and plot LXB files
dpeak_title_str = cell(1,num_reps);
for ii = 1:num_reps
    %subplot(3,num_reps,ii)
    lxb = parse_lxb(fullfile('/Volumes/cmap_obelix/pod/custom/REP/lxb/',roast_plates{ii},sprintf('%s_%s.lxb',roast_plates{ii},well)));
    x = lxb.RP1(lxb.RID == analyte);
    [~,~,dpeak_title_str{ii}] = detect_lxb_peaks_single(x,'showfig','Off','pkmethod','kmeans_jn');
    ax(ii) = gca;
    ax(ii).XLim = [0 15];
    %title(roast_plates{ii})
end

figure;

for ii = 1:num_reps
    if ii == 1
        ha = subplot(4,3,1);
    elseif ii == 2
        ha = subplot(4,3,2);
    elseif ii == 3 
        ha = subplot(4,3,3);
    end
    title(sprintf('%s \n %s',roast_plates{ii},dpeak_title_str{ii}),'Interpreter','none')
    copyobj(allchild(ax(ii)),ha)
end

%Make GEX scatter plots
for ii = 1:num_reps
    gex = parse_gctx(strtrim(ls(fullfile(basedir,roast_plates{ii},'*GEX*.gct*'))));
    this_gex = ds_subset(gex,'row','pr_analyte_id',sprintf('Analyte %d',analyte));
    well_idx = find(strcmp(well,this_gex.cdesc(:,this_gex.cdict('rna_well'))));
    %figure('Visible','Off')

    %crazy stupid hack.
    %subplot(3,3,3+ii) doesn't work
    if ii == 1
        subplot(4,3,4)
    elseif ii == 2
        subplot(4,3,5)
    elseif ii == 3 
        subplot(4,3,6)
    end
    scatter(safelog2(this_gex.mat(1,:)),safelog2(this_gex.mat(2,:)));
    hold on
    scatter(safelog2(this_gex.mat(1,well_idx)),safelog2(this_gex.mat(2,well_idx)),...
        [],'r*');
    hold on
    plot([0 15],[0 15],'r:')
    gene1 = this_gex.rdesc(1,this_gex.rdict('pr_gene_symbol'));
    gene2 = this_gex.rdesc(2,this_gex.rdict('pr_gene_symbol'));
    bset1 = this_gex.rdesc(1,this_gex.rdict('pr_bset_id'));
    bset2 = this_gex.rdesc(2,this_gex.rdict('pr_bset_id'));
    
    xlabel([gene1{1} ' - ' bset1{1}])
    ylabel([gene2{1} ' - ' bset2{1}])
    %ylabel(this_gex.rdesc(2,this_gex.rdict('pr_gene_symbol')))
    grid on
    axis([0 15 0 15])
    title('GEX')
    
end

%Make QNORM plots
for ii = 1:num_reps
    qnorm = parse_gctx(strtrim(ls(fullfile(basedir,roast_plates{ii},'*QNORM*.gct*'))));
    this_qnorm = ds_subset(qnorm,'row','pr_analyte_id',sprintf('Analyte %d',analyte));
    expr = ds_subset(this_qnorm,'column','rna_well',well);
    if ii == 1
        subplot(4,3,7)
    elseif ii == 2
        subplot(4,3,8)
    elseif ii == 3 
        subplot(4,3,9)
    end
    data = this_qnorm.mat(1,:);
    histogram(data)
    ax = gca;
    hold on
    plot([expr.mat(1,1) expr.mat(1,1)],[0 ax.YLim(2)],'r-')
    gene_name = this_qnorm.rdesc(1,this_qnorm.rdict('pr_gene_symbol'));
    title_str = sprintf('%s \n qnorm = %.2f median = %.2f MAD = %.2f zscore = %.2f',...
        gene_name{1},...
        expr.mat(1,1),median(data),mad(data,1),(expr.mat(1,1) - median(data))/(mad_factor*mad(data,1)));
    title(title_str)
    xlim([0 15])
    
    
    if ii == 1
        subplot(4,3,10)
    elseif ii == 2
        subplot(4,3,11)
    elseif ii == 3 
        subplot(4,3,12)
    end
    data = this_qnorm.mat(2,:);
    histogram(data)
    ax = gca;
    hold on
    plot([expr.mat(2,1) expr.mat(2,1)],[0 ax.YLim(2)],'r-')
    gene_name = this_qnorm.rdesc(2,this_qnorm.rdict('pr_gene_symbol'));
    title_str = sprintf('%s \n qnorm = %.2f median = %.2f MAD = %.2f zscore = %.2f',...
        gene_name{1},...
        expr.mat(2,1),median(data),mad(data,1),(expr.mat(2,1) - median(data))/(mad_factor*mad(data,1)));
    title(title_str)
    xlim([0 15])
end


end
