function debug_dpeak_call(varargin)
%Make diagnostic plots to debug a dpeak call

pnames = {'roastdir','brewdir','brew_prefix','analyte','well','dpeak_version'};
dflts =  {'','', '', 0, '', 'v2' };
arg = parse_args(pnames, dflts,varargin{:});      

%Wants
%lxb histogram
%GEX scatter plots
%QNORM platewide plots

mad_factor = 1.4826;

%Find number of replicates we have
brewdata = parse_gctx(strtrim(ls(fullfile(arg.brewdir,arg.brew_prefix,'by_rna_well','*MODZ*SCORE*LM*.gctx'))));

this_brew_data = ds_subset(brewdata,'column','rna_well',arg.well);
num_reps = cell2mat(this_brew_data.cdesc(:,this_brew_data.cdict('distil_nsample')));
temp = this_brew_data.cdesc(:,this_brew_data.cdict('det_plate'));
roast_plates = strsplit(char(temp),'|');

%Get and plot LXB files
dpeak_title_str = cell(1,num_reps);
for ii = 1:num_reps
    %subplot(3,num_reps,ii)
    lxb = parse_lxb(fullfile('/Volumes/cmap_obelix/pod/custom/REP/lxb/',roast_plates{ii},sprintf('%s_%s.lxb',roast_plates{ii},arg.well)));
    x = lxb.RP1(lxb.RID == arg.analyte);
    
    %load random seed to replicate production
    load(fullfile(arg.roastdir,roast_plates{ii},'dpeak','rndseed.mat'))
    switch arg.dpeak_version
        case 'v1'
            [~,~,dpeak_title_str{ii}] = detect_lxb_peaks_single(x,...
                'showfig','Off',...
                'mkfig',true,...
                'pkmethod','kmeans_viable',...
                'lowthresh',4,...
                'highthresh',15,...
                'use_dbscan',false);
        case 'v2'
            [~,~,dpeak_title_str{ii}] = detect_lxb_peaks_single(x,...
                'showfig','Off',...
                'mkfig',true,...
                'pkmethod','kmeans_jn',...
                'lowthresh',1,...
                'highthresh',15,...
                'use_dbscan',true);
    end
    ax(ii) = gca;
    ax(ii).XLim = [0 15];
    %title(roast_plates{ii})
end

figure;

%Copy the lxb figures into the big subplot
for ii = 1:num_reps
    if ii == 1
        ha = subplot(5,3,1);
    elseif ii == 2
        ha = subplot(5,3,2);
    elseif ii == 3 
        ha = subplot(5,3,3);
    end
    title(sprintf('%s \n %s',roast_plates{ii},dpeak_title_str{ii}),'Interpreter','none')
    copyobj(allchild(ax(ii)),ha)
end

%Make GEX scatter plots
for ii = 1:num_reps
    gex = parse_gctx(strtrim(ls(fullfile(arg.roastdir,roast_plates{ii},'*GEX*.gct*'))));
    this_gex = ds_subset(gex,'row','pr_analyte_id',sprintf('Analyte %d',arg.analyte));
    well_idx = find(strcmp(arg.well,this_gex.cdesc(:,this_gex.cdict('rna_well'))));
    
    %crazy stupid hack.
    %subplot(3,3,3+ii) doesn't work
    if ii == 1
        subplot(5,3,4)
    elseif ii == 2
        subplot(5,3,5)
    elseif ii == 3 
        subplot(5,3,6)
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
    ax = gca;
    ax.XTick = [1:15];
    ax.YTick = [1:2:15];
    title('GEX')
    
%end

    %Make QNORM plots
%for ii = 1:num_reps
    qnorm = parse_gctx(strtrim(ls(fullfile(arg.roastdir,roast_plates{ii},'*QNORM*.gct*'))));
    this_qnorm = ds_subset(qnorm,'row','pr_analyte_id',sprintf('Analyte %d',arg.analyte));
    expr = ds_subset(this_qnorm,'column','rna_well',arg.well);
    
    this_well_gex = ds_subset(gex,'column','rna_well',arg.well);
    this_well_qnorm = ds_subset(qnorm,'column','rna_well',arg.well);
    assert(isequal(this_well_gex.rid,this_well_qnorm.rid))
    if ii == 1
        subplot(5,3,7)
    elseif ii == 2
        subplot(5,3,8)
    elseif ii == 3 
        subplot(5,3,9)
    end
    scatter(safelog2(this_well_gex.mat),this_well_qnorm.mat)
    grid on
    title(sprintf('GEX vs QNORM in Well %s',arg.well))
    xlabel('log2(GEX)')
    ylabel('QNORM')
    
    if ii == 1
        subplot(5,3,10)
    elseif ii == 2
        subplot(5,3,11)
    elseif ii == 3 
        subplot(5,3,12)
    end
    data = this_qnorm.mat(1,:);
    histogram(data)
    ax = gca;
    hold on
    plot([expr.mat(1,1) expr.mat(1,1)],[0 ax.YLim(2)],'r-')
    gene_name1 = this_qnorm.rdesc(1,this_qnorm.rdict('pr_gene_symbol'));
    title_str = sprintf('%s - Platewide QNORM \n qnorm = %.2f median = %.2f MAD = %.2f zscore = %.2f',...
        gene_name1{1},...
        expr.mat(1,1),median(data),mad(data,1),(expr.mat(1,1) - median(data))/(mad_factor*mad(data,1)));
    title(title_str)
    xlim([0 15])
    ax = gca;
    ax.XTick = [1:15];
    
    
    if ii == 1
        subplot(5,3,13)
    elseif ii == 2
        subplot(5,3,14)
    elseif ii == 3 
        subplot(5,3,15)
    end
    data = this_qnorm.mat(2,:);
    histogram(data)
    ax = gca;
    hold on
    plot([expr.mat(2,1) expr.mat(2,1)],[0 ax.YLim(2)],'r-')
    gene_name2 = this_qnorm.rdesc(2,this_qnorm.rdict('pr_gene_symbol'));
    title_str = sprintf('%s - Platewide QNORM \n qnorm = %.2f median = %.2f MAD = %.2f zscore = %.2f',...
        gene_name2{1},...
        expr.mat(2,1),median(data),mad(data,1),(expr.mat(2,1) - median(data))/(mad_factor*mad(data,1)));
    title(title_str)
    xlim([0 15])
    ax = gca;
    ax.XTick = [1:15];
%end

end

% Get brew level stats
gene_brew_stats = ds_subset(this_brew_data,'row','pr_analyte_id',sprintf('Analyte %d',arg.analyte));
gene_brew_zscore = gene_brew_stats.mat;
gene1_zscore = gene_brew_zscore(1);
gene2_zscore = gene_brew_zscore(2);
this_brew_wts = gene_brew_stats.cdesc(1,gene_brew_stats.cdict('distil_wt'));


[~,l] = suplabel(sprintf('%s - Well: %s - Analyte: %d \n Brew z-scores: %s = %.2f, %s = %.2f \n modz wts = %s, dpeak_version: %s',...
    arg.brew_prefix,arg.well,arg.analyte,gene_name1{1},gene1_zscore,gene_name2{1},gene2_zscore,this_brew_wts{1},arg.dpeak_version),...
    't');
set(l,'FontSize',20)
set(l,'Interpreter','none')

end
