function [per_sample_susp,per_sample_susp_in_gene_set] = run_single_susp(brew,susp,susp_cutoff,cell_name,varargin)

params = {'mkfig'};
dflts = {true};
arg = parse_args(params,dflts,varargin{:});

%% Check IDs
common_ids = intersect(brew.cid,susp.cid);
brew = ds_slice(brew,'cid',common_ids);
susp = ds_slice(susp,'cid',common_ids);

assert(isequal(brew.rid,susp.rid))
assert(isequal(brew.cid,susp.cid))

%% Overall CDF
if arg.mkfig
    figure
    [~,~,x,f] = kde(susp.mat(:),2^10,0,50);
    plot(x,f);
    hold on
    title([cell_name ' - Suspiciousness CDF'])
    xlabel('Suspiciousness')
    grid on
    namefig(strcat(cell_name,'_susp_cdf'))
    xlim([3 10])
end

%% Susp per sample
per_sample_susp = sum(susp.mat > susp_cutoff);
med = median(per_sample_susp);
if arg.mkfig
    figure
    histogram(per_sample_susp)%,'normalization','cdf')
    hold on
    ax = gca;
    plot([med med],[0 ax.YLim(2)],'r-')
    title_str = sprintf('%s - Suspect genes per signature\n cutoff = %d \n num_sig = %d\n median = %d',cell_name,susp_cutoff,length(brew.cid),med);
    title(title_str,'Interpreter','none')
    xlabel('Number of Suspect Genes')
    ylabel('Count')
    xlim([0 200])
    namefig(strcat(cell_name,'_susp_per_sig'))
end

%% Number of suspects in top/bottom N genes
N = 50;
M = length(brew.rid);
% [~,sort_idx] = sort(brew.mat);
% sig_idx = sort_idx([1:N (M-N):M],:)
% susp.mat(sig_idx,:)
% size(susp.mat(sig_idx,:))
% per_sample_susp_in_gene_set = sum(susp.mat(sig_idx,:) > susp_cutoff)

per_sample_susp_in_gene_set = zeros(1,length(brew.cid));
for ii = 1:length(brew.cid)
    [~,sort_idx] = sort(brew.mat(:,ii));
    temp_idx = sort_idx([1:N (M-N):M]);
    per_sample_susp_in_gene_set(ii) = sum(susp.mat(temp_idx,ii) > susp_cutoff);
end

if arg.mkfig
    figure
    histogram(per_sample_susp_in_gene_set)%,'normalization','cdf')
    med = median(per_sample_susp_in_gene_set);
    hold on
    ax = gca;
    plot([med med],[0 ax.YLim(2)],'r-')
    title_str = sprintf('%s - Suspect entries in gene set per signature\n cutoff = %d \n num_sig = %d \n median = %d',cell_name,susp_cutoff,length(brew.cid),med);
    title(title_str,'Interpreter','none')
    xlabel('Number of Suspect Genes in top/bottom 50')
    ylabel('Count')
    namefig(strcat(cell_name,'_susp_per_sig_gene_set'))

    %Scatter
    figure
    scatter(per_sample_susp,per_sample_susp_in_gene_set)
    xlabel('Number of suspects per sample')
    ylabel('Number of suspects in gene set per sample')
    title(cell_name)
    ax = gca;
    temp_lim = min(horzcat(ax.YLim(2),ax.XLim(2)));
    hold on
    plot([0 temp_lim],[0 temp_lim],'r:')
    namefig(strcat(cell_name,'_susp_scatter'))
end

