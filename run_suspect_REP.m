%Compute suspect zscores on A2

cells = {'A375','YAPC','MCF7','HA1E','PC3','HELA','HT29'};

susp_cutoff = 5;
num_reps = 3;

%close all

figure(1); ax1 = gca;
figure(2); ax2 = gca;
figure(3); ax3 = gca;

basedir = '/cmap/projects/repurposing/jnwork/batch_effect/REP_brew_level_suspects';

for ii = 1:length(cells)
    
    this_susp = parse_gctx(strtrim(ls(fullfile(basedir,sprintf('*suspects_%s*.gctx',cells{ii})))));
    this_brew = parse_gctx(strtrim(ls(fullfile(basedir,sprintf('*BREW_%s*.gctx',cells{ii})))));
    
    has_num_rep = cell2mat(this_brew.cdesc(:,this_brew.cdict('distil_nsample'))) == num_reps;
    this_brew = ds_slice(this_brew,'cidx',find(has_num_rep));
    
    has_num_rep2 = cell2mat(this_susp.cdesc(:,this_susp.cdict('distil_nsample'))) == num_reps;
    this_susp = ds_slice(this_susp,'cidx',find(has_num_rep2));
    
    %[f,x] = ecdf(this_susp.mat(:));
    figure(1)
    [~,~,x,f] = kde(this_susp.mat(:),2^10,0,50);
    plot(ax1,x,f,'DisplayName',cells{ii});
    hold on
   
    %cells(ii)
    %describe(this_susp.mat)
 
    [this_susp_per_sig,this_susp_per_sig_in_gene_set] = run_single_susp(this_brew,this_susp,susp_cutoff,cells{ii},'mkfig',false);
    
    figure(2)
    [~,~,x,f] = kde(this_susp_per_sig,2^10,0,200);
    plot(ax2,x,f,'DisplayName',cells{ii});
    hold on
    
    figure(3)
    [~,~,x,f] = kde(this_susp_per_sig_in_gene_set,2^10,0,100);
    plot(ax3,x,f,'DisplayName',cells{ii});
    hold on
    
end

num_sig = length(this_susp_per_sig);

figure(1)
xlabel('Suspisciousness')
ylabel('cdf')
title(sprintf('CDF of suspect z-scores \n num_rep = %d',num_reps),...
    'interpreter','none')
ax = gca;
ax.YLim = [0 1];
ax.YTick = [.8 .85 .9:.025:1];
grid on
legend show
namefig('REP_susp_cdf')

figure(2)
title(sprintf('CDF of suspect z-scores per signature \n cutoff = %d \n num_rep = %d',susp_cutoff,num_reps),...
    'interpreter','none')
grid on
xlabel('Suspect genes per sample')
ylabel('cdf')
legend show
ax = gca;
ax.YTick = [0 .2 .4 .5 .6 .8 1];
namefig('REP_susp_per_sig_cdf')

figure(3)
title(sprintf('CDF of suspect z-scores per signature in gene set \n cutoff = %d \n num_rep = %d',susp_cutoff,num_reps),...
    'interpreter','none')
xlabel('Suspect genes per sample in gene set')
ylabel('cdf')
grid on
legend show
ax = gca;
ax.YTick = [0 .2 .4 .5 .6 .8 1];
namefig('REP_susp_per_sig_gene_set_cdf')


%savefigures('mkdir',false,'out','/cmap/users/jnasser/dpeak_v2/fig/A2_suspect','overwrite',true)

% legend show
% grid on
% title('A2 CDF of suspect Z-scores. Compounds with 3 replicates')
% ax = gca;
% ax.YTick = [0:.1:.9 .91:.01:1]