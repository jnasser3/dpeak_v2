%Compute suspect zscores on A2

cells = {'A375','A549','HT29','HEPG2','HCC515','HA1E','MCF7','PC3','VCAP'};

susp_cutoff = 5;

close all
for ii = 1:length(cells)
    
    this_susp = parse_gctx(strtrim(ls(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices/',sprintf('%s*SUSP*.gctx',cells{ii})))));
    this_brew = parse_gctx(strtrim(ls(fullfile('/cmap/users/jnasser/dpeak_v2/data/A2_suspect_matrices/',sprintf('%s*BREW*.gctx',cells{ii})))));
    
    %[f,x] = ecdf(this_susp.mat(:));
%     [~,~,x,f] = kde(this_susp.mat(:),2^10,0,50);
%     plot(x,f,'DisplayName',cells{ii});
%     hold on
%    
%     cells(ii)
%     describe(this_susp.mat)

    run_single_susp(this_brew,this_susp,susp_cutoff,cells{ii})
    
end

savefigures('mkdir',false,'out','/cmap/users/jnasser/dpeak_v2/fig/','overwrite',true)

% legend show
% grid on
% title('A2 CDF of suspect Z-scores. Compounds with 3 replicates')
% ax = gca;
% ax.YTick = [0:.1:.9 .91:.01:1]