function [susp1,susp2,brew1,brew2] = compare_brews_wrt_dpeak(brewdir1,brewdir2)

%% Load in matrices
roast1 = parse_gctx(strtrim(ls(fullfile(brewdir1,'*ZSPCQNORM*.gctx'))));
roast2 = parse_gctx(strtrim(ls(fullfile(brewdir2,'*ZSPCQNORM*.gctx'))));
brew1 = parse_gctx(strtrim(ls(fullfile(brewdir1,'by_rna_well','*MODZ_SCORE_LM*.gctx'))));
brew2 = parse_gctx(strtrim(ls(fullfile(brewdir2,'by_rna_well','*MODZ_SCORE_LM*.gctx'))));

assert(isequal(brew1.rid,brew2.rid))
assert(isequal(roast1.rid,roast2.rid))
%% Subset to common ids
common_roast_ids = intersect(roast1.cid,roast2.cid);
common_brew_ids = intersect(brew1.cid,brew2.cid);

roast1 = ds_slice(roast1,'cid',common_roast_ids);
roast2 = ds_slice(roast2,'cid',common_roast_ids);
brew1 = ds_slice(brew1,'cid',common_brew_ids);
brew2 = ds_slice(brew2,'cid',common_brew_ids);

%% Suspect Matrices
susp1 = compute_misscall_suspects_per_brew_gene(brew1,roast1);
susp2 = compute_misscall_suspects_per_brew_gene(brew2,roast2);

%% Plots
compare_susp_mat(susp1,susp2)

%% Signature strength of DMSO
figure
dmso1 = ds_subset(brew1,'column','pert_iname','DMSO');
dmso2 = ds_subset(brew2,'column','pert_iname','DMSO');
assert(isequal(dmso1.cid,dmso2.cid))

scatter(cell2mat(dmso1.cdesc(:,dmso1.cdict('distil_ss'))),...
    cell2mat(dmso2.cdesc(:,dmso2.cdict('distil_ss'))))
ax = gca;
max_lim = max(ax.YLim(2),ax.XLim(2));
min_lim = max(ax.YLim(1),ax.XLim(1));
axis([min_lim max_lim min_lim max_lim])
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:')
grid on
title('Sig Strength of DMSO')
xlabel('Data1')
ylabel('Data2')

%% S-C plots comparison
figure
idx1 = (cell2mat(brew1.cdesc(:,brew1.cdict('distil_cc_q75'))) > 0);
idx2 = (cell2mat(brew2.cdesc(:,brew2.cdict('distil_cc_q75'))) > 0);
idx = idx1 & idx2;

sc1 = ds_slice(brew1,'cidx',idx);
sc2 = ds_slice(brew2,'cidx',idx);
assert(isequal(sc1.cid,sc2.cid));

rr1 = cell2mat(sc1.cdesc(:,sc1.cdict('distil_cc_q75')));
rr2 = cell2mat(sc2.cdesc(:,sc2.cdict('distil_cc_q75')));
ss1 = cell2mat(sc1.cdesc(:,sc1.cdict('distil_ss')));
ss2 = cell2mat(sc2.cdesc(:,sc2.cdict('distil_ss')));

scatter(rr1,ss1,[],'b')
hold on
scatter(rr2,ss2,[],'r')

end

