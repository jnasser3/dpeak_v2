%%
addpath('/cmap/projects/repurposing/jnwork/code/')

%% Roast the plates
roast('plate','/cmap/users/jnasser/dpeak_v2/small_test/plate_list.grp', ...
  'plate_path','/cmap/users/jnasser/dpeak_v2/small_test', ...
  'raw_path','/cmap/obelix/pod/custom/REP/lxb', ...
  'map_path','/cmap/obelix/pod/custom/REP/maps')

%% Brew the plates
brew('plate','/cmap/users/jnasser/dpeak_v2/small_test/plate_list_B.grp',...
        'plate_path','/cmap/users/jnasser/dpeak_v2/small_test', ...
        'brew_path','/cmap/users/jnasser/dpeak_v2/small_test', ...
        'group_by', 'rna_well', ...
        'zmad_ref', 'ZSPC', ...
        'filter_vehicle', false)

%% Compute dpeak suspiciousness
brew_ds = parse_gctx('/cmap/users/jnasser/dpeak_v2/small_test/REP.A010_A375_24H/by_rna_well/REP.A010_A375_24H_COMPZ.MODZ_SCORE_LM_n375x978.gctx');
roast_ds = parse_gctx('/cmap/users/jnasser/dpeak_v2/small_test/REP.A010_A375_24H/REP.A010_A375_24H_ZSPCQNORM_n1087x978.gctx');
suspA = compute_misscall_suspects_per_brew_gene(...
    brew_ds,roast_ds);

% %% Compute dpeak suspiciousness
% brew_ds = parse_gctx('/cmap/users/jnasser/dpeak_v2/small_test/REP.B010_A375_24H/by_rna_well/REP.B010_A375_24H_COMPZ.MODZ_SCORE_LM_n375x978.gctx');
% roast_ds = parse_gctx('/cmap/users/jnasser/dpeak_v2/small_test/REP.B010_A375_24H/REP.B010_A375_24H_ZSPCQNORM_n1087x978.gctx');
% suspB = compute_misscall_suspects_per_brew_gene(...
%     brew_dsB,roast_dsB);