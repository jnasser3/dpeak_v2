

for temp = 43:10:300
    rng(1)
    temp
    x = lxb.RP1(lxb.RID == temp);
    detect_lxb_peaks_single(x,'pkmethod','kmeans_jn','showfig',true);
    detect_lxb_peaks_single(x,'pkmethod','kmeans_viable','showfig',true);
end

%%
temp = 72
%rng(1)
lxb = parse_lxb('/Volumes/cmap_obelix/pod/custom/REP/lxb/REP.B006_A375_24H_X1_B22/REP.B006_A375_24H_X1_B22_B20.lxb',...
    'allfields',false);
x = lxb.RP1(lxb.RID == temp);
t1 = detect_lxb_peaks_single(x,'pkmethod','kmeans_jn','showfig',true,...
    'remove_outlier',true,...
    'ksmethod','kspeak');
%t2 = detect_lxb_peaks_single(x,'pkmethod','kmeans_viable','showfig',true);
% colors = 'kr';
% coloridx = ((lxb.RID == temp) + 1);
% plot3(lxb.CL1,lxb.CL2,lxb.CL3);


%% Find instances in pkstats where only 1 peak called
load /cmap/users/jnasser/dpeak_v2/small_test/REP.A006_A375_24H_X2_B22/dpeak/REP.A006_A375_24H_X2_B22.mat
s = cellfun(@(c) length(c), {pkstats(:).pksupport});
idx = find(s == 1);

counter = 0;
for ii = 1:length(idx)
    if ~ismember(pkstats(idx(ii)).analyte,[1:11,499])
        counter = counter + 1;
        idx(ii)
        pkstats(idx(ii))
    end
end

%% GEX Level
% idx = 'Analyte 151';
% gex = parse_gctx('/cmap/users/jnasser/dpeak_v2/small_test/REP.A006_A375_24H_X2_B22/REP.A006_A375_24H_X2_B22_GEX_n384x978.gct');
% data = ds_subset(gex,'row','pr_analyte_id',idx);
% data.mat = safelog2(data.mat);
% figure
% scatter(data.mat(1,:),data.mat(2,:))
% xlabel(data.rdesc(1,data.rdict('pr_gene_symbol')))
% ylabel(data.rdesc(2,data.rdict('pr_gene_symbol')))

%% GEX Level Scatter
%load /cmap/users/jnasser/dpeak_v2/small_test3/REP.B006_A375_24H_X2_B22/dpeak/REP.B006_A375_24H_X2_B22.mat

for ii = 72%39:20:500
figure('position', [300, 300, 800, 400]) 
subplot(1,4,[1 2]);
analyte_idx = ii;
has2_idx = (cellfun(@(x) length(x), {pkstats(analyte_idx,:).pkexp}) == 2);
data = safelog2(vertcat(pkstats(analyte_idx,has2_idx).pkexp));
conf = vertcat(pkstats(analyte_idx,has2_idx).conf2);
high_conf_idx = (conf < 5);
scatter(data(:,1),data(:,2),[],~high_conf_idx)
colormap(cool)
axis([0 16 0 16])
hold on 
plot([0 16], [0 16], 'r:')
title(sprintf('Analyte: %s - Preflip \n Purple = low confidence calls, blue = high confidence calls \n num_high_conf = %d, num_low_conf = %d', num2str(ii), nnz(high_conf_idx),nnz(~high_conf_idx)),'Interpreter','none')
ax = gca;
ax.XTick = [0:2:16];
ax.YTick = [0:2:16];
grid on

%plot gaussian contour
GM = fitgmdist(data(high_conf_idx,:),1);
hold on
x = 0:.1:16;
[X1,X2] = meshgrid(x,x);
F = mvnpdf([X1(:) X2(:)],GM.mu,GM.Sigma);
F = reshape(F,length(x),length(x));
contour(x,x,F);

%Run flip adjust v2
subplot(1,4,[3 4])
[data2,flips] = flip_correction_v2(data,conf);
flip_idx = ismember(1:length(data2),flips)';
scatter(data2(:,1),data2(:,2),[],flip_idx)
colormap(cool)
axis([0 16 0 16])
hold on 
plot([0 16], [0 16], 'r:')
title(sprintf('Analyte: %s - Postflip \n Purple = flipped, blue = not flipped \n num_flipped = %d', num2str(ii), length(flips)),'Interpreter','none')
ax = gca;
ax.XTick = [0:2:16];
ax.YTick = [0:2:16];
grid on
end

%% Investigate OLD dpeak data
lxb = parse_lxb('/Volumes/cmap_obelix/lxb/ERG005_VCAP_24H_X1_B1_UNI5253R/ERG005_VCAP_24H_X1_B1_UNI5253R_C20.lxb');
for analyte_num = 18%41:70;
    x = lxb.RP1(lxb.RID == analyte_num);
    %t1 = detect_lxb_peaks_single(x,'pkmethod','kmeans_jn','showfig',true);
    t2 = detect_lxb_peaks_single(x,'ksmethod','median','showfig',true);
end
