function compare_susp_mat(susp1,susp2)

%% Params
susp_cutoff = 10

%% setup
common_ids = intersect(susp1.cid,susp2.cid);
susp1 = ds_slice(susp1,'cid',common_ids);
susp2 = ds_slice(susp2,'cid',common_ids);

assert(isequal(susp1.rid,susp2.rid))
assert(isequal(susp1.cid,susp2.cid))

%% Overall cdf of suspiciousness
[e1,x1] = ecdf(susp1.mat(:));
[e2,x2] = ecdf(susp2.mat(:));

figure;
plot(x1,e1,'DisplayName','Susp1')
hold on
plot(x2,e2,'DisplayName','Susp2')
title('Overall cdf of suspisciousness')
grid on
legend show
xlabel('Zscore suspiciousness')

%% Per gene
figure;
data1 = sum(susp1.mat > susp_cutoff,2);
data2 = sum(susp2.mat > susp_cutoff,2);
scatter(data1,data2)
grid on

% figure
% plot(sum(susp1.mat > 10,2))
% hold on
% plot(sum(susp2.mat > 10,2))

%% How many entries in the two susp matrices are aligned
figure;
idx1 = (susp1.mat > susp_cutoff);
idx2 = (susp2.mat > susp_cutoff);
idx = idx1 | idx2;
scatter(susp1.mat(idx),susp2.mat(idx))
xlabel('Susp1')
ylabel('Susp2')
title_str = sprintf('Suspect zscore cutoff = %d \n number susp1 = %d \n number susp2 = %d \n number in common = %d',...
    susp_cutoff,nnz(idx1),nnz(idx2),nnz(idx1 & idx2));
title(title_str);

end