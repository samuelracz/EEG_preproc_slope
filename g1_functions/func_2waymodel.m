function [rm_full, table_normality, table_pairwise] = func_2waymodel(data_eo, data_ec, specs)

% v1: scaling range (lo vs. hi)
% v2: physiological state (EO vs. EC)
% v3: EEG channel (e.g., FC2)
% v4: spatial filter (e.g., car)

% get model specifics
vnames = specs.vnames;
v1_list = specs.v1_list(:);
v2_list = specs.v2_list(:);
v3 = specs.v3{1};
v4 = specs.v4{1};

n1 = length(v1_list);
n2 = length(v2_list);

% sort variables
table_full = table();
N = n1*n2;
v1_dummy = cell(N,1);
v2_dummy = cell(N,1);
ind = 1;

% eyes open
for v1 = 1:n1
    for v2 = 1:n2
        v1i = v1_list{v1};
        v2i = v2_list{v2};
        if contains('EO',{v1i,v2i})
            table_full.([v4 '_' v1i '_' v2i '_' v3]) = data_eo.([v1i '_' v3]).(v4);
        elseif contains('EC',{v1i,v2i})
            table_full.([v4 '_' v1i '_' v2i '_' v3]) = data_ec.([v1i '_' v3]).(v4);
        end

        v1_dummy{ind} = v1i;
        v2_dummy{ind} = v2i;
        ind = ind + 1;
    end
end

% test for normality
vn = table_full.Properties.VariableNames;
nv = length(vn);

table_normality = table(cell(nv,1),zeros(nv,1),zeros(nv,1),...
    'VariableNames',{'var','h','p'});

for n = 1:nv
    [h,p] = lillietest(table_full.(vn{n}),'Alpha',0.05/N);
    table_normality.var{n} = vn{n};
    table_normality.h(n) = h;
    table_normality.p(n) = p;
end

% 2-way repeated measures model
WithinDesign = table(...
    v2_dummy, v1_dummy,...
    'VariableNames', vnames);

WithinDesign.(vnames{2}) = categorical(WithinDesign.(vnames{2}));
WithinDesign.(vnames{1}) = categorical(WithinDesign.(vnames{1}));

rm_exp = [v4 '_' v1_list{1} '_' v2_list{1} '_' v3 '-' v4 '_' v1_list{end} '_' v2_list{end} '_' v3 ' ~ 1'];
rm_full_model = fitrm(table_full, rm_exp, 'WithinDesign', WithinDesign);
rm_full = ranova(rm_full_model, 'WithinModel', [vnames{1} '*' vnames{2}]);

% pairwise comparisons
table_pairwise = table(...
    cell(N,1),...
    cell(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    zeros(N,1),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    'VariableNames',{'var1','var2','E1','E2','pc','hc','t','df','CohenD','v1','v2','ttest','ES'});

comp_pairs = [1,3; 2,4; 1,2; 3,4];


for n = 1:N
    i1 = comp_pairs(n,1);
    i2 = comp_pairs(n,2);
    vname1 = vn{i1};
    vname2 = vn{i2};
    var1 = table_full.(vname1);
    var2 = table_full.(vname2);

    [h,p,~,stat] = ttest(var1,var2);
    Cd = meanEffectSize(var1,var2,Effect='cohen');

    table_pairwise.var1{n} = vname1;
    table_pairwise.var2{n} = vname2;
    table_pairwise.E1(n) = mean(var1);
    table_pairwise.E2(n) = mean(var2);
    table_pairwise.pc(n) = p*N;
    table_pairwise.hc(n) = p*N < 0.05;
    table_pairwise.t(n) = stat.tstat;
    table_pairwise.df(n) = stat.df;
    table_pairwise.CohenD(n) = abs(Cd.Effect('CohensD'));
    table_pairwise.v1{n} = var1;
    table_pairwise.v2{n} = var2;
    table_pairwise.ttest{n} = stat;
    table_pairwise.ES{n} = Cd;
end

end