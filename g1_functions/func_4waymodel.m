function [rm_full, table_normality] = func_4waymodel(data_eo, data_ec, specs)

% v1: scaling range (lo vs. hi)
% v2: physiological state (EO vs. EC)
% v3: EEG channel list (e.g., FC2 vs. POz)
% v4: spatial filter list (e.g., car vs. lap)

% get model specifics
vnames = specs.vnames;
v1_list = specs.v1_list(:);
v2_list = specs.v2_list(:);
v3_list = specs.v3_list(:);
v4_list = specs.v4_list(:);

n1 = length(v1_list);
n2 = length(v2_list);
n3 = length(v3_list);
n4 = length(v4_list);

% sort variables
table_full = table();
N = n1*n2*n3*n4;
v1_dummy = cell(N,1);
v2_dummy = cell(N,1);
v3_dummy = cell(N,1);
v4_dummy = cell(N,1);
ind = 1;

% eyes open
for v1 = 1:n1
    for v2 = 1:n2
        for v3 = 1:n3
            for v4 = 1:n4
                v1i = v1_list{v1};
                v2i = v2_list{v2};
                v3i = v3_list{v3};
                v4i = v4_list{v4};
                if contains('EO',{v1i,v2i,v3i,v4i})
                    table_full.([v4i '_' v1i '_' v2i '_' v3i]) = data_eo.([v1i '_' v3i]).(v4i);
                elseif contains('EC',{v1i,v2i,v3i,v4i})
                    table_full.([v4i '_' v1i '_' v2i '_' v3i]) = data_ec.([v1i '_' v3i]).(v4i);
                end
                
                v1_dummy{ind} = v1i;
                v2_dummy{ind} = v2i;
                v3_dummy{ind} = v3i;
                v4_dummy{ind} = v4i;
                ind = ind + 1;
            end
        end
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

% 4-way repeated measures model
WithinDesign = table(...
    v4_dummy, v3_dummy, v2_dummy, v1_dummy,...
    'VariableNames', vnames);

WithinDesign.(vnames{4}) = categorical(WithinDesign.(vnames{4}));
WithinDesign.(vnames{3}) = categorical(WithinDesign.(vnames{3}));
WithinDesign.(vnames{2}) = categorical(WithinDesign.(vnames{2}));
WithinDesign.(vnames{1}) = categorical(WithinDesign.(vnames{1}));

rm_exp = [v4_list{1} '_' v1_list{1} '_' v2_list{1} '_' v3_list{1} '-' v4_list{end} '_' v1_list{end} '_' v2_list{end} '_' v3_list{end} ' ~ 1'];
rm_full_model = fitrm(table_full, rm_exp, 'WithinDesign', WithinDesign);
rm_full = ranova(rm_full_model, 'WithinModel', [vnames{1} '*' vnames{2} '*' vnames{3} '*' vnames{4}]);

end