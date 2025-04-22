function [f, corr_struct] = func_plot_corr_EO_EC_bb(data_eo, data_ec, filt_list, ch_list, atype, mc_alpha)

% NOTE: using 'lao' for log-transformed oscillatory alpha power 

if ~strcmp(atype,'ar') && ~strcmp(atype,'ao') && ~strcmp(atype,'lar') && ~strcmp(atype,'lao')
    error('Wrong input for alpha power: use raw mixd alpha [ar], raw oscillatory alpha [ao], log-transformed mixd alpha [lar] or log-transformed oscillatory alpha [lao].')
end

params = struct(...
    'filt_list', {filt_list},...
    'ch_list', {ch_list},...
    'atype', atype,...
    'ch', [],...
    'ftype', []);

nf = length(filt_list);
nch = length(ch_list);
ind_sp = 1;

corr_struct = struct(...
    'ctype', [],...
    'r', [],...
    'p', [],...
    'c', []);


if length(filt_list) == 2
    f = figure('color','w','units','normalized','outerposition',[0.15 0 0.7 0.95]);
elseif length(filt_list) == 3
    f = figure('color','w','units','normalized','outerposition',[0 0 1 0.95]);
end

for ch = 1:nch
    for ftype = 1:nf
        tmp_dsp_bb = (data_ec.(['bb_' ch_list{ch}]).(filt_list{ftype}) - data_eo.(['bb_' ch_list{ch}]).(filt_list{ftype}));
        tmp_dlao = (data_ec.([atype '_' ch_list{ch}]).(filt_list{ftype}) - data_eo.([atype '_' ch_list{ch}]).(filt_list{ftype}));

        params.ch = ch_list{ch};
        params.ftype = filt_list{ftype};

        subplot(nch,nf,ind_sp)
        corr_struct_tmp = func_plot_correlations(tmp_dlao, tmp_dsp_bb, params, mc_alpha);
        corr_struct(ind_sp).ctype = [ch_list{ch} '_' filt_list{ftype} '_' atype];
        corr_struct(ind_sp).r = corr_struct_tmp.r;
        corr_struct(ind_sp).p = corr_struct_tmp.p;
        corr_struct(ind_sp).c = corr_struct_tmp.c;
        ind_sp = ind_sp + 1;
    end
end


end

%% helper functions
function [corr_struct] = func_plot_correlations(data_alpha, data_bb, params, mc_alpha)

if strcmp(params.atype,'ar')
    xlab = '\DeltaBLP^{\alpha}_{raw} (EC-EO)';
elseif strcmp(params.atype,'ao')
    xlab = '\DeltaBLP^{\alpha}_{osci} (EC-EO)';
elseif strcmp(params.atype,'lar')
    xlab = '\DeltaLog(BLP^{\alpha}_{raw}) (EC-EO)';
elseif strcmp(params.atype,'lao')
    xlab = '\DeltaLog(BLP^{\alpha}_{osci}) (EC-EO)';
end

chlab = params.ch;
ftype = params.ftype;

pbb = plot(data_alpha, data_bb, 'o','color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12, 'LineWidth', 2);
hold on
ls = lsline;
ls = fliplr(ls);

% broadband beta
if lillietest(data_alpha) || lillietest(data_bb)
    [r,p] = corr(data_alpha, data_bb, 'type','spearman');
    ct = 'spearman';
else
    [r,p] = corr(data_alpha, data_bb, 'type','pearson');
    ct = 'pearson';
end

p = p*mc_alpha;

% multiple comparisons adjustment
if p<0.05
    ls(1).LineWidth = 3;
else
    ls(1).LineWidth = 1.5;
    ls(1).LineStyle = ':';
end

% plot specifics
set(gca,'FontSize',24,'LineWidth',1,'box','on')
xlabel(xlab,'FontSize',28)
ylabel('\Delta\beta (EC-EO)','FontSize',28)

if strcmp(ftype,'raw')
    str_filt = 'Raw EEG';
elseif strcmp(ftype,'car')
    str_filt = 'Common Average Reference';
elseif strcmp(ftype,'lap')
    str_filt = 'Small Laplacian Filtering';
end

title({[str_filt ' - ' chlab], '\DeltaPower^{\alpha} vs. \Delta\beta'},'FontSize',32)
legend(ls,{'\Delta\beta_{bb}'},...
    'FontSize',24,'location','northwest')

corr_struct = struct(...
    'r', r,...
    'p', p,...
    'c', ct);
end
