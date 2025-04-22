function [f, corr_struct] = func_plot_corr_EO_EC_sa(data_eo, data_ec, filt_list, ch_list, atype, mc_alpha)

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
    'r1', [],...
    'p1', [],...
    'c1', [],...
    'r2', [],...
    'p2', [],...
    'c2', []);


if length(filt_list) == 2
    f = figure('color','w','units','normalized','outerposition',[0.15 0 0.7 0.95]);
elseif length(filt_list) == 3
    f = figure('color','w','units','normalized','outerposition',[0 0 1 0.95]);
end

for ch = 1:nch
    for ftype = 1:nf
        tmp_dsp_lo = (data_ec.(['lo_' ch_list{ch}]).(filt_list{ftype}) - data_eo.(['lo_' ch_list{ch}]).(filt_list{ftype}));
        tmp_dsp_hi = (data_ec.(['hi_' ch_list{ch}]).(filt_list{ftype}) - data_eo.(['hi_' ch_list{ch}]).(filt_list{ftype}));
        tmp_dlao = (data_ec.([atype '_' ch_list{ch}]).(filt_list{ftype}) - data_eo.([atype '_' ch_list{ch}]).(filt_list{ftype}));

        params.ch = ch;
        params.ftype = ftype;

        subplot(nch,nf,ind_sp)
        corr_struct_tmp = func_plot_both_correlations(tmp_dlao, tmp_dsp_lo, tmp_dsp_hi, params, mc_alpha);
        corr_struct(ind_sp).ctype = [ch_list{ch} '_' filt_list{ftype} '_' atype];
        corr_struct(ind_sp).r1 = corr_struct_tmp.r1;
        corr_struct(ind_sp).p1 = corr_struct_tmp.p1;
        corr_struct(ind_sp).c1 = corr_struct_tmp.c1;
        corr_struct(ind_sp).r2 = corr_struct_tmp.r2;
        corr_struct(ind_sp).p2 = corr_struct_tmp.p2;
        corr_struct(ind_sp).c2 = corr_struct_tmp.c2;
        ind_sp = ind_sp + 1;
    end
end


end

%% helper functions
function [corr_struct] = func_plot_both_correlations(data_alpha, data_lo, data_hi, params, mc_alpha)

if strcmp(params.atype,'ar')
    xlab = '\DeltaBLP^{\alpha}_{raw} (EC-EO)';
elseif strcmp(params.atype,'ao')
    xlab = '\DeltaBLP^{\alpha}_{osci} (EC-EO)';
elseif strcmp(params.atype,'lar')
    xlab = '\DeltaLog(BLP^{\alpha}_{raw}) (EC-EO)';
elseif strcmp(params.atype,'lao')
    xlab = '\DeltaLog(BLP^{\alpha}_{osci}) (EC-EO)';
end

chlab = params.ch_list{params.ch};
ftype = params.filt_list{params.ftype};

plo = plot(data_alpha, data_lo, 'bo','MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', 12, 'LineWidth', 2);
hold on
phi = plot(data_alpha, data_hi, 'ro','MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 2);
ls = lsline;
ls = fliplr(ls);

% low-range beta
if lillietest(data_alpha) || lillietest(data_lo)
    [r1,p1] = corr(data_alpha, data_lo, 'type','spearman');
    ct1 = 'spearman';
else
    [r1,p1] = corr(data_alpha, data_lo, 'type','pearson');
    ct1 = 'pearson';
end

p1 = p1*mc_alpha;

% multiple comparisons adjustment
if p1<0.05
    leg_r1 = ['{\itr}=' num2str(r1,'%.4f')];
    leg_p1 =['{\itp}=' num2str(p1,'%.4f')];
    ls(1).LineWidth = 3;
else
    leg_r1 = ['{\itr}=' num2str(r1,'%.4f')];
    leg_p1 = ['{\itp}=' num2str(p1,'%.4f') ' (ns.)'];
    ls(1).LineWidth = 1.5;
    ls(1).LineStyle = ':';
end

% high-range beta
if lillietest(data_alpha) || lillietest(data_hi)
    [r2,p2] = corr(data_alpha, data_hi, 'type','spearman');
    ct2 = 'spearman';
else
    [r2,p2] = corr(data_alpha, data_hi, 'type','pearson');
    ct2 = 'pearson';
end

p2 = p2*mc_alpha;

% multiple comparisons adjustment
if p2<0.05
    leg_r2 = ['{\itr}=' num2str(r2,'%.4f')];
    leg_p2 =['{\itp}=' num2str(p2,'%.4f')];
    ls(2).LineWidth = 3;
else
    leg_r2 = ['{\itr}=' num2str(r2,'%.4f')];
    leg_p2 = ['{\itp}=' num2str(p2,'%.4f') ' (ns.)'];
    ls(2).LineWidth = 1.5;
    ls(2).LineStyle = ':';
end

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
legend(ls,{'\Delta\beta_{lo}', '\Delta\beta_{hi}'},...
    'FontSize',24,'location','northwest')

corr_struct = struct(...
    'r1', r1,...
    'r2', r2,...
    'c1', ct1,...
    'p1', p1,...
    'p2', p2,...
    'c2', ct2);
end