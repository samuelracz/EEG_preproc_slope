function [ f, rm_cell] = func_plot_EO_EC_bbslopes(data_eo, data_ec, filt_list, ch_list)

nf = length(filt_list);
nch = length(ch_list);
ind_sp = 1;

rm_cell = cell(nch*nf,1);

if length(filt_list) == 2
    f = figure('color','w','units','normalized','outerposition',[0.15 0.2 0.7 0.55]);
elseif length(filt_list) == 3
    f = figure('color','w','units','normalized','outerposition',[0 0 1 0.95]);
end

n_multcomp = length(ch_list);

for ftype = 1:nf
    subplot(1,2,ftype)

    data_table = table();
    for ch = 1:nch
        tmp_bb_eo = data_eo.(['bb_' ch_list{ch}]).(filt_list{ftype});
        tmp_bb_ec = data_ec.(['bb_' ch_list{ch}]).(filt_list{ftype});

        data_table.([ch_list{ch} '_bb_eo']) = tmp_bb_eo;
        data_table.([ch_list{ch} '_bb_ec']) = tmp_bb_ec;

        N = length(tmp_bb_eo);

        % create boxplots
        plot((2+(ch-1)*10)*ones(N,1),tmp_bb_eo,'r*','MarkerSize',12)
        hold on
        plot((5+(ch-1)*10)*ones(N,1),tmp_bb_ec,'b*','MarkerSize',12)

        % draw boxes
        pm_eo = fill([1+(ch-1)*10 3+(ch-1)*10 3+(ch-1)*10 1+(ch-1)*10],[quantile(tmp_bb_eo,0.25), quantile(tmp_bb_eo,0.25), quantile(tmp_bb_eo,0.75) quantile(tmp_bb_eo,0.75)],[1 0.5 0.5],'LineWidth',1.5,'FaceAlpha',0.4);
        pm_ec = fill([4+(ch-1)*10 6+(ch-1)*10 6+(ch-1)*10 4+(ch-1)*10],[quantile(tmp_bb_ec,0.25), quantile(tmp_bb_ec,0.25), quantile(tmp_bb_ec,0.75) quantile(tmp_bb_ec,0.75)],[0.5 0.5 1],'LineWidth',1.5,'FaceAlpha',0.4);

        % draw percentiles
        plot([2+(ch-1)*10 2+(ch-1)*10],[quantile(tmp_bb_eo,0.75) quantile(tmp_bb_eo,0.95)],'k-','LineWidth',1.5)
        plot([1.5+(ch-1)*10 2.5+(ch-1)*10],[quantile(tmp_bb_eo,0.95) quantile(tmp_bb_eo,0.95)],'k-','LineWidth',1.5)
        plot([2+(ch-1)*10 2+(ch-1)*10],[quantile(tmp_bb_eo,0.25) quantile(tmp_bb_eo,0.05)],'k-','LineWidth',1.5)
        plot([1.5+(ch-1)*10 2.5+(ch-1)*10],[quantile(tmp_bb_eo,0.05) quantile(tmp_bb_eo,0.05)],'k-','LineWidth',1.5)

        plot([5+(ch-1)*10 5+(ch-1)*10],[quantile(tmp_bb_ec,0.75) quantile(tmp_bb_ec,0.95)],'k-','LineWidth',1.5)
        plot([4.5+(ch-1)*10 5.5+(ch-1)*10],[quantile(tmp_bb_ec,0.95) quantile(tmp_bb_ec,0.95)],'k-','LineWidth',1.5)
        plot([5+(ch-1)*10 5+(ch-1)*10],[quantile(tmp_bb_ec,0.25) quantile(tmp_bb_ec,0.05)],'k-','LineWidth',1.5)
        plot([4.5+(ch-1)*10 5.5+(ch-1)*10],[quantile(tmp_bb_ec,0.05) quantile(tmp_bb_ec,0.05)],'k-','LineWidth',1.5)

        % plot means and medians
        plot([1+(ch-1)*10 3+(ch-1)*10],[mean(tmp_bb_eo) mean(tmp_bb_eo)],'r-', 'LineWidth',3)
        plot([4+(ch-1)*10 6+(ch-1)*10],[mean(tmp_bb_ec) mean(tmp_bb_ec)],'b-', 'LineWidth',3)

        plot([1+(ch-1)*10 3+(ch-1)*10],[median(tmp_bb_eo) median(tmp_bb_eo)],'r--', 'LineWidth',3)
        plot([4+(ch-1)*10 6+(ch-1)*10],[median(tmp_bb_ec) median(tmp_bb_ec)],'b--', 'LineWidth',3)

        % plot trends
        if lillietest(tmp_bb_eo,'Alpha',0.0125) || lillietest(tmp_bb_ec,'Alpha',0.0125)
            [p,~] = signrank(tmp_bb_eo,tmp_bb_ec);
            p = p*2;
            if p<0.05
                plot([2+(ch-1)*10 5+(ch-1)*10],[median(tmp_bb_eo) median(tmp_bb_ec)],'ko-', 'LineWidth', 3)
            else
                plot([2+(ch-1)*10 5+(ch_list-1)*10],[median(tmp_bb_eo) median(tmp_bb_ec)],'ko:', 'LineWidth', 3)
            end
        else
            [~,p] = ttest(tmp_bb_eo,tmp_bb_ec);
            p = p*2;
            if p<0.05
                plot([2+(ch-1)*10 5+(ch-1)*10],[mean(tmp_bb_eo) mean(tmp_bb_ec)],'ko-', 'LineWidth', 3)
            else
                plot([2+(ch-1)*10 5+(ch-1)*10],[mean(tmp_bb_eo) mean(tmp_bb_ec)],'ko:', 'LineWidth', 3)
            end
        end

    end

    % get ranges
    ymin = min(get(gca,'ylim'));
    ymax = max(get(gca,'ylim'));
    yrange = ymax - ymin;

    ylim([ymin-0.1*yrange ymax+0.1*yrange])
    

    % WithinDesign = table(repmat({'EO';'EC'},nch,1), {ch_list{1}; ch_list{2}; ch_list{1}; ch_list{2}}, 'VariableNames', {'state','ch'});
    WithinDesign = table(repmat({'EO';'EC'},nch,1), repmat(ch_list',2,1), 'VariableNames', {'state','ch'});
    WithinDesign.state = categorical(WithinDesign.state);
    WithinDesign.ch = categorical(WithinDesign.ch);

    rm_model = fitrm(data_table, [ch_list{1} '_bb_eo-' ch_list{end} '_bb_ec ~ 1'],'WithinDesign',WithinDesign);
    rm = ranova(rm_model,'WithinModel','state*ch');

    rm_cell{ind_sp} = rm;

    % statistical tests - EO vs EC
    for ch = 1:nch
        if lillietest(data_table.([ch_list{ch} '_bb_eo']),'Alpha',0.05/n_multcomp) || lillietest(data_table.([ch_list{ch} '_bb_ec']),'Alpha',0.05/n_multcomp)
            [p,~] = signrank(data_table.([ch_list{ch} '_bb_eo']),data_table.([ch_list{ch} '_bb_ec']));
        else
            [~,p] = ttest(data_table.([ch_list{ch} '_bb_eo']),data_table.([ch_list{ch} '_bb_ec']));
        end

        p = p*n_multcomp;

        if p<0.05
            plot([(ch-1)*10+1 (ch-1)*10+6],[ymax ymax],'k-','LineWidth',4)
            if p>=0.0001
                text((ch-1)*10+3.5,ymax+0.1*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text((ch-1)*10+3.5,ymax+0.1*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end
    end

    % plot specifics
    ymin = min(get(gca,'ylim'));
    ymax = max(get(gca,'ylim'));
    yrange = ymax - ymin;
    ylim([ymin-0.3*yrange ymax+0.4*yrange])
    xlim([0 (nch-1)*10+7])
    if nch == 2
        set(gca,'FontSize',20,'LineWidth',2,'XTick',[3.5 13.5],'XTickLabel',ch_list)
        xlabel(['\beta_E_O vs. \beta_E_C over ' ch_list{1} ' and ' ch_list{2}],'FontSize',24)
    elseif nch == 3
        set(gca,'FontSize',20,'LineWidth',2,'XTick',[3.5 13.5 23.5],'XTickLabel',ch_list)
        xlabel(['\beta_E_O vs. \beta_E_C over ' ch_list{1} ', ' ch_list{2} ' and ' ch_list{3}],'FontSize',24)
    end
    ylabel('Spectral Slope','FontSize',24)
    legend([pm_eo,pm_ec],{'EO','EC'},'location','southeast','FontSize',24)
    p = rm.pValueGG('(Intercept):state:ch');
    if p < 0.0001
        p_str = '{\itp}<0.0001';
    else
        p_str = ['{\itp}=' num2str(p,'%.4f')];
    end

    if strcmp(filt_list{ftype},'raw')
        str_filt = 'Raw EEG';
    elseif strcmp(filt_list{ftype},'car')
        str_filt = 'Common Average Reference';
    elseif strcmp(filt_list{ftype},'lap')
        str_filt = 'Small Laplacian Filtering';
    end
    
    if nch == 2
        title({[str_filt ' - ' ch_list{1} ' and ' ch_list{2}],['state\timesch: ' p_str]}, 'FontSize',24)
    elseif nch == 3
        title({[str_filt ' - ' ch_list{1} ', ' ch_list{2} ' and ' ch_list{3}],['state\timesch: ' p_str]}, 'FontSize',24)
    end

    ind_sp = ind_sp + 1;
end




end