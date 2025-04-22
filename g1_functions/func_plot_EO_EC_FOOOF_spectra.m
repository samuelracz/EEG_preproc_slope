function [ f ] = func_plot_EO_EC_FOOOF_spectra(data_eo, data_ec, filt_list, ch_list)

nf = length(filt_list);
nch = length(ch_list);
ind_sp = 1;

if length(filt_list) == 2
    f = figure('color','w','units','normalized','outerposition',[0.15 0 0.7 0.95]);
elseif length(filt_list) == 3
    f = figure('color','w','units','normalized','outerposition',[0 0 1 0.95]);
end

for ch = 1:nch
    for ftype = 1:nf
        tmp_eo = data_eo.(ch_list{ch}).(filt_list{ftype});
        tmp_ec = data_ec.(ch_list{ch}).(filt_list{ftype});

        % ctype = tmp_eo.ctype;

        freq_eo = tmp_eo.freq(:);
        freq_ec = tmp_ec.freq(:);
        freq_lo_eo = tmp_eo.Freq_lo;
        freq_hi_eo = tmp_eo.Freq_hi;
        freq_lo_ec = tmp_ec.Freq_lo;
        freq_hi_ec = tmp_ec.Freq_hi;

        mu_mixd_eo = mean(tmp_eo.mixd,2);
        mu_frac_eo = mean(tmp_eo.frac,2);
        std_pos_eo = mean(tmp_eo.mixd,2) + std(tmp_eo.mixd,[],2);
        std_neg_eo = mean((tmp_eo.mixd),2) - std(tmp_eo.mixd,[],2);
        mu_plaw_lo_eo = mean(tmp_eo.Plaw_lo,2);
        mu_plaw_hi_eo = mean(tmp_eo.Plaw_hi,2);
        mu_mixd_ec = mean(tmp_ec.mixd,2);
        mu_frac_ec = mean(tmp_ec.frac,2);
        std_pos_ec = mean(tmp_ec.mixd,2) + std(tmp_ec.mixd,[],2);
        std_neg_ec = mean(tmp_ec.mixd,2) - std(tmp_ec.mixd,[],2);
        mu_plaw_lo_ec = mean(tmp_ec.Plaw_lo,2);
        mu_plaw_hi_ec = mean(tmp_ec.Plaw_hi,2);

        subplot(nch,nf,ind_sp)
        plot(freq_eo, mu_mixd_eo, 'r')
        hold on
        plot(freq_ec, mu_mixd_ec, 'b')

        % error ranges
        fill([freq_eo; flipud(freq_eo)],[mu_mixd_eo; flipud(std_pos_eo)],[1 0.5 0.5],'FaceAlpha',0.2,'EdgeAlpha',0)
        fill([freq_eo; flipud(freq_eo)],[mu_mixd_eo; flipud(std_neg_eo)],[1 0.5 0.5],'FaceAlpha',0.2,'EdgeAlpha',0)
        fill([freq_ec; flipud(freq_ec)],[mu_mixd_ec; flipud(std_pos_ec)],[0.5 0.5 1],'FaceAlpha',0.2,'EdgeAlpha',0)
        fill([freq_ec; flipud(freq_ec)],[mu_mixd_ec; flipud(std_neg_ec)],[0.5 0.5 1],'FaceAlpha',0.2,'EdgeAlpha',0)

        % means
        pm_eo = plot(freq_eo, mu_mixd_eo, ':', 'color', [1 0.5 0.5], 'LineWidth' ,3);
        plot(freq_eo, mu_frac_eo, '-', 'color', [0.7 0 0], 'LineWidth' ,3);
        plot(freq_lo_eo, mu_plaw_lo_eo, 'r--', 'LineWidth',3);
        plot(freq_hi_eo, mu_plaw_hi_eo, 'r--', 'LineWidth',3);
        pm_ec = plot(freq_ec, mu_mixd_ec, ':', 'color', [0.5 0.5 1], 'LineWidth' ,3);
        plot(freq_ec, mu_frac_ec, '-', 'color', [0 0 0.7], 'LineWidth' ,3);
        plot(freq_lo_ec, mu_plaw_lo_ec, 'b--', 'LineWidth',3);
        plot(freq_hi_ec, mu_plaw_hi_ec, 'b--', 'LineWidth',3);

        % plot scaling range limits
        xmin_l = min(freq_lo_eo);
        xmax_l = max(freq_lo_eo);
        xmin_h = min(freq_hi_eo);
        xmax_h = max(freq_hi_eo);

        xlim([xmin_l xmax_h])
        ymin_tmp = min(get(gca,'ylim'));
        ymax_tmp = max(get(gca,'ylim'));

        ymin = ymin_tmp;
        ymax = ymax_tmp;
        ylim([ymin ymax])

        % plot exclusion range
        plot([xmax_l xmax_l],[ymin ymax],'k--')
        plot([xmin_h xmin_h],[ymin ymax],'k--')
        fe = fill([xmax_l xmin_h xmin_h xmax_l],[ymin ymin ymax ymax],[0.4 0.4 0.4],'FaceAlpha',0.2);

        % plot specifics
        set(gca,'FontSize',20,'LineWidth',2,'xscale','log')
        legend([pm_eo,pm_ec,fe],{'Eyes Open','Eyes Closed','Exclusion Range'},'location','southwest','FontSize',24)
        xlabel('log frequency','FontSize',24)
        ylabel('log power','FontSize',24)
        if strcmp(filt_list{ftype},'raw')
            str_filt = 'Raw EEG';
        elseif strcmp(filt_list{ftype},'car')
            str_filt = 'Common Average Reference';
        elseif strcmp(filt_list{ftype},'lap')
            str_filt = 'Small Laplacian Filtering';
        end
        title([str_filt ' - ' ch_list{ch}], 'FontSize',24)

        ind_sp = ind_sp + 1;
    end
end




end