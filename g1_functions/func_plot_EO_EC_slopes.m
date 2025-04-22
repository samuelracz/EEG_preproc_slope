function [ f, rm_table] = func_plot_EO_EC_slopes(data_eo, data_ec, filt_list, ch_list)

nf = length(filt_list);
nch = length(ch_list);
ind_sp = 1;

rm_table = table(cell(nch*nf,1), cell(nch*nf,1), cell(nch*nf,1), 'VariableNames', {'ch','ftype','rm'});
rm_ind = 1;

n_multcomp = 4; % EO vs. EC and Beta_lo vs. Beta_hi

% set figure size
if length(filt_list) == 1
    if length(ch_list) == 1
        f_x = 0.35;
        f_w = 0.3;
        f_y = 0.3;
        f_h = 0.55;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 2
        f_x = 0.35;
        f_w = 0.3;
        f_y = 0;
        f_h = 0.95;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 3
        f_x = 0;
        f_w = 1;
        f_y = 0.3;
        f_h = 0.55;
        sp1 = nf;
        sp2 = nch;
    end
elseif length(filt_list) == 2
    if length(ch_list) == 1
        f_x = 0.35;
        f_w = 0.3;
        f_y = 0.3;
        f_h = 0.55;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 2
        f_x = 0.15;
        f_w = 0.7;
        f_y = 0;
        f_h = 0.95;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 3
        f_x = 0;
        f_w = 1;
        f_y = 0;
        f_h = 0.95;
        sp1 = nf;
        sp2 = nch;
    end
elseif length(filt_list) == 3
    if length(ch_list) == 1
        f_x = 0;
        f_w = 1;
        f_y = 0.3;
        f_h = 0.55;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 2
        f_x = 0;
        f_w = 1;
        f_y = 0;
        f_h = 0.95;
        sp1 = nch;
        sp2 = nf;
    elseif length(ch_list) == 3
        f_x = 0;
        f_w = 1;
        f_y = 0;
        f_h = 0.95;
        sp1 = nch;
        sp2 = nf;
    end
end

f = figure('color','w','units','normalized','outerposition',[f_x f_y f_w f_h]);

for ch = 1:nch
    for ftype = 1:nf
        tmp_lo_eo = data_eo.(['lo_' ch_list{ch}]).(filt_list{ftype});
        tmp_lo_ec = data_ec.(['lo_' ch_list{ch}]).(filt_list{ftype});
        tmp_hi_eo = data_eo.(['hi_' ch_list{ch}]).(filt_list{ftype});
        tmp_hi_ec = data_ec.(['hi_' ch_list{ch}]).(filt_list{ftype});

        N = length(tmp_lo_eo);

        subplot(sp1,sp2,ind_sp)
        % create boxplots
        plot(2*ones(N,1),tmp_lo_eo,'r*','MarkerSize',12)
        hold on
        plot(5*ones(N,1),tmp_lo_ec,'b*','MarkerSize',12)
        plot(12*ones(N,1),tmp_hi_eo,'rd','MarkerSize',12)
        plot(15*ones(N,1),tmp_hi_ec,'bd','MarkerSize',12)

        % draw boxes
        pm_eo = fill([1 3 3 1],[quantile(tmp_lo_eo,0.25), quantile(tmp_lo_eo,0.25), quantile(tmp_lo_eo,0.75) quantile(tmp_lo_eo,0.75)],[1 0.5 0.5],'LineWidth',1.5,'FaceAlpha',0.4);
        pm_ec = fill([4 6 6 4],[quantile(tmp_lo_ec,0.25), quantile(tmp_lo_ec,0.25), quantile(tmp_lo_ec,0.75) quantile(tmp_lo_ec,0.75)],[0.5 0.5 1],'LineWidth',1.5,'FaceAlpha',0.4);
        fill([11 13 13 11],[quantile(tmp_hi_eo,0.25), quantile(tmp_hi_eo,0.25), quantile(tmp_hi_eo,0.75) quantile(tmp_hi_eo,0.75)],[1 0.5 0.5],'LineWidth',1.5,'FaceAlpha',0.4)
        fill([14 16 16 14],[quantile(tmp_hi_ec,0.25), quantile(tmp_hi_ec,0.25), quantile(tmp_hi_ec,0.75) quantile(tmp_hi_ec,0.75)],[0.5 0.5 1],'LineWidth',1.5,'FaceAlpha',0.4)

        % draw percentiles
        plot([2 2],[quantile(tmp_lo_eo,0.75) quantile(tmp_lo_eo,0.95)],'k-','LineWidth',1.5)
        plot([1.5 2.5],[quantile(tmp_lo_eo,0.95) quantile(tmp_lo_eo,0.95)],'k-','LineWidth',1.5)
        plot([2 2],[quantile(tmp_lo_eo,0.25) quantile(tmp_lo_eo,0.05)],'k-','LineWidth',1.5)
        plot([1.5 2.5],[quantile(tmp_lo_eo,0.05) quantile(tmp_lo_eo,0.05)],'k-','LineWidth',1.5)

        plot([5 5],[quantile(tmp_lo_ec,0.75) quantile(tmp_lo_ec,0.95)],'k-','LineWidth',1.5)
        plot([4.5 5.5],[quantile(tmp_lo_ec,0.95) quantile(tmp_lo_ec,0.95)],'k-','LineWidth',1.5)
        plot([5 5],[quantile(tmp_lo_ec,0.25) quantile(tmp_lo_ec,0.05)],'k-','LineWidth',1.5)
        plot([4.5 5.5],[quantile(tmp_lo_ec,0.05) quantile(tmp_lo_ec,0.05)],'k-','LineWidth',1.5)

        plot([12 12],[quantile(tmp_hi_eo,0.75) quantile(tmp_hi_eo,0.95)],'k-','LineWidth',1.5)
        plot([11.5 12.5],[quantile(tmp_hi_eo,0.95) quantile(tmp_hi_eo,0.95)],'k-','LineWidth',1.5)
        plot([12 12],[quantile(tmp_hi_eo,0.25) quantile(tmp_hi_eo,0.05)],'k-','LineWidth',1.5)
        plot([11.5 12.5],[quantile(tmp_hi_eo,0.05) quantile(tmp_hi_eo,0.05)],'k-','LineWidth',1.5)

        plot([15 15],[quantile(tmp_hi_ec,0.75) quantile(tmp_hi_ec,0.95)],'k-','LineWidth',1.5)
        plot([14.5 15.5],[quantile(tmp_hi_ec,0.95) quantile(tmp_hi_ec,0.95)],'k-','LineWidth',1.5)
        plot([15 15],[quantile(tmp_hi_ec,0.25) quantile(tmp_hi_ec,0.05)],'k-','LineWidth',1.5)
        plot([14.5 15.5],[quantile(tmp_hi_ec,0.05) quantile(tmp_hi_ec,0.05)],'k-','LineWidth',1.5)

        % plot means and medians
        plot([1 3],[mean(tmp_lo_eo) mean(tmp_lo_eo)],'r-', 'LineWidth',3)
        plot([4 6],[mean(tmp_lo_ec) mean(tmp_lo_ec)],'b-', 'LineWidth',3)
        plot([11 13],[mean(tmp_hi_eo) mean(tmp_hi_eo)],'r-', 'LineWidth',3)
        plot([14 16],[mean(tmp_hi_ec) mean(tmp_hi_ec)],'b-', 'LineWidth',3)

        plot([1 3],[median(tmp_lo_eo) median(tmp_lo_eo)],'r--', 'LineWidth',3)
        plot([4 6],[median(tmp_lo_ec) median(tmp_lo_ec)],'b--', 'LineWidth',3)
        plot([11 13],[median(tmp_hi_eo) median(tmp_hi_eo)],'r--', 'LineWidth',3)
        plot([14 16],[median(tmp_hi_ec) median(tmp_hi_ec)],'b--', 'LineWidth',3)

        % plot trends
        if lillietest(tmp_lo_eo,'Alpha',0.05/n_multcomp) || lillietest(tmp_hi_eo,'Alpha',0.05/n_multcomp)
            [p,~] = signrank(tmp_lo_eo,tmp_hi_eo);
            p = p*n_multcomp;
            if p<0.05
                plot([2 12],[median(tmp_lo_eo) median(tmp_hi_eo)],'ro-', 'LineWidth', 3)
            else
                plot([2 12],[median(tmp_lo_eo) median(tmp_hi_eo)],'ro:', 'LineWidth', 3)
            end
        else
            [~,p] = ttest(tmp_lo_eo,tmp_hi_eo);
            p = p*n_multcomp;
            if p<0.05
                plot([2 12],[mean(tmp_lo_eo) mean(tmp_hi_eo)],'ro-', 'LineWidth', 3)
            else
                plot([2 12],[mean(tmp_lo_eo) mean(tmp_hi_eo)],'ro:', 'LineWidth', 3)
            end
        end

        if lillietest(tmp_lo_ec,'Alpha',0.05/n_multcomp) || lillietest(tmp_hi_ec,'Alpha',0.05/n_multcomp)
            [p,~] = signrank(tmp_lo_ec,tmp_hi_ec);
            p = p*n_multcomp;
            if p<0.05
                plot([5 15],[median(tmp_lo_ec) median(tmp_hi_ec)],'bo-', 'LineWidth', 3)
            else
                plot([5 15],[median(tmp_lo_ec) median(tmp_hi_ec)],'bo:', 'LineWidth', 3)
            end
        else
            [~,p] = ttest(tmp_lo_ec, tmp_hi_ec);
            p = p*n_multcomp;
            if p<0.05
                plot([5 15],[mean(tmp_lo_ec) mean(tmp_hi_ec)],'bo-', 'LineWidth', 3)
            else
                plot([5 15],[mean(tmp_lo_ec) mean(tmp_hi_ec)],'bo:', 'LineWidth', 3)
            end
        end

        ymin = min(get(gca,'ylim'));
        ymax = max(get(gca,'ylim'));
        yrange = ymax - ymin;

        ylim([ymin-0.3*yrange ymax+0.4*yrange])

        % statistical tests - full rmANOVA model
        data_table = table();
        data_table.lo_eo = tmp_lo_eo;
        data_table.hi_eo = tmp_hi_eo;
        data_table.lo_ec = tmp_lo_ec;
        data_table.hi_ec = tmp_hi_ec;

        WithinDesign = table({'EO','EO','EC','EC'}',{'low','high','low','high'}','VariableNames',{'state','range'});
        WithinDesign.state = categorical(WithinDesign.state);
        WithinDesign.range = categorical(WithinDesign.range);

        rm_model = fitrm(data_table,'lo_eo-hi_ec ~ 1','WithinDesign',WithinDesign);
        rm = ranova(rm_model,'WithinModel','state*range');

        % rm_cell{ind_sp} = rm;
        rm_table.ch{rm_ind} = ch_list{ch};
        rm_table.ftype{rm_ind} = filt_list{ftype};
        rm_table.rm{rm_ind} = rm;
        rm_ind = rm_ind + 1;

        % statistical tests - low- vs. high-range
        if lillietest(tmp_lo_eo,'Alpha',0.05/n_multcomp) || lillietest(tmp_hi_eo,'Alpha',0.05/n_multcomp)
            [p,~] = signrank(tmp_lo_eo,tmp_hi_eo);
        else
            [~,p] = ttest(tmp_lo_eo,tmp_hi_eo);
        end

        p = p*n_multcomp;

        if p<0.05
            plot([1 13],[ymin ymin],'k-','LineWidth',4)
            if p>=0.0001
                text(7,ymin-0.1*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text(7,ymin-0.1*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end

        if lillietest(tmp_lo_ec,'Alpha',0.0125) || lillietest(tmp_hi_ec,'Alpha',0.0125)
            [p,~] = signrank(tmp_lo_ec,tmp_hi_ec);
        else
            [~,p] = ttest(tmp_lo_ec,tmp_hi_ec);
        end

        p = p*n_multcomp;

        if p<0.05
            plot([4 16],[ymax+0.2*yrange ymax+0.2*yrange],'k-','LineWidth',4)
            if p>=0.0001
                text(10,ymax+0.28*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text(10,ymax+0.28*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end

        % statistical tests - EO vs EC
        if lillietest(tmp_lo_eo,'Alpha',0.0125) || lillietest(tmp_lo_ec,'Alpha',0.0125)
            [p,~] = signrank(tmp_lo_eo,tmp_lo_ec);
        else
            [~,p] = ttest(tmp_lo_eo,tmp_lo_ec);
        end

        p = p*n_multcomp;

        if p<0.05
            plot([1 6],[ymax ymax],'k-','LineWidth',4)
            if p>=0.0001
                text(3.5,ymax+0.08*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text(3.5,ymax-0.08*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end

        if lillietest(tmp_hi_eo,'Alpha',0.0125) || lillietest(tmp_hi_ec,'Alpha',0.0125)
            [p,~] = signrank(tmp_hi_eo,tmp_hi_ec);
        else
            [~,p] = ttest(tmp_hi_eo,tmp_hi_ec);
        end

        p = p*n_multcomp;

        if p<0.05
            plot([11 16],[ymax ymax],'k-','LineWidth',4)
            if p>=0.0001
                text(13.5,ymax+0.08*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text(13.5,ymax+0.08*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end

        % plot specifics
        xlim([0 17])
        set(gca,'FontSize',20,'LineWidth',2,'XTick',[2 5 12 15],'XTickLabel',{'EO','EC','EO','EC'})
        legend([pm_eo,pm_ec],{'EO','EC'},'location','southeast','FontSize',24)
        xlabel('\beta_{low} (left) vs. \beta_{high} (right)','FontSize',24)
        ylabel('Spectral Slope','FontSize',24)
        p = rm.pValueGG('(Intercept):state:range');
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
        
        title({[str_filt ' - ' ch_list{ch}],['state\timesrange: ' p_str]}, 'FontSize',24)

        ind_sp = ind_sp + 1;
    end
end




end