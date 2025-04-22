function [ f, rm_cell, table_pairwise] = func_plot_EO_EC_alpha(data_eo, data_ec, filt_list, ch_list, atype)

% NOTE: using 'lao' for log-transformed oscillatory alpha power 

nf = length(filt_list);
nch = length(ch_list);
ind_sp = 1;

if ~contains(atype,{'ar','ao','lar','lao'})
    error('Wrong alpha power type (use ar, ao, lar or lao)')
end

if strcmp(atype,'ar')
    a_str = '\alpha_r_a_w';
elseif strcmp(atype,'ao')
    a_str = '\alpha_o_s_c_i';
elseif strcmp(atype,'lar')
    a_str = 'log(\alpha_r_a_w)';
elseif strcmp(atype,'lao')
    a_str = 'log(\alpha_o_s_c_i)';
end

rm_cell = cell(nch*nf,1);


if length(filt_list) == 2
    f = figure('color','w','units','normalized','outerposition',[0.15 0 0.7 0.95]);
elseif length(filt_list) == 3
    f = figure('color','w','units','normalized','outerposition',[0 0 1 0.95]);
end

% for pairwise comparisons
N = nch*nf;
var_cell = cell(N,4);

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
    zeros(N,1),...
    zeros(N,1),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    'VariableNames',{'var1','var2','E_eo','E_ec','p','h','pc','hc','st','df','ES','v1','v2','ttype','stat','ES_cell'});

for ch = 1:nch
    for ftype = 1:nf
        tmp_eo = data_eo.([atype '_' ch_list{ch}]).(filt_list{ftype});
        tmp_ec = data_ec.([atype '_' ch_list{ch}]).(filt_list{ftype});

        var_cell{(ch-1)*nch+ftype,1} = filt_list{ftype};
        var_cell{(ch-1)*nch+ftype,2} = ch_list{ch};
        var_cell{(ch-1)*nch+ftype,3} = tmp_eo;
        var_cell{(ch-1)*nch+ftype,4} = tmp_ec;

        Ns = length(tmp_eo);

        % ctype = data_eo.(['lo_' ch_list{ch}]).ctype{1};

        subplot(nch,nf,ind_sp)
        % create boxplots
        plot(2*ones(Ns,1),tmp_eo,'r*','MarkerSize',12)
        hold on
        plot(6*ones(Ns,1),tmp_ec,'b*','MarkerSize',12)

        % draw boxes
        pm_eo = fill([1 3 3 1],[quantile(tmp_eo,0.25), quantile(tmp_eo,0.25), quantile(tmp_eo,0.75) quantile(tmp_eo,0.75)],[1 0.5 0.5],'LineWidth',1.5,'FaceAlpha',0.4);
        pm_ec = fill([5 7 7 5],[quantile(tmp_ec,0.25), quantile(tmp_ec,0.25), quantile(tmp_ec,0.75) quantile(tmp_ec,0.75)],[0.5 0.5 1],'LineWidth',1.5,'FaceAlpha',0.4);

        % draw percentiles
        plot([2 2],[quantile(tmp_eo,0.75) quantile(tmp_eo,0.95)],'k-','LineWidth',1.5)
        plot([1.5 2.5],[quantile(tmp_eo,0.95) quantile(tmp_eo,0.95)],'k-','LineWidth',1.5)
        plot([2 2],[quantile(tmp_eo,0.25) quantile(tmp_eo,0.05)],'k-','LineWidth',1.5)
        plot([1.5 2.5],[quantile(tmp_eo,0.05) quantile(tmp_eo,0.05)],'k-','LineWidth',1.5)

        plot([6 6],[quantile(tmp_ec,0.75) quantile(tmp_ec,0.95)],'k-','LineWidth',1.5)
        plot([5.5 6.5],[quantile(tmp_ec,0.95) quantile(tmp_ec,0.95)],'k-','LineWidth',1.5)
        plot([6 6],[quantile(tmp_ec,0.25) quantile(tmp_ec,0.05)],'k-','LineWidth',1.5)
        plot([5.5 6.5],[quantile(tmp_ec,0.05) quantile(tmp_ec,0.05)],'k-','LineWidth',1.5)

        % plot means and medians
        plot([1 3],[mean(tmp_eo) mean(tmp_eo)],'r-', 'LineWidth',3)
        plot([5 7],[mean(tmp_ec) mean(tmp_ec)],'b-', 'LineWidth',3)

        plot([1 3],[median(tmp_eo) median(tmp_eo)],'r--', 'LineWidth',3)
        plot([5 7],[median(tmp_ec) median(tmp_ec)],'b--', 'LineWidth',3)


        ymin = min(get(gca,'ylim'));
        ymax = max(get(gca,'ylim'));
        yrange = ymax - ymin;

        % statistical tests - full rmANOVA model
        data_table = table();
        data_table.eo = tmp_eo;
        data_table.ec = tmp_ec;

        WithinDesign = table({'EO','EC'}','VariableNames',{'state'});
        WithinDesign.state = categorical(WithinDesign.state);

        rm_model = fitrm(data_table,'eo-ec ~ 1','WithinDesign',WithinDesign);
        rm = ranova(rm_model);

        rm_cell{ind_sp} = rm;


        % statistical tests - EO vs EC interaction
        if lillietest(tmp_eo,'Alpha',0.05/N) || lillietest(tmp_ec,'Alpha',0.05/N)
            [p,~] = signrank(tmp_eo,tmp_ec);
        else
            [~,p] = ttest(tmp_eo,tmp_ec);
        end

        p = p*N;

        if p<0.05
            plot([1 7],[ymax+0.1*yrange ymax+0.1*yrange],'k-','LineWidth',4)
            if p>=0.0001
                text(4,ymax+0.18*yrange,['\itp\rm=' num2str(p,'%.4f')],'FontSize',24,'HorizontalAlignment','center')
            else
                text(4,ymax+0.18*yrange,'\itp\rm<0.0001','FontSize',24,'HorizontalAlignment','center')
            end
        end

        % plot trends
        if lillietest(tmp_eo,'Alpha',0.0125) || lillietest(tmp_ec,'Alpha',0.0125)
            if p<0.05
                plot([2 6],[median(tmp_eo) median(tmp_ec)],'o-', 'color', [0.4 0.4 0.4], 'LineWidth', 3)
            else
                plot([2 6],[median(tmp_eo) median(tmp_ec)],'o:', 'color', [0.4 0.4 0.4], 'LineWidth', 3)
            end
        else
            if p<0.05
                plot([2 6],[mean(tmp_eo) mean(tmp_ec)],'o-', 'color', [0.4 0.4 0.4], 'LineWidth', 3)
            else
                plot([2 6],[mean(tmp_eo) mean(tmp_ec)],'o:', 'color', [0.4 0.4 0.4], 'LineWidth', 3)
            end
        end


        % plot specifics
        xlim([-1 9])
        ylim([ymin-0.3*yrange ymax+0.4*yrange])
        set(gca,'FontSize',20,'LineWidth',2,'XTick',[2 6],'XTickLabel',{'Eyes Open','Eyes Closed'})
        legend([pm_eo,pm_ec],{'EO','EC'},'location','southeast','FontSize',24)
        xlabel('Eyes Open vs. Eyes Closed','FontSize',24)
        ylabel(a_str,'FontSize',24)

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

% pairwise comparisons
for n = 1:N
    vname1 = [var_cell{n,1} '_' var_cell{n,2} '_EO'];
    vname2 = [var_cell{n,1} '_' var_cell{n,2} '_EC'];
    var1 = var_cell{n,3};
    var2 = var_cell{n,4};

    if lillietest(var1,'Alpha',0.05/N) || lillietest(var2,'Alpha',0.05/N)
        [p,h,s] = signrank(var1, var2);
        E_eo = median(var1);
        E_ec = median(var2);
        ttype = 'signrank';
        statvalue = s.zval;
        df = length(var1)-1;
        ES = abs(statvalue)/sqrt(length(var1)+length(var2));
        ES_tmp = ES;
    else
        [h,p,~,s] = ttest(var1, var2);
        E_eo = mean(var1);
        E_ec = mean(var2);
        ttype = 'ttest';
        statvalue = s.tstat;
        df = s.df;
        ES_tmp = meanEffectSize(var1,var2,'Effect','cohen');
        ES = abs(ES_tmp.Effect('CohensD'));
    end

    table_pairwise.var1{n} = vname1;
    table_pairwise.var2{n} = vname2;
    table_pairwise.E_eo(n) = E_eo;
    table_pairwise.E_ec(n) = E_ec;
    table_pairwise.p(n) = p;
    table_pairwise.h(n) = h;
    table_pairwise.pc(n) = p*N;
    table_pairwise.hc(n) = p*N < 0.05;
    table_pairwise.st(n) = statvalue;
    table_pairwise.df(n) = df;
    table_pairwise.ES(n) = ES;
    table_pairwise.v1{n} = var1;
    table_pairwise.v2{n} = var2;
    table_pairwise.ttype{n} = ttype;
    table_pairwise.stat{n} = s;
    table_pairwise.ES_cell{n} = ES_tmp;
end




end