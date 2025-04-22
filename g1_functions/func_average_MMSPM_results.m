function [ data_out ] = func_average_MMSPM_results(input_struct, chlab, f_list)

struct_tmp = struct(...
    'freq', [],...
    'srate', [],...
    'mixd', [],...
    'frac', [],...
    'osci', [],...
    'Beta_bb', [],...
    'Cons_bb', [],...
    'Plaw_bb', [],...
    'Freq_bb', [],...
    'Beta_lo', [],...
    'Cons_lo', [],...
    'Plaw_lo', [],...
    'Freq_lo', [],...
    'Beta_hi', [],...
    'Cons_hi', [],...
    'Plaw_hi', [],...
    'Freq_hi', [],...
    'Alpha_raw', [],...
    'Alpha_osci', [],...
    'logAlpha_raw', [],...
    'logAlpha_osci', []);

data_out = struct();

for f = 1:length(f_list)
    data_out.(f_list{f}) = struct_tmp;
end

data_in = input_struct;
N = length(data_in);
nch = length(chlab);

for f = 1:length(f_list)
    mat_mixd = [];
    mat_frac = [];
    mat_osci = [];
    mat_beta_bb = zeros(N,nch);
    mat_cons_bb = zeros(N,nch);
    mat_beta_lo = zeros(N,nch);
    mat_cons_lo = zeros(N,nch);
    mat_beta_hi = zeros(N,nch);
    mat_cons_hi = zeros(N,nch);
    mat_beta_bb_fooof = zeros(N,nch);
    mat_cons_bb_fooof = zeros(N,nch);
    mat_beta_lo_fooof = zeros(N,nch);
    mat_cons_lo_fooof = zeros(N,nch);
    mat_beta_hi_fooof = zeros(N,nch);
    mat_cons_hi_fooof = zeros(N,nch);
    mat_alpha_raw = zeros(N,nch);
    mat_alpha_osci = zeros(N,nch);
    mat_log_alpha_raw = zeros(N,nch);
    mat_log_alpha_osci = zeros(N,nch);
    mat_plaw_bb = [];
    mat_plaw_lo = [];
    mat_plaw_hi = [];
    for n = 1:N
        tmp = data_in(n).(f_list{f});
        if isempty(mat_mixd) && isempty(mat_frac) && isempty(mat_osci)
            ind_bb = tmp.MMSPM{1,1}.freq>=1 & tmp.MMSPM{1,1}.freq<=45;
            ind_lo = tmp.MMSPM{1,1}.freq>=1 & tmp.MMSPM{1,1}.freq<=4;
            ind_hi = tmp.MMSPM{1,1}.freq>=20 & tmp.MMSPM{1,1}.freq<=45;
            mat_mixd = zeros([sum(ind_bb),nch,N]);
            mat_frac = zeros([sum(ind_bb),nch,N]);
            mat_osci = zeros([sum(ind_bb),nch,N]);
        end

        if isempty(mat_plaw_bb) && isempty(mat_plaw_lo) && isempty(mat_plaw_hi)
            mat_plaw_bb = zeros([sum(ind_bb),nch,N]);
            mat_plaw_lo = zeros([sum(ind_lo),nch,N]);
            mat_plaw_hi = zeros([sum(ind_hi),nch,N]);
        end

        for ch = 1:nch
            mat_mixd(:,ch,n) = tmp.MMSPM{ch,1}.Y;
            mat_frac(:,ch,n) = tmp.MMSPM{ch,1}.A;
            mat_osci(:,ch,n) = tmp.MMSPM{ch,1}.Y_mmspm - tmp.MMSPM{ch,1}.A;
            mat_plaw_bb(:,ch,n) = tmp.MMSPM{ch,1}.A;
            mat_plaw_lo(:,ch,n) = tmp.MMSPM{ch,1}.A(ind_lo);
            mat_plaw_hi(:,ch,n) = tmp.MMSPM{ch,1}.A(ind_hi);
            freq_bb = tmp.MMSPM{ch,1}.freq;
            freq_lo = tmp.MMSPM{ch,1}.freq(ind_lo);
            freq_hi = tmp.MMSPM{ch,1}.freq(ind_hi);
        end

        mat_beta_bb(n,:) = tmp.Beta_bb_MMSPM;
        mat_cons_bb(n,:) = tmp.Cons_bb_MMSPM;
        mat_beta_lo(n,:) = tmp.Beta_lo_MMSPM;
        mat_cons_lo(n,:) = tmp.Cons_lo_MMSPM;
        mat_beta_hi(n,:) = tmp.Beta_hi_MMSPM;
        mat_cons_hi(n,:) = tmp.Cons_hi_MMSPM;

        mat_alpha_raw(n,:) = exp(tmp.Alpha_bb_MMSPM);
        mat_alpha_osci(n,:) = exp(tmp.Aosci_bb_MMSPM);
        mat_log_alpha_raw(n,:) = tmp.Alpha_bb_MMSPM;
        mat_log_alpha_osci(n,:) = tmp.Aosci_bb_MMSPM;
    end

    data_out.(f_list{f}).freq = tmp.MMSPM{1,1}.freq;
    data_out.(f_list{f}).srate = tmp.srate;
    data_out.(f_list{f}).mixd = mean(mat_mixd(ind_bb,:,:),3);
    data_out.(f_list{f}).frac = mean(mat_frac(ind_bb,:,:),3);
    data_out.(f_list{f}).osci = mean(mat_osci(ind_bb,:,:),3);
    data_out.(f_list{f}).Beta_bb = mean(mat_beta_bb,1);
    data_out.(f_list{f}).Cons_bb = mean(mat_cons_bb,1);
    data_out.(f_list{f}).Plaw_bb = mean(mat_plaw_bb,3);
    data_out.(f_list{f}).Freq_bb = freq_bb;
    data_out.(f_list{f}).Beta_lo = mean(mat_beta_lo,1);
    data_out.(f_list{f}).Cons_lo = mean(mat_cons_lo,1);
    data_out.(f_list{f}).Plaw_lo = mean(mat_plaw_lo,3);
    data_out.(f_list{f}).Freq_lo = freq_lo;
    data_out.(f_list{f}).Beta_hi = mean(mat_beta_hi,1);
    data_out.(f_list{f}).Cons_hi = mean(mat_cons_hi,1);
    data_out.(f_list{f}).Plaw_hi = mean(mat_plaw_hi,3);
    data_out.(f_list{f}).Freq_hi = freq_hi;
    data_out.(f_list{f}).Alpha_raw = mean(mat_alpha_raw,1);
    data_out.(f_list{f}).Alpha_osci = mean(mat_alpha_osci,1);
    data_out.(f_list{f}).logAlpha_raw = mean(mat_log_alpha_raw,1);
    data_out.(f_list{f}).logAlpha_osci = mean(mat_log_alpha_osci,1);
end

end