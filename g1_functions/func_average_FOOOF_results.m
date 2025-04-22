function [ data_out ] = func_average_FOOOF_results(input_struct, chlab, f_list)

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
            mat_mixd = zeros([length(tmp.FOOOF_bb{1,1}.freqs),length(tmp.FOOOF_bb{1,1}),N]);
            mat_frac = zeros([length(tmp.FOOOF_bb{1,1}.freqs),length(tmp.FOOOF_bb{1,1}),N]);
            mat_osci = zeros([length(tmp.FOOOF_bb{1,1}.freqs),length(tmp.FOOOF_bb{1,1}),N]);
        end

        if isempty(mat_plaw_bb) && isempty(mat_plaw_lo) && isempty(mat_plaw_hi)
            mat_plaw_bb = zeros([length(tmp.FOOOF_bb{1,1}.freqs),length(tmp.FOOOF_bb{1,1}),N]);
            mat_plaw_lo = zeros([length(tmp.FOOOF_lo{1,1}.freqs),length(tmp.FOOOF_lo{1,1}),N]);
            mat_plaw_hi = zeros([length(tmp.FOOOF_hi{1,1}.freqs),length(tmp.FOOOF_hi{1,1}),N]);
        end
        
        for ch = 1:nch
            mat_mixd(:,ch,n) = tmp.FOOOF_bb{ch,1}.power_spectrum;
            mat_frac(:,ch,n) = tmp.FOOOF_bb{ch,1}.ap_fit;
            mat_osci(:,ch,n) = tmp.FOOOF_bb{ch,1}.fooofed_spectrum - tmp.FOOOF_bb{ch,1}.ap_fit;
            mat_plaw_bb(:,ch,n) = tmp.FOOOF_bb{ch,1}.ap_fit;
            mat_plaw_lo(:,ch,n) = tmp.FOOOF_lo{ch,1}.ap_fit;
            mat_plaw_hi(:,ch,n) = tmp.FOOOF_hi{ch,1}.ap_fit;
            freq_bb = tmp.FOOOF_bb{ch,1}.freqs;
            freq_lo = tmp.FOOOF_lo{ch,1}.freqs;
            freq_hi = tmp.FOOOF_hi{ch,1}.freqs;
        end

        freq = tmp.FOOOF_bb{1,1}.freqs;
        mixd = tmp.mixd;
        osci = tmp.osci;
        log_mixd = log(tmp.mixd);
        log_osci = log(tmp.mixd) - log(tmp.frac);

        mat_beta_bb(n,:) = tmp.Beta_bb_FOOOF;
        mat_cons_bb(n,:) = tmp.Cons_bb_FOOOF;
        mat_beta_lo(n,:) = tmp.Beta_lo_FOOOF;
        mat_cons_lo(n,:) = tmp.Cons_lo_FOOOF;
        mat_beta_hi(n,:) = tmp.Beta_hi_FOOOF;
        mat_cons_hi(n,:) = tmp.Cons_hi_FOOOF;

        mat_alpha_raw(n,:) = exp(tmp.Alpha_bb_FOOOF);
        mat_alpha_osci(n,:) = exp(tmp.Aosci_bb_FOOOF);
        mat_log_alpha_raw(n,:) = tmp.Alpha_bb_FOOOF;
        mat_log_alpha_osci(n,:) = tmp.Aosci_bb_FOOOF;
    end


    data_out.(f_list{f}).freq = freq;
    data_out.(f_list{f}).srate = tmp.srate;
    data_out.(f_list{f}).mixd = mean(mat_mixd,3);
    data_out.(f_list{f}).frac = mean(mat_frac,3);
    data_out.(f_list{f}).osci = mean(mat_osci,3);
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