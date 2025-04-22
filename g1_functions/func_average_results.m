function [ data_out ] = func_average_results(input_struct, chlab, f_list)

% preallocate
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

% iterate over filter types (raw, car, lap)
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
            mat_mixd = zeros([size(tmp.mixd),N]);
            mat_frac = zeros([size(tmp.frac),N]);
            mat_osci = zeros([size(tmp.osci),N]);
        end

        if isempty(mat_plaw_bb) && isempty(mat_plaw_lo) && isempty(mat_plaw_hi)
            mat_plaw_bb = zeros([size(tmp.Plaw_bb),N]);
            mat_plaw_lo = zeros([size(tmp.Plaw_lo),N]);
            mat_plaw_hi = zeros([size(tmp.Plaw_hi),N]);
        end

        mat_mixd(:,:,n) = tmp.mixd;
        mat_frac(:,:,n) = tmp.frac;
        mat_osci(:,:,n) = tmp.osci;
        mat_plaw_bb(:,:,n) = tmp.Plaw_bb;
        mat_plaw_lo(:,:,n) = tmp.Plaw_lo;
        mat_plaw_hi(:,:,n) = tmp.Plaw_hi;

        freq = tmp.freq;
        mixd = tmp.mixd;
        osci = tmp.osci;
        log_mixd = log(tmp.mixd);
        log_osci = log(tmp.mixd) - log(tmp.frac);

        mat_beta_bb(n,:) = tmp.Beta_bb;
        mat_cons_bb(n,:) = tmp.Cons_bb;
        mat_beta_lo(n,:) = tmp.Beta_lo;
        mat_cons_lo(n,:) = tmp.Cons_lo;
        mat_beta_hi(n,:) = tmp.Beta_hi;
        mat_cons_hi(n,:) = tmp.Cons_hi;

        mat_alpha_raw(n,:) = sum(mixd(freq>=8 & freq<=12,:),1);
        mat_alpha_osci(n,:) = sum(osci(freq>=8 & freq<=12,:),1);
        mat_log_alpha_raw(n,:) = sum(log_mixd(freq>=8 & freq<=12,:),1);
        mat_log_alpha_osci(n,:) = sum(log_osci(freq>=8 & freq<=12,:),1);
    end

    data_out.(f_list{f}).freq = tmp.freq;
    data_out.(f_list{f}).srate = tmp.srate;
    data_out.(f_list{f}).mixd = mean(mat_mixd,3);
    data_out.(f_list{f}).frac = mean(mat_frac,3);
    data_out.(f_list{f}).osci = mean(mat_osci,3);
    data_out.(f_list{f}).Beta_bb = mean(mat_beta_bb,1);
    data_out.(f_list{f}).Cons_bb = mean(mat_cons_bb,1);
    data_out.(f_list{f}).Plaw_bb = mean(mat_plaw_bb,3);
    data_out.(f_list{f}).Freq_bb = tmp.Freq_bb;
    data_out.(f_list{f}).Beta_lo = mean(mat_beta_lo,1);
    data_out.(f_list{f}).Cons_lo = mean(mat_cons_lo,1);
    data_out.(f_list{f}).Plaw_lo = mean(mat_plaw_lo,3);
    data_out.(f_list{f}).Freq_lo = tmp.Freq_lo;
    data_out.(f_list{f}).Beta_hi = mean(mat_beta_hi,1);
    data_out.(f_list{f}).Cons_hi = mean(mat_cons_hi,1);
    data_out.(f_list{f}).Plaw_hi = mean(mat_plaw_hi,3);
    data_out.(f_list{f}).Freq_hi = tmp.Freq_hi;
    data_out.(f_list{f}).Alpha_raw = mean(mat_alpha_raw,1);
    data_out.(f_list{f}).Alpha_osci = mean(mat_alpha_osci,1);
    data_out.(f_list{f}).logAlpha_raw = mean(mat_log_alpha_raw,1);
    data_out.(f_list{f}).logAlpha_osci = mean(mat_log_alpha_osci,1);
end

end