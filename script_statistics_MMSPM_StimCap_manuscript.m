OD = pwd;
fnames_stimcap = dir('results_MMSPM_StimCap/*.mat');
addpath(genpath('g1_functions'))

ns_stimcap = length(fnames_stimcap);
chlab = {'FC1','FC2','CP1','CP2','POz'};
filt_list_stat = {'raw','car','lap'};
filt_list = {'car','lap'};


%% load data

try
    load('ws_stat_MMSPM_StimCap_EC.mat')
catch

    % iterate over subjects
    subj_stimcap_eo = cell(ns_stimcap,1);
    subj_stimcap_ec = cell(ns_stimcap,1);

    for subj = 1:ns_stimcap
        load(['results_MMSPM_StimCap/' fnames_stimcap(subj).name]);
        data_eo = func_average_MMSPM_results(results_eo, chlab, filt_list_stat);
        data_ec = func_average_MMSPM_results(results_ec, chlab, filt_list_stat);

        subj_stimcap_eo{subj,1} = data_eo;
        subj_stimcap_ec{subj,1} = data_ec;
    end

    % synthesize data
    [data_spec_stimcap_eo] = func_structure_spectra(subj_stimcap_eo, chlab, 'EO', 'StimCap', filt_list_stat);
    [data_spec_stimcap_ec] = func_structure_spectra(subj_stimcap_ec, chlab, 'EC', 'StimCap', filt_list_stat);

    [data_mmspm_stimcap_eo] = func_structure_data(subj_stimcap_eo, chlab, 'EO', 'StimCap', filt_list_stat);
    [data_mmspm_stimcap_ec] = func_structure_data(subj_stimcap_ec, chlab, 'EC', 'StimCap', filt_list_stat);

    DT = datetime;
    save('ws_stat_MMSPM_StimCap_EC.mat','DT',...
        'data_spec_stimcap_eo', 'data_spec_stimcap_ec',...
        'data_mmspm_stimcap_eo', 'data_mmspm_stimcap_ec');
end


%% statistics and plotting
fsp2 = func_plot_EO_EC_FOOOF_spectra(data_spec_stimcap_eo,data_spec_stimcap_ec,filt_list,{'FC2','POz'});
exportgraphics(fsp2,'h1_figures/MMSPM/fig_EO_EC_spectra_StimCap_MMSPM.png','Resolution',300)
close(fsp2)

[fsl2, stat_slope_stimcap] = func_plot_EO_EC_slopes(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'FC2','POz'});
exportgraphics(fsl2,'h1_figures/MMSPM/fig_EO_EC_slopes_StimCap_MMSPM.png','Resolution',300)
close(fsl2)

[fslfc1, stat_slope_stimcap_fc1] = func_plot_EO_EC_slopes(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'FC1','POz'});
exportgraphics(fslfc1,'h1_figures/MMSPM/fig_EO_EC_slopes_StimCap_FC1ctrl_MMSPM.png','Resolution',300)
close(fslfc1)

[fslr, stat_slope_stimcap_raw] = func_plot_EO_EC_slopes(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,{'raw'},{'FC1','POz','FC2'});
exportgraphics(fslr,'h1_figures/MMSPM/fig_EO_EC_slopes_StimCap_raw_MMSPM.png','Resolution',300)
close(fslr)

[fslcp, stat_slope_stimcap_cp] = func_plot_EO_EC_slopes(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'CP1','CP2'});
exportgraphics(fslcp,'h1_figures/MMSPM/fig_EO_EC_slopes_StimCap_CP1CP2ctrl_MMSPM.png','Resolution',300)
close(fslcp)

[flo2, ~, pw_lao_fc2] = func_plot_EO_EC_alpha(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'FC2','POz'},'lao');
exportgraphics(flo2,'h1_figures/MMSPM/fig_EO_EC_lao_StimCap_MMSPM.png','Resolution',300)
close(flo2)

[flo3, ~, pw_lao_fc3] = func_plot_EO_EC_alpha(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'FC1','POz'},'lao');
exportgraphics(flo3,'h1_figures/MMSPM/fig_EO_EC_lao_StimCap_FC1ctrl_MMSPM.png','Resolution',300)
close(flo3)

% broadband
fbl2 = func_plot_EO_EC_bbslopes(data_mmspm_stimcap_eo,data_mmspm_stimcap_ec,filt_list,{'FC1','POz','FC2'});
exportgraphics(fbl2,'h1_figures/MMSPM/fig_EO_EC_broadband_slopes_StimCap_MMSPM.png','Resolution',300)
close(fbl2)

%% plotting correlations between alpha power and slope - StimCap
[fcorr_lao,corr_stat_lao_stimcap]=func_plot_corr_EO_EC_sa(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, filt_list, {'FC2','POz'},'lao', 4*length(filt_list));
exportgraphics(fcorr_lao,'h1_figures/MMSPM/fig_alpha_vs_beta_lao_stimcap_MMSPM.png','Resolution',300)
close(fcorr_lao)

[fcorr_lao_bb,corr_stat_lao_bb_stimcap]=func_plot_corr_EO_EC_bb(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, filt_list, {'FC2','POz'},'lao',2*length(filt_list));
exportgraphics(fcorr_lao_bb,'h1_figures/MMSPM/fig_alpha_vs_beta_bb_lao_stimcap_MMSPM.png','Resolution',300)
close(fcorr_lao_bb)

%% Statistics: full model range*condition*filter*location (4-way rmANOVA) for StimCap
% Within effect 1: spatial filter - raw (none) vs. CAR vs. Small Laplacian
% Within effect 2: frequency range - low vs. high
% Within effect 3: physiological condition - EO vs. EC
% Within effect 4: cortical location - FC2 (frontal) vs. POz (occipital)

specs = struct(...
    'vnames', {{'filter','ch','condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3_list', {{'FC2','POz'}},...
    'v4_list', {{'car','lap'}});

[rm_stim_full, nt_stim_full] = func_4waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 3-way models for each filtering scheme - CAR

% Fixed condition: CAR-filtered data
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC
% Within effect 3: cortical location - FC2 vs. POz

specs = struct(...
    'vnames', {{'ch','condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3_list', {{'FC2','POz'}},...
    'v4', {{'car'}});

[rm_stim_car, nt_stim_car] = func_3waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 3-way models for each filtering scheme - Small Laplacian

% Fixed condition: Small Laplacian-filtered data
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC
% Within effect 3: cortical location - FC2 vs. POz

specs = struct(...
    'vnames', {{'ch','condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3_list', {{'FC2','POz'}},...
    'v4', {{'lap'}});

[rm_stim_lap, nt_stim_lap] = func_3waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);



%% 3-way models for each filtering scheme - Raw EEG

% Fixed condition: Raw EEG data (no spatial filtering)
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC
% Within effect 3: cortical location - FC2 vs. POz

specs_3w_raw = struct(...
    'vnames', {{'ch','condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3_list', {{'FC2','POz'}},...
    'v4', {{'raw'}});

[rm_stim_raw, nt_stim_raw] = func_3waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 2-way models for each filtering, over Fc2 - CAR
% Fixed condition: CAR data over Fc2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC2'}},...
    'v4', {{'car'}});

[rm_stim_car_fc2, nt_stim_car_fc2, pw_stim_car_fc2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Fc2 - Small Laplacian
% Fixed condition: Laplacian data over Fc2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC2'}},...
    'v4', {{'lap'}});

[rm_stim_lap_fc2, nt_stim_lap_fc2, pw_stim_lap_fc2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Fc2 - RAW
% Fixed condition: raw data over Fc2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC2'}},...
    'v4', {{'raw'}});

[rm_stim_raw_fc2, nt_stim_raw_fc2, pw_stim_raw_fc2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 2-way models for each filtering, over Fc1 - CAR
% Fixed condition: CAR data over Fc1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC1'}},...
    'v4', {{'car'}});

[rm_stim_car_fc1, nt_stim_car_fc1, pw_stim_car_fc1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Fc1 - Small Laplacian
% Fixed condition: Laplacian data over Fc1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC1'}},...
    'v4', {{'lap'}});

[rm_stim_lap_fc1, nt_stim_lap_fc1, pw_stim_lap_fc1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Fc1 - RAW
% Fixed condition: raw data over Fc1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'FC1'}},...
    'v4', {{'raw'}});

[rm_stim_raw_fc1, nt_stim_raw_fc1, pw_stim_raw_fc1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Cp2 - CAR
% Fixed condition: CAR data over Cp2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP2'}},...
    'v4', {{'car'}});

[rm_stim_car_cp2, nt_stim_car_cp2, pw_stim_car_cp2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Cp2 - Small Laplacian
% Fixed condition: Laplacian data over Cp2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP2'}},...
    'v4', {{'lap'}});

[rm_stim_lap_cp2, nt_stim_lap_cp2, pw_stim_lap_cp2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Cp2 - RAW
% Fixed condition: raw data over Cp2
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP2'}},...
    'v4', {{'raw'}});

[rm_stim_raw_cp2, nt_stim_raw_cp2, pw_stim_raw_cp2] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 2-way models for each filtering, over Cp1 - CAR
% Fixed condition: CAR data over Cp1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP1'}},...
    'v4', {{'car'}});

[rm_stim_car_cp1, nt_stim_car_cp1, pw_stim_car_cp1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Cp1 - Small Laplacian
% Fixed condition: Laplacian data over Cp1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP1'}},...
    'v4', {{'lap'}});

[rm_stim_lap_cp1, nt_stim_lap_cp1, pw_stim_lap_cp1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over Cp1 - RAW
% Fixed condition: raw data over Cp1
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'CP1'}},...
    'v4', {{'raw'}});

[rm_stim_raw_cp1, nt_stim_raw_cp1, pw_stim_raw_cp1] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%% 2-way models for each filtering, over POz - CAR
% Fixed condition: CAR data over POz
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'POz'}},...
    'v4', {{'car'}});

[rm_stim_car_poz, nt_stim_car_poz, pw_stim_car_poz] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over POz - Small Laplacian
% Fixed condition: Laplacian data over POz
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'POz'}},...
    'v4', {{'lap'}});

[rm_stim_lap_poz, nt_stim_lap_poz, pw_stim_lap_poz] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);

%% 2-way models for each filtering, over POz - RAW
% Fixed condition: raw data over POz
% Within effect 1: frequency range - low vs. high
% Within effect 2: physiological condition - EO vs. EC

specs = struct(...
    'vnames', {{'condition','range'}},...
    'v1_list', {{'lo','hi'}},...
    'v2_list', {{'EO','EC'}},...
    'v3', {{'POz'}},...
    'v4', {{'raw'}});

[rm_stim_raw_poz, nt_stim_raw_poz, pw_stim_raw_poz] = func_2waymodel(data_mmspm_stimcap_eo, data_mmspm_stimcap_ec, specs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% display all results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-way model
disp(rm_stim_full)

%% 3-way model: CAR
disp(rm_stim_car)

%% 3-way model: Small Laplacian
disp(rm_stim_lap)

%% 2-way model: CAR, Fc2
disp(rm_stim_car_fc2)
disp(pw_stim_car_fc2)

%% 2-way model: CAR, POz
disp(rm_stim_car_poz)
disp(pw_stim_car_poz)

%% 2-way model: SL, Fc2
disp(rm_stim_lap_fc2)
disp(pw_stim_lap_fc2)

%% 2-way model: SL, POz
disp(rm_stim_lap_poz)
disp(pw_stim_lap_poz)

%% Change in oscillatory alpha power, Fc2 & POz
disp(pw_lao_fc2)

%% Correlation between oscillatory alpha power and slope, Fc2 & POz
disp(struct2table(corr_stat_lao_stimcap))

%% 2-way model confirmation: CAR, Fc1
disp(rm_stim_car_fc1)
disp(pw_stim_car_fc1)

%% 2-way model confirmation: SL, Fc1
disp(rm_stim_lap_fc1)
disp(pw_stim_lap_fc1)

%% 2-way model confirmation: CAR, Cp1
disp(rm_stim_car_cp1)
disp(pw_stim_car_cp1)

%% 2-way model confirmation: CAR, Cp2
disp(rm_stim_car_cp2)
disp(pw_stim_car_cp2)

%% 2-way model confirmation: SL, Cp1
disp(rm_stim_lap_cp1)
disp(pw_stim_lap_cp1)

%% 2-way model confirmation: SL, Cp2
disp(rm_stim_lap_cp2)
disp(pw_stim_lap_cp2)

%% broadband slope over POz
[~,p1,~,bb_stat_car_poz]=ttest(data_mmspm_stimcap_eo.bb_POz.car,data_mmspm_stimcap_ec.bb_POz.car);
disp(struct2table(bb_stat_car_poz))
disp(['p=' num2str(p1,'%.4f')])
disp(meanEffectSize(data_mmspm_stimcap_eo.bb_POz.car,data_mmspm_stimcap_ec.bb_POz.car,Effect='cohen'))

[~,p2,~,bb_stat_lap_poz]=ttest(data_mmspm_stimcap_eo.bb_POz.lap,data_mmspm_stimcap_ec.bb_POz.lap);
disp(struct2table(bb_stat_lap_poz))
disp(['p=' num2str(p2,'%.4f')])
disp(meanEffectSize(data_mmspm_stimcap_eo.bb_POz.lap,data_mmspm_stimcap_ec.bb_POz.lap,Effect='cohen'))

%% Correlation between oscillatory alpha power and broadband slope
disp(struct2table(corr_stat_lao_bb_stimcap))
