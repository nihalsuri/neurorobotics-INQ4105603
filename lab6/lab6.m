%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAB 06 — ERD/ERS on band power
%
%  PURPOSE:
%  --------
%  This script shows how raw EEG is transformed into usable features for
%  Brain–Machine Interfaces (BMI) — focusing on μ (10–12 Hz) and β
%  (18–24 Hz) rhythms from motor imagery (hands, feet, rest).
%
%  WHAT YOU'LL SEE:
%   1. Load & concatenate .GDF files (EEG + event markers)
%   2. (Optional) Apply spatial filtering: None, CAR, or Laplacian
%   3. Filter in μ/β bands and compute log band power per channel
%   4. Extract trials based on fixation–feedback periods
%   5. Plot:
%       - One example trial (raw, μ, β)
%       - Mean ± SE curves across all trials per class (hands/feet/rest)
%
%  HOW TO USE:
%   1. Edit the folder path below
%   2. Choose the spatial filter mode ('none', 'car', or 'lap')
%   3. Run section-by-section (Ctrl + Enter)
%
%  AUTHOR: Nihal Suri — Neurorobotics, University of Padova (2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath(genpath(pwd)); % ensure helper functions are in path

%% -------------------- SETUP --------------------
gdfFolder    = 'C:\Users\nihal\OneDrive\Documents\unipd\Semester3\neurorobotics-INQ4105603';
spatialMode  = 'lap';              % 'none' | 'car' | 'lap'
lapMaskFile  = 'laplacian16.mat';  % used only for 'lap'

[allS, fs, labels, EVENT, filesUsed] = concatGdfDropLast(gdfFolder);

fprintf('\nLoaded %d GDF file(s). Sampling rate = %.1f Hz, %d channels.\n', ...
        numel(filesUsed), fs, size(allS,2));
disp('Channel labels:'); disp(labels');

% Define event codes (standard BCI2000 conventions)
EV_FIX = hex2dec('0312');  % fixation ON
EV_FB  = hex2dec('030D');  % feedback ON
EV_HND = hex2dec('0305');  % both hands
EV_FT  = hex2dec('0303');  % both feet
EV_RST = hex2dec('030F');  % rest

%% -------------------- SPATIAL FILTER --------------------
Sx = applySpatialFilter(allS, spatialMode, lapMaskFile);
fprintf('Applied spatial filter: %s\n', upper(spatialMode));

%% -------------------- BANDPASS FILTERS --------------------
[bu_mu, au_mu]     = butter(5, [10 12]*2/fs, 'bandpass');
[bu_beta, au_beta] = butter(5, [18 24]*2/fs, 'bandpass');

%% -------------------- COMPUTE LOG BAND POWER --------------------
[muPow, betaPow] = computeBandPower(Sx, fs, bu_mu, au_mu, bu_beta, au_beta);

%% -------------------- EXTRACT TRIALS --------------------
trials = extractTrials(EVENT, fs, EV_FIX, EV_FB, EV_HND, EV_FT, EV_RST, size(allS, 1), 8);
fprintf('Extracted %d trials.\n', numel(trials));

%% -------------------- SELECT CHANNELS --------------------
wantedLabels = {'eeg:7','eeg:9','eeg:11'};
plotIdx = find(ismember(labels, wantedLabels));
if numel(plotIdx) < 3
    plotIdx = 1:min(3, size(allS,2)); % fallback
end
disp('Plotting channels:'); disp(labels(plotIdx));

%% -------------------- VISUALIZE A SINGLE TRIAL --------------------
resultsDir = fullfile(pwd, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

%% Visualization — single figure, 3 rows × 1 column (Raw, mu-band, beta-band)
if ~isempty(trials)
    k  = 2;                                        % choose trial index
    seg = trials(k).start : trials(k).stop;
    tt  = (seg - seg(1)) / fs;                     % time axis (s)

    % --- Channel selection ---
    if numel(plotIdx) > 3, plotIdx = plotIdx(1:3); end
    chNames = labels(plotIdx);                     % chNames has length == numel(plotIdx)

    % --- Create one figure with 3 stacked subplots ---
    figure('Name', sprintf('Trial %d (label=%d) — Raw, μ, β bands', k, trials(k).label));

    % ---------- Subplot 1: Raw ----------
    subplot(3,1,1); hold on;
    for c = 1:numel(plotIdx)
        plot(tt, allS(seg, plotIdx(c)), 'DisplayName', chNames{c});   
    end
    grid on; ylabel('Amplitude [\muV]');
    title(sprintf('Raw signal — Trial %d', k));
    legend('show','Location','best'); xlim([tt(1), tt(end)]);

    % ---------- Subplot 2: μ-band (log-power) ----------
    subplot(3,1,2); hold on;
    for c = 1:numel(plotIdx)
        plot(tt, muPow(seg, plotIdx(c)), 'DisplayName', chNames{c}); 
    end
    grid on; ylabel('power [dB]');
    title('Filtered \mu band (10–12 Hz)');
    legend('show','Location','best'); xlim([tt(1), tt(end)]);

    % ---------- Subplot 3: β-band (log-power) ----------
    subplot(3,1,3); hold on;
    for c = 1:numel(plotIdx)
        plot(tt, betaPow(seg, plotIdx(c)), 'DisplayName', chNames{c}); 
    end
    grid on; xlabel('Time [s]'); ylabel('power [dB]');
    title('Filtered \beta band (18–24 Hz)');
    legend('show','Location','best'); xlim([tt(1), tt(end)]);

    linkaxes(findall(gcf,'Type','axes'), 'x');
    sgtitle(sprintf('Trial %d (label=%d) — Spatial mode: %s', ...
        k, trials(k).label, upper(spatialMode)));
else
    warning('No trials to plot.');
end


%% -------------------- MEAN ± SEM PLOTS --------------------
classOrder = [1 2 0];                 % 1=hands, 2=feet, 0=rest
classNames = {'Hands','Feet','Rest'};
classCols  = [1 0 0; 0 0.6 0; 0 0 0];

winSec  = 8;
winSamp = round(winSec * fs);
t = (0:winSamp-1)/fs;

mean_sem_ts = @(TS, ch, L) local_mean_sem_ts(TS, ch, L, trials, winSamp);

% --- μ band ---
figure('Name','μ band — Mean ± SE');
for ci = 1:numel(plotIdx)
    subplot(3,1,ci); hold on; grid on;

    hMean = gobjects(1, numel(classOrder));  % store handles to mean lines
    for k = 1:numel(classOrder)
        [m,s] = mean_sem_ts(muPow, plotIdx(ci), classOrder(k));

        % mean (goes in legend)
        hMean(k) = plot(t, m, 'Color', classCols(k,:), ...
                        'LineWidth', 1.5, 'DisplayName', classNames{k});

        % ± SEM (kept out of legend)
        plot(t, m+s, ':', 'Color', classCols(k,:), 'HandleVisibility','off');
        plot(t, m-s, ':', 'Color', classCols(k,:), 'HandleVisibility','off');
    end

    title(sprintf('\\mu band — %s', labels{plotIdx(ci)}));
    legend(hMean, classNames, 'Location','best');   % legend by class
    xlim([0, winSec]);
    ylabel('power [dB]');
end
xlabel('Time [s]');
sgtitle(sprintf('\\mu band | Mean ± SE — %s', upper(spatialMode)));
saveas(gcf, fullfile(resultsDir, sprintf('MuBand_%s.png', spatialMode)));

% --- β band ---
figure('Name','β band — Mean ± SE');
for ci = 1:numel(plotIdx)
    subplot(3,1,ci); hold on; grid on;

    hMean = gobjects(1, numel(classOrder));
    for k = 1:numel(classOrder)
        [m,s] = mean_sem_ts(betaPow, plotIdx(ci), classOrder(k));

        hMean(k) = plot(t, m, 'Color', classCols(k,:), ...
                        'LineWidth', 1.5, 'DisplayName', classNames{k});

        plot(t, m+s, ':', 'Color', classCols(k,:), 'HandleVisibility','off');
        plot(t, m-s, ':', 'Color', classCols(k,:), 'HandleVisibility','off');
    end

    title(sprintf('\\beta band — %s', labels{plotIdx(ci)}));
    legend(hMean, classNames, 'Location','best');
    xlim([0, winSec]);
    ylabel('power [dB]');
end
xlabel('Time [s]');
sgtitle(sprintf('\\beta band | Mean ± SE — %s', upper(spatialMode)));
saveas(gcf, fullfile(resultsDir, sprintf('BetaBand_%s.png', spatialMode)));
disp('Figures saved to /results');