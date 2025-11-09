%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAB 06 — ERD/ERS on band power
%
%  PURPOSE:
%  --------
%  This script shows how raw EEG is transformed into usable features for
%  Brain–Machine Interfaces (BMI) — focusing on μ (10–12 Hz) and β
%  (18–24 Hz) rhythms from motor imagery (hands, feet).
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

%% -------------------- PART 1 --------------------
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

%% -------------------- SPATIAL FILTER --------------------
Sx = applySpatialFilter(allS, spatialMode, lapMaskFile);
fprintf('Applied spatial filter: %s\n', upper(spatialMode));

%% -------------------- BANDPASS FILTERS --------------------
[bu_mu, au_mu]     = butter(5, [10 12]*2/fs, 'bandpass');
[bu_beta, au_beta] = butter(5, [18 24]*2/fs, 'bandpass');

%% -------------------- COMPUTE LOG BAND POWER --------------------
[muPow, betaPow] = computeBandPower(Sx, fs, bu_mu, au_mu, bu_beta, au_beta);

%% -------------------- EXTRACT TRIALS --------------------
trials = extractTrials(EVENT, fs, EV_FIX, EV_FB, EV_HND, EV_FT, size(allS, 1), 8);
fprintf('Extracted %d trials.\n', numel(trials));

%% -------------------- PART 2: ERD/ERS --------------------
winSec = 8;                   % fixed window length (s) from fixation
T      = round(winSec*fs);    % samples
C      = size(muPow,2);
K      = numel(trials);

TYP = EVENT.TYP(:); 
POS = EVENT.POS(:);

% Preallocate tensors
X_mu   = nan(T, C, K);
X_beta = nan(T, C, K);
FixMask = false(T, K);    % fixation (reference) samples per trial
ActMask = false(T, K);    % activity (feedback) samples per trial

for k = 1:K
    st  = trials(k).start;
    en  = min(trials(k).stop, st + T - 1);
    len = en - st + 1;

    if len > 0
        X_mu(1:len,:,k)   = muPow(st:en,   :);
        X_beta(1:len,:,k) = betaPow(st:en, :);
    end

    % --- find first FEEDBACK ON after this fixation ---
    idx_fb = find( (TYP == EV_FB) & (POS > st), 1, 'first' );
    if isempty(idx_fb)
        fb_on = NaN;
    else
        fb_on = POS(idx_fb);
    end

    % --- build fixation vs activity masks in 1..T ---
    if ~isnan(fb_on)
        fb_idx = max(1, min(T, fb_on - st + 1));   % index within [1..T]
        FixMask(1:max(1,fb_idx-1), k) = true;      % [1 .. fb_idx-1]
        ActMask(fb_idx:min(T,len), k)  = true;     % [fb_idx .. len]
    else
        % if no feedback marker, split the available window in half
        half = max(1, floor(len/2));
        FixMask(1:half, k)      = true;
        ActMask(half+1:len, k)  = true;
    end
end

% ---- ERD/ERS computation: 100 * (Trial - RefMean)/RefMean ----
% reference mean per (channel x trial)
FixMean_mu   = nan(1, C, K);
FixMean_beta = nan(1, C, K);

for k = 1:K
    ref_rows = FixMask(:,k);
    FixMean_mu(1,:,k)   = mean(X_mu  (ref_rows, :, k), 1, 'omitnan');
    FixMean_beta(1,:,k) = mean(X_beta(ref_rows, :, k), 1, 'omitnan');
end

Ref_mu   = repmat(FixMean_mu,   [T 1 1]);   % T x C x K
Ref_beta = repmat(FixMean_beta, [T 1 1]);

ERD_mu   = 100 * (X_mu   - Ref_mu)   ./ Ref_mu;
ERD_beta = 100 * (X_beta - Ref_beta) ./ Ref_beta;


%% -------------------- PART 3: Temporal visualization (mean ± SEM) + feedback marker --------------------
% Pick one channel to visualize - eeg:7, eeg:11 - useful
wantedLabels = {'eeg:9'};                         
chIdx = find(ismember(labels, wantedLabels), 1); 
if isempty(chIdx), chIdx = 1; end

% Class trial indices
handsIdx = find([trials.label] == 1);                 % 1 = Hands
feetIdx  = find([trials.label] == 2);                 % 2 = Feet

% Time axis (relative to fixation)
t   = (0:T-1) / fs;

% --- Compute the feedback onset time (in seconds) for each trial, relative to fixation ---
TYP = EVENT.TYP(:); POS = EVENT.POS(:);
fb_on_abs = nan(1, numel(trials));                    % absolute sample index of FB ON for each trial
for k = 1:numel(trials)
    st = trials(k).start;                             % Fixation ON (absolute sample)
    idx_fb = find((TYP==EV_FB) & (POS > st), 1, 'first');  % strictly after fixation
    if ~isempty(idx_fb), fb_on_abs(k) = POS(idx_fb); end
end

% Convert to "seconds from fixation" for each trial (NaN if missing)
fb_t_rel = (fb_on_abs - [trials.start]) / fs;         % vector length K with NaNs where no FB

% Robust per-class summary of feedback onset (median)
fb_t_hands = median(fb_t_rel(handsIdx), 'omitnan');   % seconds
fb_t_feet  = median(fb_t_rel(feetIdx ), 'omitnan');   % seconds

% --- helper: SEM along trials for a time course (ignores NaNs) ---
sem = @(X,dim) std(X,0,dim,'omitnan') ./ max(1, sqrt(sum(~isnan(X),dim)));

% Helper to get mean±SEM time course for a given ERD tensor and trial set
get_curve = @(ERD, ch, idx) deal( ...
    mean(ERD(:, ch, idx), 3, 'omitnan'), ...
    sem (ERD(:, ch, idx), 3) );

%% ---- μ band plot ----
figure('Name','ERD/ERS — μ band'); hold on; grid on;

[mH, sH] = get_curve(ERD_mu, chIdx, handsIdx);
[mF, sF] = get_curve(ERD_mu, chIdx, feetIdx);

p1 = plot(t, mH, 'r', 'LineWidth', 1.5, 'DisplayName','Hands');
plot(t, mH+sH, 'r:', 'HandleVisibility','off'); 
plot(t, mH-sH, 'r:', 'HandleVisibility','off');

p2 = plot(t, mF, 'Color',[0 0.6 0], 'LineWidth', 1.5, 'DisplayName','Feet');
plot(t, mF+sF, 'Color',[0 0.6 0], 'LineStyle',':', 'HandleVisibility','off'); 
plot(t, mF-sF, 'Color',[0 0.6 0], 'LineStyle',':', 'HandleVisibility','off');

% --- vertical feedback markers (per class) ---
vh = gobjects(0);
if ~isnan(fb_t_hands)
    vh(end+1) = xline(fb_t_hands, '--r', 'DisplayName','FB start (Hands)', ...
        'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
end
if ~isnan(fb_t_feet)
    vh(end+1) = xline(fb_t_feet,  '--', 'Color',[0 0.6 0], 'DisplayName','FB start (Feet)', ...
        'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
end

xlabel('Time (s) from fixation'); ylabel('ERD/ERS (%)');
title(sprintf('ERD/ERS — \\mu band @ %s', labels{chIdx}));
legend([p1 p2 vh(:)'], 'Location','best');    % include lines in legend
xlim([t(1), t(end)]);

%% ---- β band plot ----
figure('Name','ERD/ERS — β band'); hold on; grid on;

[mH, sH] = get_curve(ERD_beta, chIdx, handsIdx);
[mF, sF] = get_curve(ERD_beta, chIdx, feetIdx);

p1 = plot(t, mH, 'r', 'LineWidth', 1.5, 'DisplayName','Hands');
plot(t, mH+sH, 'r:', 'HandleVisibility','off'); 
plot(t, mH-sH, 'r:', 'HandleVisibility','off');

p2 = plot(t, mF, 'Color',[0 0.6 0], 'LineWidth', 1.5, 'DisplayName','Feet');
plot(t, mF+sF, 'Color',[0 0.6 0], 'LineStyle',':', 'HandleVisibility','off'); 
plot(t, mF-sF, 'Color',[0 0.6 0], 'LineStyle',':', 'HandleVisibility','off');

vh = gobjects(0);
if ~isnan(fb_t_hands)
    vh(end+1) = xline(fb_t_hands, '--r', 'DisplayName','FB start (Hands)', ...
        'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
end
if ~isnan(fb_t_feet)
    vh(end+1) = xline(fb_t_feet,  '--', 'Color',[0 0.6 0], 'DisplayName','FB start (Feet)', ...
        'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
end

xlabel('Time (s) from fixation'); ylabel('ERD/ERS (%)');
title(sprintf('ERD/ERS — \\beta band @ %s', labels{chIdx}));
legend([p1 p2 vh(:)'], 'Location','best');
xlim([t(1), t(end)]);


%% -------------------- PART 4: Spatial visualization (topoplots) --------------------
% Requires EEGLAB's topoplot and chanlocs16.mat (electrode positions)
load('chanlocs16.mat');   % gives you chanlocs16

% We’ll compute the *average ERD* across time windows and trials, then plot scalp maps.
% Periods (indices in 1..T):
FixPeriod = false(T,1);  % we’ll define as "where most trials are in fixation"
ActPeriod = false(T,1);  % and "where most trials are in activity"

% Build "majority" masks so periods are common across trials (nice for comparability)
fixCount = sum(FixMask, 2);   % how many trials consider each time sample as fixation
actCount = sum(ActMask, 2);   % how many trials consider each time sample as activity
FixPeriod(fixCount >= round(0.5*K)) = true;   % samples that are fixation in ≥50% trials
ActPeriod(actCount >= round(0.5*K)) = true;   % samples that are activity in ≥50% trials

% Class logicals
isHands = ([trials.label] == 1);
isFeet  = ([trials.label] == 2);

% Helper: average ERD over a time window and across selected trials → [1 x C] vector
avg_erd_topo = @(ERD, timeMask, trialMask) squeeze( ...
    mean( mean( ERD(timeMask, :, trialMask), 1, 'omitnan' ), 3, 'omitnan') ).';

% ---- μ band topoplots ----
ERD_Ref_Hands_mu = avg_erd_topo(ERD_mu,   FixPeriod, isHands);   % [1 x C]
ERD_Act_Hands_mu = avg_erd_topo(ERD_mu,   ActPeriod, isHands);
ERD_Ref_Feet_mu  = avg_erd_topo(ERD_mu,   FixPeriod, isFeet);
ERD_Act_Feet_mu  = avg_erd_topo(ERD_mu,   ActPeriod, isFeet);

figure('Name','Topoplots — μ band ERD/ERS');
subplot(2,2,1); topoplot(ERD_Ref_Hands_mu, chanlocs16); title('μ: Hands — Reference');
subplot(2,2,2); topoplot(ERD_Act_Hands_mu, chanlocs16); title('μ: Hands — Activity');
subplot(2,2,3); topoplot(ERD_Ref_Feet_mu,  chanlocs16); title('μ: Feet — Reference');
subplot(2,2,4); topoplot(ERD_Act_Feet_mu,  chanlocs16); title('μ: Feet — Activity');
colorbar;

% ---- β band topoplots ----
ERD_Ref_Hands_bt = avg_erd_topo(ERD_beta, FixPeriod, isHands);
ERD_Act_Hands_bt = avg_erd_topo(ERD_beta, ActPeriod, isHands);
ERD_Ref_Feet_bt  = avg_erd_topo(ERD_beta, FixPeriod, isFeet);
ERD_Act_Feet_bt  = avg_erd_topo(ERD_beta, ActPeriod, isFeet);

figure('Name','Topoplots — β band ERD/ERS');
subplot(2,2,1); topoplot(ERD_Ref_Hands_bt, chanlocs16); title('β: Hands — Reference');
subplot(2,2,2); topoplot(ERD_Act_Hands_bt, chanlocs16); title('β: Hands — Activity');
subplot(2,2,3); topoplot(ERD_Ref_Feet_bt,  chanlocs16); title('β: Feet — Reference');
subplot(2,2,4); topoplot(ERD_Act_Feet_bt,  chanlocs16); title('β: Feet — Activity');
colorbar;





