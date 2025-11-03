%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAB 04 and 05 — MI BMI Logarithmic Band Power and Spatial Filters
%
%  GOAL:
%  -----
%  This script guides you through the complete processing chain for
%  EEG-based motor-imagery experiments, using GDF files and the BioSig /
%  EEGLAB toolboxes.  It is meant for first-time neurorobotics students
%  who want to understand, step-by-step, how raw EEG becomes usable
%  control features for a Brain–Machine Interface (BMI).
%
%  WHAT IT DOES:
%  --------------
%   1.  Recursively loads and concatenates all *.gdf files from a folder,
%       dropping the last channel (trigger line) automatically.
%   2.  Extracts event markers (fixation, cue, feedback) from the GDF
%       headers and builds trial windows:
%           [fixation ON  →  feedback ON  +  duration]
%       Each trial is labeled as:
%           1 = both hands imagery
%           2 = both feet imagery
%           0 = rest / idle
%   3.  Designs 5th-order Butterworth filters for the μ (10–12 Hz)
%       and β (18–24 Hz) rhythms — the two main motor-related frequency
%       bands of interest.
%   4.  Computes a time-series of logarithmic band-power for every
%       channel, sample-by-sample:
%           band-pass → square → 1 s moving average → log10(power)
%   5.  Visualizes data at three complementary levels:
%         • Within one trial  →  stacked plots of raw, μ, and β signals
%         • Across all trials →  mean ± SEM power curves per class
%           (hands / feet / rest) for three representative channels
%   6.  Demonstrates how feature extraction and averaging work before
%       classification (later labs will use these μ/β powers for
%       machine-learning decoders).
%
%  REQUIRED TOOLBOXES:
%  -------------------
%   • BioSig toolbox   – for SLOAD / GDF I/O
%   • EEGLAB toolbox   – for EEG utilities (optional, but recommended)
%
%  INPUT:
%   • Folder path containing one or more GDF recordings.
%
%  OUTPUT / FIGURES:
%   • Prints a summary of files, events, and extracted trials.
%   • Produces:
%        – 1 figure (3 rows) showing Raw, μ, β bands for a sample trial
%        – 1 figure (μ band)  with mean ± SEM time courses (3 subplots)
%        – 1 figure (β band)  with mean ± SEM time courses (3 subplots)
%
%  HOW TO USE:
%  -----------
%   1. Edit the variable  gdfFolder  below to point to your local data.
%   2. Run the script section-by-section (Ctrl+Enter) and read the comments.
%   3. Inspect the console output to verify event codes and trial counts.
%   4. Study the plots to see how the motor-imagery patterns appear in
%      μ and β bands across time and across channels.
%
%  NOTES:
%  ------
%   • All signals are in microvolts (raw) or log-power (μ/β).
%   • The code is intentionally verbose and modular for learning purposes.
%   • Once you understand each step, you can wrap it into reusable
%     functions for your own BMI pipelines.
%
%  AUTHOR:
%   NIHAL SURI — Neurorobotics, University of Padova, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; 

%% Provide gdf rootfolder
gdfFolder = 'C:\Users\nihal\OneDrive\Documents\unipd\Semester3\neurorobotics-INQ4105603'; 

%% Run concatenation and show a quick summary
[allS, fs, labels, EVENT, filesUsed] = concatGdfDropLast(gdfFolder);

%% Spatial filtering mode: 'none' | 'car' | 'lap'
spatialMode = 'none';         % <— change to 'car' or 'lap' to switch
lapMaskFile = 'laplacian16.mat';   % your provided mask in working directory

%% For basic info 
fprintf('Files used (%d):\n', numel(filesUsed));
for i = 1:numel(filesUsed), fprintf('  %s\n', filesUsed{i}); end

fprintf('Concatenated shape: %d samples x %d channels (after dropping last channel)\n', ...
    size(allS,1), size(allS,2));
fprintf('Fs = %.1f Hz\n', fs);

disp('Channel labels (after drop):'); disp(labels(:)');

printEventHex(EVENT);

%% Event codes and quick inspection

EV_TRIAL_START = hex2dec('0001');  %# not strictly needed here
EV_FEET        = hex2dec('0303');  % "both feet"
EV_HANDS       = hex2dec('0305');  % "both hands"
EV_FEEDBACK    = hex2dec('030D');  % continuous feedback ON  (no OFF present)
EV_REST        = hex2dec('030F');  % rest (if used)
EV_FIXATION    = hex2dec('0312');  % fixation cross ON


% show which of these exist in your data
present = unique(EVENT.TYP)';
fprintf('Present event types (hex): ');
fprintf('0x%04X ', present);
fprintf('\n');

%% Apply spatial filter to raw EEG (before band-power pipeline)
% Input: allS  [samples x channels]
% Output: Sx   [samples x channels] (spatially filtered)

Sx = allS;  % default (no spatial filtering)

switch lower(spatialMode)
    case 'none'
        % No change; baseline condition for comparison
        fprintf('Spatial mode: NONE (baseline)\n');

    case 'car'
        % Common Average Reference: subtract per-sample channel mean
        % CAR works row-wise: for each time sample, subtract the mean across channels
        m = mean(allS, 2, 'omitnan');      % [samples x 1]
        Sx = allS - m;                     % implicit broadcast to all channels
        fprintf('Spatial mode: CAR (common average reference)\n');

    case 'lap'
        % Laplacian: multiply by a spatial mask (channels x channels)
        fprintf('Spatial mode: LAPLACIAN (using %s)\n', lapMaskFile);

        % Load the mask; try to find a square matrix in the MAT-file
        L = load(lapMaskFile);
        % Heuristic: pick the first square numeric field as the mask
        mask = [];
        fn = fieldnames(L);
        for i = 1:numel(fn)
            val = L.(fn{i});
            if isnumeric(val) && ismatrix(val) && size(val,1)==size(val,2)
                mask = val; break;
            end
        end
        assert(~isempty(mask), 'No square numeric mask found in %s', lapMaskFile);

        % Sanity checks on dimensions & channel order
        assert(size(mask,1) == size(allS,2), ...
            'Laplacian mask size (%d) does not match #channels (%d)', ...
            size(mask,1), size(allS,2));

        % If the mask expects a different channel order than 'labels', reorder here.
        % (Most course masks match the recording montage after dropping trigger.)

        % Apply Laplacian: each sample row times mask → filtered channels
        % If your mask is defined as channels×channels *column* mixing, use Sx = allS * mask.
        % If defined as *row* mixing, use Sx = mask * allS' then transpose.
        % Most Laplacian masks are applied as right-multiply:
        Sx = allS * mask;

    otherwise
        error('Unknown spatialMode: %s', spatialMode);
end


%% Filter design for mu (10–12 Hz) and beta (18–24 Hz)
muBand   = [10 12];      % Hz (narrow, as many labs specify)
betaBand = [18 24];      % Hz
butterOrder = 5;         % 

Wmu   = muBand   * 2 / fs;   % normalized to Nyquist Freq
Wbeta = betaBand * 2 / fs;

[bu_mu,   au_mu]   = butter(butterOrder, Wmu,   'bandpass');
[bu_beta, au_beta] = butter(butterOrder, Wbeta, 'bandpass');

% optional visual check:
% fvtool(bu_mu, au_mu, 'Fs', fs);  % inspect frequency response
% fvtool(bu_beta, au_beta, 'Fs', fs);


%% Compute mu/beta logarithmic band power time series
avgWinSec = 1.0;                 % 1-second moving average
avgN = round(avgWinSec * fs);
avgKernel = ones(avgN,1) / avgN; % uniform window
useLog10 = true;                 % log base 10

muPow   = zeros(size(allS));
betaPow = zeros(size(allS));

for ch = 1:size(allS,2)
    x = allS(:,ch);

    % band-pass (zero-phase)
    x_mu   = filtfilt(bu_mu,   au_mu,   x);
    x_beta = filtfilt(bu_beta, au_beta, x);

    % instantaneous power
    x_mu2   = x_mu.^2;
    x_beta2 = x_beta.^2;

    % 1-s moving average (causal smoothing is fine for feature building)
    x_mu_avg   = filter(avgKernel, 1, x_mu2);
    x_beta_avg = filter(avgKernel, 1, x_beta2);

    % log transform to stabilize variance / compress scale
    if useLog10
        muPow(:,ch)   = log10(x_mu_avg + eps);
        betaPow(:,ch) = log10(x_beta_avg + eps);
    else
        muPow(:,ch)   = log(x_mu_avg + eps);
        betaPow(:,ch) = log(x_beta_avg + eps);
    end
end

%% Extract trials (FIXATION ON -> FEEDBACK ON + DUR), label by cue
% Uses DUR from FEEDBACK ON (0x030D). If DUR missing/zero, falls back to
% next FIXATION or a capped window so you still get trials.


% --- Pull and sanitize event vectors ---
TYP = EVENT.TYP(:);
POS = EVENT.POS(:);
DUR = EVENT.DUR(:);
if isempty(DUR), DUR = zeros(size(POS)); end

% --- Ensure events are in chronological order ---
if ~issorted(POS)
    [POS, ord] = sort(POS);
    TYP = TYP(ord);
    DUR = DUR(ord);
end

% --- Indices of the key on-events we care about ---
idx_fix_on = find(TYP == EV_FIXATION);
idx_fb_on  = find(TYP == EV_FEEDBACK);

% Helper: first index in array 'posArray' whose POS > t
first_after = @(posArray, t) find(posArray > t, 1, 'first');

% Fallback cap when no clear stop can be inferred (seconds)
maxWinSec = 8;
maxWinSamp = round(maxWinSec * fs);

trials = struct('start',{},'stop',{},'label',{});
tcount = 0;

for ii = 1:numel(idx_fix_on)
    t_start = POS(idx_fix_on(ii));

    % 1) primary stop = first FEEDBACK ON after fixation + its DUR
    j_fb = first_after(POS(idx_fb_on), t_start);
    if ~isempty(j_fb)
        fb_ev_idx = idx_fb_on(j_fb);
        fb_on_pos = POS(fb_ev_idx);
        fb_on_dur = DUR(fb_ev_idx);

        if fb_on_dur > 0
            t_stop = fb_on_pos + fb_on_dur;
        else
            % 2) no DUR → try to bound by the *next* fixation (as a safe end)
            j_fix_next = first_after(POS(idx_fix_on), t_start);
            if ~isempty(j_fix_next) && j_fix_next < numel(idx_fix_on)
                t_stop = POS(idx_fix_on(j_fix_next+1)) - 1;
            else
                % 3) final fallback: cap window length
                t_stop = min(t_start + maxWinSamp, size(allS,1));
            end
        end
    else
        % no feedback after this fixation → next fixation or cap
        j_fix_next = first_after(POS(idx_fix_on), t_start);
        if ~isempty(j_fix_next) && j_fix_next < numel(idx_fix_on)
            t_stop = POS(idx_fix_on(j_fix_next+1)) - 1;
        else
            t_stop = min(t_start + maxWinSamp, size(allS,1));
        end
    end

    % Guards
    if ~(t_stop > t_start)
        continue;
    end
    t_stop = min(t_stop, size(allS,1));  % don’t exceed signal length

    % --- Determine class label from cues inside [t_start, t_stop] ---
    inwin = (POS >= t_start) & (POS <= t_stop);
    types_in = TYP(inwin);

    % Priority: HANDS > FEET > REST (tweak if your lab says otherwise)
    label = NaN;          % 1=hands, 2=feet, 0=rest
    if any(types_in == EV_HANDS)
        label = 1;
    elseif any(types_in == EV_FEET)
        label = 2;
    elseif any(types_in == EV_REST)
        label = 0;
    else
        % No recognizable cue → skip this trial
        continue;
    end

    % Record
    tcount = tcount + 1;
    trials(tcount).start = t_start;
    trials(tcount).stop  = t_stop;
    trials(tcount).label = label;
end

fprintf('Trials extracted (fixation->feedback+dur): %d\n', numel(trials));

%% Pick 3 channels to visualize 
wantedLabels = {'eeg:3','eeg:6','eeg:8'};
plotIdx = [];

for nm = wantedLabels
    j = find(labels == nm, 1);
    if ~isempty(j), plotIdx(end+1) = j; end 
end
if numel(plotIdx) < 3
    plotIdx = 1:min(3, size(allS,2));  % fallback to first three
end
disp('Plotting channels:');
disp(labels(plotIdx));

%% Visualization — single figure, 3 rows × 1 column (Raw, mu-band, beta-band)
if ~isempty(trials)
    k  = 1;                                        % choose trial index
    seg = trials(k).start : trials(k).stop;
    tt  = (seg - seg(1)) / fs;                     % time axis (s)

    % --- Channel selection ---
    if numel(plotIdx) > 3, plotIdx = plotIdx(1:3); end
    chNames = cellstr(labels(plotIdx));

    % --- Create one figure with 3 stacked subplots ---
    figure('Name', sprintf('Trial %d (label=%d) — Raw, μ, β bands', k, trials(k).label));

    % ---------- Subplot 1: Raw ----------
    subplot(3,1,1); hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, allS(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    ylabel('Amplitude [\muV]');
    title(sprintf('Raw signal — Trial %d (%d samples) | Channels %s', ...
        k, numel(seg), strjoin(chNames, ', ')));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

    % ---------- Subplot 2: mu-band filtered ----------
    subplot(3,1,2); hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, muPow(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    ylabel('Amplitude [\muV]');
    title(sprintf('Filtered μ band (10–12 Hz) — Trial %d', k));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

    % ---------- Subplot 3: beta-band filtered ----------
    subplot(3,1,3); hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, betaPow(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    xlabel('Time [s]'); ylabel('Amplitude [\muV]');
    title(sprintf('Filtered β band (18–24 Hz) — Trial %d', k));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

    % Optional: link x-axes for synchronized zoom/pan
    linkaxes(findall(gcf,'Type','axes'), 'x');
    sgtitle(sprintf('Trial %d (label = %d): Raw, μ and β bands', k, trials(k).label));
else
    warning('No trials to plot.');
end



%% Average time-course plots per class (mu and beta), one figure with 3 subplots

% Classes and display names
classOrder = [1 2 0];                      % 1=hands, 2=feet, 0=rest
classNames = {'both hands','both feet','rest'};

% Colors (match your mockup: hands=red, feet=green, rest=black)
classCols = [ ...
    1 0 0;   % hands - red
    0 0.6 0; % feet  - green
    0 0 0];  % rest  - black

% Window length aligned to fixation (seconds)
winSec  = 8;                          
winSamp = round(winSec * fs);
t = (0:winSamp-1)/fs;

% Ensure exactly three channels chosen
if numel(plotIdx) > 3, plotIdx = plotIdx(1:3); end
chNames = cellstr(labels(plotIdx));

% Helper for mean/SEM time series
mean_sem_ts = @(TS, ch, L) local_mean_sem_ts(TS, ch, L, trials, winSamp);

%% ==================== mu-band figure ====================
figure('Name','μ band | Mean ± SE | all channels');
for ci = 1:numel(plotIdx)
    ch = plotIdx(ci);
    subplot(numel(plotIdx), 1, ci); hold on; grid on;

    for k = 1:numel(classOrder)
        L = classOrder(k);
        [m, s] = mean_sem_ts(muPow, ch, L);
        % mean
        plot(t, m, 'Color', classCols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', classNames{k});
        % ±SEM as dotted
        plot(t, m + s, ':', 'Color', classCols(k,:));
        plot(t, m - s, ':', 'Color', classCols(k,:));
    end

    xlabel('Time [s]'); ylabel('power [dB]');
    title(sprintf('\\mu band | %s', chNames{ci}));
    legend('Location','best');
    xlim([0, winSec]); ylim auto;
end
sgtitle('\mu band | Mean ± SE across classes and channels');

%% ==================== beta-band figure ====================
figure('Name','β band | Mean ± SE | all channels');
for ci = 1:numel(plotIdx)
    ch = plotIdx(ci);
    subplot(numel(plotIdx), 1, ci); hold on; grid on;

    for k = 1:numel(classOrder)
        L = classOrder(k);
        [m, s] = mean_sem_ts(betaPow, ch, L);
        % mean
        plot(t, m, 'Color', classCols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', classNames{k});
        % ±SEM as dotted
        plot(t, m + s, ':', 'Color', classCols(k,:));
        plot(t, m - s, ':', 'Color', classCols(k,:));
    end

    xlabel('Time [s]'); ylabel('power [dB]');
    title(sprintf('\\beta band | %s', chNames{ci}));
    legend('Location','best');
    xlim([0, winSec]); ylim auto;
end
sgtitle('\beta band | Mean ± SE across classes and channels');
















