clear; clc; 

%% Provide gdf rootfolder
gdfFolder = 'C:\Users\nihal\OneDrive\Documents\unipd\Semester3\NR'; 

%% Run concatenation and show a quick summary
[allS, fs, labels, EVENT, filesUsed] = concatGdfDropLast(gdfFolder);

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

%% Filter design for mu (10–12 Hz) and beta (18–24 Hz)
muBand   = [10 12];      % Hz (narrow, as many labs specify)
betaBand = [18 24];      % Hz
butterOrder = 5;         % commonly used in teaching labs

Wmu   = muBand   * 2 / fs;   % normalized to Nyquist
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

%% Visualization — three figures (Raw, μ-band, β-band)
if ~isempty(trials)
    k  = 1;                                        % choose trial index
    seg = trials(k).start : trials(k).stop;
    tt  = (seg - seg(1)) / fs;                     % time axis (s)

    % --- Channel selection ---
    if numel(plotIdx) > 3, plotIdx = plotIdx(1:3); end
    chNames = cellstr(labels(plotIdx));

    % ---------- Figure 1: Raw ----------
    figure('Name', sprintf('Trial %d (label=%d) — Raw', k, trials(k).label));
    hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, allS(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    xlabel('Time [s]'); ylabel('Amplitude [\muV]');
    title(sprintf('Raw signal — Trial %d (%d samples) | Channels %s', ...
        k, numel(seg), strjoin(chNames, ', ')));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

    % ---------- Figure 2: mu-band filtered ----------
    figure('Name', sprintf('Trial %d (label=%d) — μ-band', k, trials(k).label));
    hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, muPow(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    xlabel('Time [s]'); ylabel('Amplitude [\muV]');
    title(sprintf('Filtered μ band (10–12 Hz) — Trial %d | Channels %s', ...
        k, strjoin(chNames, ', ')));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

    % ---------- Figure 3: beta-band filtered ----------
    figure('Name', sprintf('Trial %d (label=%d) — β-band', k, trials(k).label));
    hold on;
    for ii = 1:numel(plotIdx)
        plot(tt, betaPow(seg, plotIdx(ii)), 'DisplayName', sprintf('Channel %s', chNames{ii}));
    end
    hold off; grid on;
    xlabel('Time [s]'); ylabel('Amplitude [\muV]');
    title(sprintf('Filtered β band (18–24 Hz) — Trial %d | Channels %s', ...
        k, strjoin(chNames, ', ')));
    legend('show','Location','best');
    xlim([tt(1), tt(end)]);

else
    warning('No trials to plot.');
end


%% Average time-course plots per class (μ and β), one figure per channel

% Classes and display names
classOrder = [1 2 0];                      % 1=hands, 2=feet, 0=rest
classNames = {'both hands','both feet','rest'};

% Colors (match your mockup: hands=red, feet=green, rest=black)
classCols = [ ...
    1 0 0;   % hands - red
    0 0.6 0; % feet  - green
    0 0 0];  % rest  - black

% Window length aligned to fixation (seconds)
winSec  = 8;                          % change if your lab uses a different duration
winSamp = round(winSec * fs);
t = (0:winSamp-1)/fs;

% Ensure exactly three channels chosen
if numel(plotIdx) > 3, plotIdx = plotIdx(1:3); end
chNames = cellstr(labels(plotIdx));

% Helper to compute mean/SEM time series for a band TS and one channel
% Returns [1 x winSamp] mean and SEM across all trials for a given class label L
mean_sem_ts = @(TS, ch, L) local_mean_sem_ts(TS, ch, L, trials, winSamp);

% ==================== μ-band (one figure per channel) ====================
for ci = 1:numel(plotIdx)
    ch = plotIdx(ci);
    figure('Name', sprintf('mu band | Mean +/- SE | channel %s', chNames{ci}));
    hold on; grid on;

    for k = 1:numel(classOrder)
        L = classOrder(k);
        [m, s] = mean_sem_ts(muPow, ch, L);
        % mean
        plot(t, m, 'Color', classCols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', classNames{k});
        % +/- SEM as dotted
        plot(t, m + s, ':', 'Color', classCols(k,:));
        plot(t, m - s, ':', 'Color', classCols(k,:));
    end

    xlabel('Time [s]'); ylabel('power [dB]');
    title(sprintf('mu band | Mean +/- SE | channel %s', chNames{ci}));
    legend('Location','best');
    xlim([0, winSec]); ylim auto;
end

% ==================== β-band (one figure per channel) ====================
for ci = 1:numel(plotIdx)
    ch = plotIdx(ci);
    figure('Name', sprintf('beta band | Mean +/- SE | channel %s', chNames{ci}));
    hold on; grid on;

    for k = 1:numel(classOrder)
        L = classOrder(k);
        [m, s] = mean_sem_ts(betaPow, ch, L);
        % mean
        plot(t, m, 'Color', classCols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', classNames{k});
        % +/- SEM as dotted
        plot(t, m + s, ':', 'Color', classCols(k,:));
        plot(t, m - s, ':', 'Color', classCols(k,:));
    end

    xlabel('Time [s]'); ylabel('power [dB]');
    title(sprintf('beta band | Mean +/- SE | channel %s', chNames{ci}));
    legend('Location','best');
    xlim([0, winSec]); ylim auto;
end

% ----------------------- local helper (place below script code) -----------------------
function [m, s] = local_mean_sem_ts(TS, ch, L, trials, winSamp)
% TS:      time series (e.g., muPow or betaPow), size [N x C]
% ch:      channel index
% L:       class label to aggregate (1=hands, 2=feet, 0=rest)
% trials:  struct array with fields .start .stop .label
% winSamp: number of samples to take from fixation (t=0) forward
%
% Returns:
%   m: [1 x winSamp] mean across trials at each time point (NaN-aware)
%   s: [1 x winSamp] SEM  across trials at each time point (NaN-aware)

    trIdx = find([trials.label] == L);
    if isempty(trIdx)
        m = nan(1, winSamp);
        s = nan(1, winSamp);
        return;
    end

    M = nan(numel(trIdx), winSamp);  % trials x time

    for r = 1:numel(trIdx)
        st = trials(trIdx(r)).start;
        en = trials(trIdx(r)).stop;
        % segment from fixation for up to winSamp samples (cap by trial end)
        segEnd = min(st + winSamp - 1, en);
        len = segEnd - st + 1;
        if len > 0
            M(r,1:len) = TS(st:segEnd, ch);
        end
    end

    m = mean(M, 1, 'omitnan');
    % SEM with NaNs handled per time point
    n = sum(~isnan(M), 1);
    s = std(M, 0, 1, 'omitnan') ./ max(1, sqrt(n));
end













