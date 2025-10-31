%% ================================================================
%  This script follows the Lab03 brief for concatenation and label
%  vectors (Tk, Fk, Ak, CFk, Xk, Rk) and includes Lab02 essentials.
%  (Reference: Lab03 _ GDF data concatenation.pdf)
%  ---------------------------------------------------------------
%% ================================================================

clear; clc; close all;

%% ------------------- Paths & parameters (EDIT) ----------------
% Point these to your local folders
data_dir = fullfile('C:\Users\nihal\OneDrive\Documents\unipd\Semester3\neurorobotics-INQ4105603', 'data');                  

% List OFFLINE .gdf files for concatenation (Lab03 requirement)
gdf_list = dir(fullfile(data_dir, '*.gdf'));
assert(~isempty(gdf_list), 'No .gdf files found in data_dir.');

% Trial window relative to CUE onset (seconds). Adjust to your protocol.
winSec = [0, 4];   % from cue onset to +4 s

% Channel to visualize for examples/averages
vizChannelName = "eeg:7";   % change if your dataset has different labels

%% ------------------- Load & concatenate GDFs ------------------
% Concatenate EEG and EVENTS (with POS offsets). Also track run index Rk.

EEG_all      = [];               % concatenated EEG (samples x channels)
events_all   = struct('TYP',[],'POS',[],'DUR',[],'CHN',[]);
Fs_ref       = [];               % reference sampling rate
labels_ref   = [];               % reference channel label set
cumSamples   = 0;                % cumulative sample count for POS offset
run_boundaries = [];             % [startSample endSample] per run

fprintf('Loading and concatenating %d files ...\n', numel(gdf_list));
for f = 1:numel(gdf_list)
    fname = fullfile(gdf_list(f).folder, gdf_list(f).name);
    fprintf('  -> %s\n', gdf_list(f).name);
    [X, H] = sload(fname);  % BioSig loader (OK if it says slower M-function)

    % Basic consistency checks across files
    if isempty(Fs_ref), Fs_ref = H.SampleRate; else
        assert(abs(H.SampleRate - Fs_ref) < 1e-9, 'Sampling rate mismatch.');
    end
    if isempty(labels_ref), labels_ref = string(H.Label); else
        assert(isequal(string(H.Label), labels_ref), 'Channel labels mismatch.');
    end

    % Append EEG
    EEG_all = [EEG_all; X]; 
    nS = size(X,1);

    % Append events with position offset
    if isfield(H,'EVENT') && ~isempty(H.EVENT)
        TYP = double(H.EVENT.TYP(:));
        POS = double(H.EVENT.POS(:)) + cumSamples;  % offset POS
        DUR = double(H.EVENT.DUR(:));
        CHN = isfield(H.EVENT,'CHN');
        if CHN, CHN = double(H.EVENT.CHN(:)); else, CHN = zeros(size(TYP)); end

        events_all.TYP = [events_all.TYP; TYP];
        events_all.POS = [events_all.POS; POS];
        events_all.DUR = [events_all.DUR; DUR];
        events_all.CHN = [events_all.CHN; CHN];
    end

    % Track run (file) boundaries for Rk later
    run_boundaries = [run_boundaries; [cumSamples+1, cumSamples+nS]]; %#ok<AGROW>

    cumSamples = cumSamples + nS;
end

nSamples = size(EEG_all, 1);
nCh      = size(EEG_all, 2);
fprintf('Done. Total samples: %d | Channels: %d | Fs = %.2f Hz\n', nSamples, nCh, Fs_ref);

%% ------------------- Event summary (discover codes) ----------
% This prints a table of unique event codes, counts, and median duration.
% Use this to fill in the EVENT CODE MAP in the next section.

if ~isempty(events_all.TYP)
    TYP = events_all.TYP;
    POS = events_all.POS;
    DUR = events_all.DUR;

    [u,~,ix] = unique(TYP);
    medDur = accumarray(ix, DUR, [], @(x) median(x)./Fs_ref);
    cnt    = accumarray(ix, 1);

    fprintf('\nEvent summary (discover your codes):\n');
    disp(table(u, cnt, medDur, 'VariableNames', {'Code','Count','MedianDur_s'}));
else
    error('No events found in concatenated data. Cannot proceed.');
end

%% ------------------- EVENT CODE MAP----------------------
% IMPORTANT: After you read the summary above, set the code groups below.
% Example placeholders are commented; replace with your actual codes.
% If you are unsure, run once, inspect the table, edit, and run again.

codeMap = struct();
codeMap.FIX  = 786;                 % e.g., fixation / trial start
codeMap.CUE  = 773;               % e.g., cue classes (left/right)
codeMap.FEED = 781;                 % e.g., continuous feedback
codeMap.HIT  = 897;                    % if available
codeMap.MISS = 898;
codeMap.FEET = 771;
codeMap.START = 1;
codeMap.REST = 773;


% Minimal safety: ensure the fields exist (empty allowed)
defaults = {'FIX','CUE','FEED','HIT','MISS', 'REST', 'START', 'FEET'};
for k = 1:numel(defaults)
    if ~isfield(codeMap, defaults{k}), codeMap.(defaults{k}) = []; end
end

%% ------------------- Build label vectors (Tk, Fk, Ak, CFk, Xk, Rk)
% Lab03 requirement: create per-sample label vectors for concatenated data.
% Rk marks which run/file each sample belongs to (1..#files).

% Helper: make per-sample labels from event table + mapping
S = make_labels_from_events_conc(events_all, nSamples, codeMap);

% Build Rk (run index) from run boundaries
Rk = zeros(nSamples,1,'int16');
for r = 1:size(run_boundaries,1)
    Rk(run_boundaries(r,1):run_boundaries(r,2)) = r;
end

% Quick plot of label vectors
figure('Name','Label vectors','Color','w');
subplot(6,1,1); plot(S.Tk,'k');  ylabel('Tk');   title('Per-sample label vectors');
subplot(6,1,2); plot(S.Fk,'b');  ylabel('Fk');
subplot(6,1,3); plot(S.Ak,'m');  ylabel('Ak');
subplot(6,1,4); plot(S.CFk,'g'); ylabel('CFk');
subplot(6,1,5); plot(S.Xk,'r');  ylabel('Xk');
subplot(6,1,6); plot(Rk,'c');    ylabel('Rk');   xlabel('Samples');

%% ------------------- Trial extraction (using cues) ------------
% We epoch around cue onsets (rising edges of Ak != 0) and create:
%   - trials: [samples x channels x trials]
%   - Ck: one cue code per trial
% Adjust winSec above if your protocol is different.

preS  = round(winSec(1) * Fs_ref);
postS = round(winSec(2) * Fs_ref);
winN  = postS - preS;

Ak = S.Ak ~= 0;
cueOn = find(diff([0; Ak]) > 0);           % sample indices of cue onsets
cueOn = cueOn(cueOn + postS - 1 <= nSamples);  % ensure full window fits

nTrials = numel(cueOn);
trials  = zeros(winN, nCh, nTrials, 'double');
Ck      = zeros(nTrials,1,'int32');

for i = 1:nTrials
    a = cueOn(i) + preS;
    b = cueOn(i) + postS - 1;
    trials(:,:,i) = EEG_all(a:b, :);
    % Store the raw cue code at onset (keeps class identity)
    Ck(i) = S.Ak(cueOn(i));
end

fprintf('\nTrials built: %d | Window [%.2f %.2f] s | Samples/window: %d\n', ...
    nTrials, winSec(1), winSec(2), winN);

%% ------------------- Visualization ---------------------------
% a) One example trial per unique cue code, plotting selected channel
% b) Grand average per cue code for the selected channel

% Resolve viz channel index (fallback to 1 if not found)
chLabels = labels_ref;
chIdx = find(chLabels == vizChannelName, 1);
if isempty(chIdx), chIdx = 1; end

t_axis = (0:winN-1)/Fs_ref + winSec(1);

uCues = unique(Ck);
figure('Name','Example trials by cue','Color','w');
for u = 1:numel(uCues)
    idx = find(Ck == uCues(u), 1); % first trial of this cue
    if ~isempty(idx)
        subplot(numel(uCues),1,u);
        plot(t_axis, trials(:, chIdx, idx));
        grid on; ylabel('\muV');
        title(sprintf('Cue %d — example trial — %s', uCues(u), chLabels(chIdx)));
        if u == numel(uCues), xlabel('Time (s)'); end
    end
end

figure('Name','Grand averages by cue','Color','w');
for u = 1:numel(uCues)
    sel = find(Ck == uCues(u));
    if ~isempty(sel)
        ga = mean(trials(:, chIdx, sel), 3);
        subplot(numel(uCues),1,u);
        plot(t_axis, ga);
        grid on; ylabel('\muV');
        title(sprintf('Cue %d — grand average — %s', uCues(u), chLabels(chIdx)));
        if u == numel(uCues), xlabel('Time (s)'); end
    end
end


