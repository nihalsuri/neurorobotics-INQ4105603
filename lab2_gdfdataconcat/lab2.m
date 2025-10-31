%% Load EEG Data
clear;

% Adjust the path/filename - ACCORDING TO USER
gdf_file = fullfile('C:\Users\nihal\OneDrive\Documents\unipd\Semester3\neurorobotics-INQ4105603\data','ah7.20170613.161402.offline.mi.mi_bhbf.gdf');

% BioSig sload returns [data, header]
[EEG, HDR] = sload(gdf_file);

Fs   = HDR.SampleRate;          % samples per second (Hz)
[chN, chLabels] = deal(HDR.NS, string(HDR.Label));  % number of channels, labels
fprintf('Loaded %d channels at %.2f Hz\n', chN, Fs);

%% Event table (if present)
if isfield(HDR,'EVENT')
    TYP = double(HDR.EVENT.TYP(:));
    POS = double(HDR.EVENT.POS(:));
    DUR = double(HDR.EVENT.DUR(:));
    if isfield(HDR.EVENT,'CHN'), CHN = double(HDR.EVENT.CHN(:)); else, CHN = zeros(size(TYP)); end

    % Quick summary of unique event codes and counts
    [uT,~,idx] = unique(TYP);
    counts = accumarray(idx,1);
    disp(table(uT,counts,'VariableNames',{'EventCode','Count'}));
else
    error('No events found in HDR.EVENT — but they are required for the lab steps.');
end

%% First Plots
t = (0:size(EEG,1)-1)/Fs;  % time vector (seconds)

% Choose a channel by label or index:
chanName = "eeg:3";            
chanIdx  = find(chLabels==chanName,1);
if isempty(chanIdx), chanIdx = 1; end    % fallback

winSec   = 5;
i0       = 1;                            % start at 0 s
i1       = i0 + round(winSec*Fs)-1;      % 5s worth of samples

figure; 
plot(t(i0:i1), EEG(i0:i1, chanIdx));
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title(sprintf('Channel %s — first %d s', chLabels(chanIdx), winSec)); grid on;

% Now 3 channels with same y-limits
chanSet  = ["eeg:4","eeg:5","eeg:6"];            
idxSet   = arrayfun(@(nm)find(chLabels==nm,1), chanSet);
idxSet(isnan(idxSet)) = [];              % drop missing
Y = EEG(i0:i1, idxSet);

yl = [min(Y,[],'all') max(Y,[],'all')];  % common scale

figure;
for k = 1:numel(idxSet)
    subplot(numel(idxSet),1,k);
    plot(t(i0:i1), EEG(i0:i1, idxSet(k)));
    ylabel('\muV'); grid on;
    title(sprintf('%s (common scale)', chLabels(idxSet(k))));
    ylim(yl + 0.05*[-1 1]*range(yl));
end
xlabel('Time (s)');

 

%% Summarize codes + median duration (s)
TYP = double(HDR.EVENT.TYP(:));
POS = double(HDR.EVENT.POS(:));
DUR = double(HDR.EVENT.DUR(:));

[u,~,ix] = unique(TYP);
medDur = accumarray(ix, DUR, [], @(x)median(x)/Fs);
cnt    = accumarray(ix, 1);

disp(table(u, cnt, medDur, 'VariableNames',{'Code','Count','MedianDur_s'}));
if isfield(HDR.EVENT,'Desc')
    try
        D = string(HDR.EVENT.Desc(:));
        % Show first occurrence per code
        firstIdx = accumarray(ix, (1:numel(ix))', [], @(x) x(1));
        disp(table(u, D(firstIdx), 'VariableNames',{'Code','ExampleDesc'}));
    catch
        disp('Could not display Desc cleanly.');
    end
end


%% Added from summary 
codeMap.FIX  = [786];             % fixation code(s)
codeMap.CUE  = [773];             % cue codes 
codeMap.FEED = [781];             % continuous feedback
codeMap.HIT  = [897];             % hit
codeMap.MISS = [898];             % miss

S = make_labels_from_events(HDR, size(EEG,1), Fs, codeMap);

% Quick plot of label vectors (downsampled for speed if needed)
figure; 
subplot(5,1,1); plot(S.Tk,'k');    ylabel('Tk');   title('Label vectors');
subplot(5,1,2); plot(S.Fk,'b');    ylabel('Fk');
subplot(5,1,3); plot(S.Ak,'m');    ylabel('Ak');
subplot(5,1,4); plot(S.CFk,'g');   ylabel('CFk');
subplot(5,1,5); plot(S.Xk,'r');    ylabel('Xk'); xlabel('Samples');


%% Concatenate 
dataDir = 'C:\Users\nihal\OneDrive\Documents\unipd\Semester3\NR\data';       % folder with all offline .gdf files
flist   = dir(fullfile(dataDir, '*.gdf'));

EEG_all = [];
events  = struct('TYP',[],'POS',[],'DUR',[],'CHN',[]);
cumN    = 0;   % cumulative samples

Fs_ref  = [];
labels_ref = [];

for f = 1:numel(flist)
    fname = fullfile(flist(f).folder, flist(f).name);
    [X, H] = sload(fname);

    % Consistency checks
    if isempty(Fs_ref), Fs_ref = H.SampleRate; else
        assert(abs(H.SampleRate - Fs_ref)<1e-6,'Inconsistent Fs across files');
    end
    if isempty(labels_ref), labels_ref = string(H.Label); else
        assert(isequal(string(H.Label), labels_ref), 'Channel labels mismatch');
    end

    EEG_all = [EEG_all; X];  
    nS = size(X,1);

    if isfield(H,'EVENT')
        TYP = double(H.EVENT.TYP(:));
        POS = double(H.EVENT.POS(:)) + cumN;    % <- Offset by cumN
        DUR = double(H.EVENT.DUR(:));
        if isfield(H.EVENT,'CHN'), CHN = double(H.EVENT.CHN(:)); else, CHN = zeros(size(TYP)); end

        events.TYP = [events.TYP; TYP];
        events.POS = [events.POS; POS];
        events.DUR = [events.DUR; DUR];
        events.CHN = [events.CHN; CHN];
    end

    cumN = cumN + nS;
end

fprintf('Concatenated %d files, total samples: %d\n', numel(flist), size(EEG_all,1));

HDRc = struct('EVENT', events);
Sc = make_labels_from_events(HDRc, size(EEG_all,1), Fs_ref, codeMap);

%% Trial Extraction
% Parameters (edit as per your task/protocol)
winSec     = [0, 4];           % from cue onset to +4 s
Ns_win     = round(diff(winSec)*Fs_ref);
preOffset  = round(winSec(1)*Fs_ref);
postOffset = round(winSec(2)*Fs_ref);

Ak = Sc.Ak;
cueOn = find(diff([0; Ak~=0])>0);  % sample indices of cue onsets

% Only keep cues that can accommodate the full window
cueOn = cueOn(cueOn+postOffset <= size(EEG_all,1));

nTrials = numel(cueOn);
trials  = zeros(postOffset-preOffset, size(EEG_all,2), nTrials, 'double');
Ck      = zeros(nTrials,1,'int16');

for k = 1:nTrials
    a = cueOn(k) + preOffset;
    b = cueOn(k) + postOffset - 1;
    trials(:,:,k) = EEG_all(a:b, :);

    % The cue code at onset:
    Ck(k) = Ak(cueOn(k));
end

fprintf('Built trials: %d (window %.2f–%.2f s)\n', nTrials, winSec(1), winSec(2));







