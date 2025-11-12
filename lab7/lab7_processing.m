%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAB 07 - ERD/ERS on spectogram
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

files = getAllGdfFiles(gdfFolder);
files = files(1:7); 
% % %% ----------------COMPUTE PSD--------------------
wlength = 0.5; % seconds. Length of the external window
pshift = 0.25; % seconds. Shift of the internal windows
wshift = 0.0625; % seconds. Shift of the external window -red window - time to process everything
mlength = 1; % seconds

outputDir = fullfile(pwd, 'data'); % local 'data' folder in current working directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% -------- PROCESS ALL GDFs: save ONLY PSD, f, POS, DUR --------
for iFile = 1:numel(files)
    gdfFile = files{iFile};
    fprintf('\n=== Processing: %s ===\n', gdfFile);

    % Load
    [S, H] = sload(gdfFile);
    % Drop the last col
    S(:,end) = []; 
    fs = H.SampleRate;

    % Spatial filter (uses your earlier settings)
    Sx = applySpatialFilter(S, spatialMode, lapMaskFile);

    % Spectrogram / PSD
    % PSD dims: [nWin x nFreq x nChan], f in Hz
    [PSD, f] = proc_spectrogram(Sx, wlength, wshift, pshift, fs, mlength);

    % Save selected freqs
    f   = f(3:25);
    PSD = PSD(:, 3:25, :);

    % --- Align events to window indices (POS) and durations to windows (DUR) ---
    winconv = 'backward';
    POS = proc_pos2win(H.EVENT.POS, wshift*fs, winconv, wlength*fs);
    % Convert event duration (samples) to number of spectrogram windows
    DUR = H.EVENT.DUR / (wshift*fs);
    TYP = H.EVENT.TYP;

    % --- Save inside local "data" folder with same filename BUT .mat extension
    [~, base, ext] = fileparts(gdfFile);    % base='subject01', ext='.gdf'
    outMat = fullfile(outputDir, [base '.mat']);

    % extra safety in case something weird sneaks in:
    if endsWith(outMat, '.gdf', 'IgnoreCase', true)
     outMat = replace(outMat, '.gdf', '.mat', 'IgnoreCase', true);
    end

    fprintf('Saving to: %s\n', outMat);
    save(outMat, 'PSD', 'f', 'POS', 'DUR', 'TYP', '-v7.3');
    
end

fprintf('\nDone. Saved PSD, f, POS, DUR, TYP for all files.\n');







% % load file by file
% [S, H] = sload(files{1,1});
% % Drop the last col
% S(:,end) = [];
% 
% fprintf('\nSampling rate = %.1f Hz, %d channels.\n', ...
%          H.SampleRate, size(S,2));
% disp('Channel labels:'); disp(H.Label');
% 
% % Define event codes (standard BCI2000 conventions)
% EV_FIX = hex2dec('0312');  % fixation ON
% EV_FB  = hex2dec('030D');  % feedback ON
% EV_HND = hex2dec('0305');  % both hands
% EV_FT  = hex2dec('0303');  % both feet
% 
% % %% -------------------- SPATIAL FILTER --------------------
% Sx = applySpatialFilter(S, spatialMode, lapMaskFile);
% fprintf('Applied spatial filter: %s\n', upper(spatialMode));
% 
% % %% ----------------COMPUTE PSD--------------------
% wlength = 0.5; % seconds. Length of the external window
% pshift = 0.25; % seconds. Shift of the internal windows
% wshift = 0.0625; % seconds. Shift of the external window -red window - time to process everything
% mlength = 1; % seconds
% % PSD : [windows X freq X channels]
% 
% [PSD, f] = proc_spectrogram(Sx, wlength, wshift, pshift, H.SampleRate, mlength);
% 
% % 23 frequencies from 4Hz to 48Hz, with a step of 2 -> informative features
% % 368 (16 X 23)
% 
% % % Select the same subset for the PSD 
% f = f(3:25); 
% PSD = PSD(:, 3:25, :);
% winconv = 'backward';
% % wshift*H.SampleRate = 32 Hz
% POS = proc_pos2win(H.EVENT.POS, wshift*H.SampleRate, winconv, wlength*H.SampleRate);
% DUR = H.EVENT.DUR / wshift; 
























