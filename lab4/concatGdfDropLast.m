%% 3) Load & concatenate all GDFs (drop last channel), merge events
function [allS, fs, labels, EVENT, filesUsed] = concatGdfDropLast(rootFolder)
    files = getAllGdfFiles(rootFolder);

    if isempty(files)
        error('No .gdf files found under: %s', rootFolder);
    end

    allS = [];
    EVENT = struct('TYP',[],'POS',[],'DUR',[]);
    labels = [];
    fs = [];
    cumN = 0;           % cumulative samples for event offsets
    filesUsed = {};

    for f = 1:numel(files)
        thisFile = files{f};
        try
            [s, h] = sload(thisFile);  % requires BioSig on path
        catch ME
            warning('Skipping %s (sload failed): %s', thisFile, ME.message);
            continue;
        end

        % --- drop last channel (trigger/status) ---
        if size(s,2) < 2
            warning('File %s has <2 channels; cannot drop last channel. Skipping.', thisFile);
            continue;
        end
        s = s(:, 1:end-1);
        lbl = string(h.Label);
        if numel(lbl) >= 2
            lbl = lbl(1:end-1);
        else
            lbl = strings(1, size(s,2)); % fallback if labels missing
        end

        % --- sampling rate & label consistency checks ---
        if isempty(fs)
            fs = h.SampleRate;
            labels = lbl;
        else
            if h.SampleRate ~= fs
                warning('Fs mismatch in %s (%.1f vs %.1f). Skipping file.', thisFile, h.SampleRate, fs);
                continue;
            end
            if numel(lbl) ~= numel(labels)
                warning('Channel count mismatch in %s. Skipping file.', thisFile);
                continue;
            end
            % If you want strict label match, enforce below:
            % if ~isequal(lower(labels), lower(lbl))
            %     warning('Channel labels mismatch in %s. Skipping file.', thisFile);
            %     continue;
            % end
        end

        % --- append signal ---
        allS = [allS; s]; %#ok<AGROW>

        % --- append events with position offset ---
        if isfield(h, 'EVENT') && ~isempty(h.EVENT) && ~isempty(h.EVENT.TYP)
            EVENT.TYP = [EVENT.TYP; h.EVENT.TYP(:)];
            EVENT.POS = [EVENT.POS; h.EVENT.POS(:) + cumN];
            if isfield(h.EVENT,'DUR') && ~isempty(h.EVENT.DUR)
                EVENT.DUR = [EVENT.DUR; h.EVENT.DUR(:)];
            else
                EVENT.DUR = [EVENT.DUR; zeros(numel(h.EVENT.TYP),1)];
            end
        end

        % --- update counters, record file used ---
        cumN = size(allS,1);
        filesUsed{end+1,1} = thisFile; %#ok<AGROW>
    end

    if isempty(allS)
        error('No valid signals concatenated. Check folders/files.');
    end
end

