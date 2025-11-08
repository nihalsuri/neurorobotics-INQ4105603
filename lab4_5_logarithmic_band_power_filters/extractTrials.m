function trials = extractTrials(EVENT, fs, EV_FIX, EV_FB, EV_HND, EV_FT, EV_RST, N, winSec)
% Build trials from fixation to (feedback+dur) when present; otherwise a fixed window.
% Inputs:
%   EVENT: struct with fields TYP, POS, DUR (samples)
%   fs:    sampling rate
%   EV_*:  event codes (hex2dec)
%   N:     total number of samples in the EEG (size(allS,1))
%   winSec: fallback trial length in seconds (e.g., 8)

    if nargin < 9 || isempty(winSec), winSec = 8; end
    winSamp = round(winSec * fs);

    TYP = EVENT.TYP(:);
    POS = EVENT.POS(:);
    DUR = EVENT.DUR(:);
    if isempty(DUR), DUR = zeros(size(POS)); end

    % make sure time-ordered
    if ~issorted(POS)
        [POS,ord] = sort(POS);
        TYP = TYP(ord); DUR = DUR(ord);
    end

    % event times we care about
    fixPos  = POS(TYP == EV_FIX);
    fbPos   = POS(TYP == EV_FB);
    fbDur   = DUR(TYP == EV_FB);

    trials = struct('start',{},'stop',{},'label',{});

    for i = 1:numel(fixPos)
        t_start = fixPos(i);

        % primary: first feedback after fixation
        j = find(fbPos > t_start, 1, 'first');
        if ~isempty(j) && fbDur(j) > 0
            t_stop = fbPos(j) + fbDur(j) - 1;   % inclusive
        else
            % fallback: fixed window
            t_stop = t_start + winSamp - 1;
        end

        % cap to signal length (THIS is the key fix)
        t_stop = min(t_stop, N);
        if t_stop <= t_start, continue; end

        % label by the cue that occurs inside this window
        inwin = (POS >= t_start) & (POS <= t_stop);
        types = TYP(inwin);
        if any(types == EV_HND)
            L = 1;
        elseif any(types == EV_FT)
            L = 2;
        elseif any(types == EV_RST)
            L = 0;
        else
            continue;   % no recognizable cue
        end

        trials(end+1) = struct('start', t_start, 'stop', t_stop, 'label', L); %#ok<AGROW>
    end
end

