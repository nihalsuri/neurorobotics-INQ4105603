function S = make_labels_from_events(HDR, nSamples, Fs, codeMap)
% codeMap is a struct with fields: FIX, CUE, FEED, HIT, MISS
% Each field contains an array of event codes (numeric) for that period.

TYP = double(HDR.EVENT.TYP(:));
POS = double(HDR.EVENT.POS(:));   % 1-based
DUR = double(HDR.EVENT.DUR(:));   % samples

S.Tk  = zeros(nSamples,1,'int32');   % trial counter (0 if not in a trial)
S.Fk  = zeros(nSamples,1,'int16');   % fixation (0 or event code)
S.Ak  = zeros(nSamples,1,'int16');   % cue      (0 or event code)
S.CFk = zeros(nSamples,1,'int16');   % feedback (0 or event code)
S.Xk  = zeros(nSamples,1,'int16');   % hit/miss (0 or code)

% 1) Assign periods by event code sets
periods = fieldnames(codeMap);
for p = 1:numel(periods)
    fld = periods{p};
    codes = codeMap.(fld);
    for c = codes(:)'
        idx = find(TYP==c);
        for j = 1:numel(idx)
            a = POS(idx(j));
            b = a + max(DUR(idx(j))-1,0);
            a = max(a,1); b = min(b,nSamples);
            if a<=b
                switch fld
                    case 'FIX',  S.Fk(a:b)  = c;
                    case 'CUE',  S.Ak(a:b)  = c;
                    case 'FEED', S.CFk(a:b) = c;
                    case {'HIT','MISS'}
                        S.Xk(a:b) = c;  % keep code; you can remap later to +1/-1
                end
            end
        end
    end
end

% 2) Build a simple trial counter Tk:
%    Let each CUE event start a new trial and continue it until the next CUE.
cueMask = S.Ak ~= 0;
cueOn   = find(diff([0; cueMask])>0);   % rising edges
cueOn(end+1) = nSamples+1;              % sentinel
for k = 1:numel(cueOn)-1
    S.Tk(cueOn(k):cueOn(k+1)-1) = k;
end
end

