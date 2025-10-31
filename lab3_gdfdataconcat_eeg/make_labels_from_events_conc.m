function S = make_labels_from_events_conc(E, nSamples, codeMap)
% Build per-sample label vectors from a concatenated event table E
% Fields in E: TYP, POS, DUR, CHN (optional; treated as global if missing)

S.Tk  = zeros(nSamples,1,'int32');  % trial index (set below)
S.Fk  = zeros(nSamples,1,'int16');  % fixation
S.Ak  = zeros(nSamples,1,'int16');  % cue
S.CFk = zeros(nSamples,1,'int16');  % feedback
S.Xk  = zeros(nSamples,1,'int16');  % outcome (hit/miss)

TYP = double(E.TYP(:));
POS = double(E.POS(:));
DUR = double(E.DUR(:));
if isfield(E,'CHN') && ~isempty(E.CHN)
    CHN = double(E.CHN(:)); 
end

% Assign labels period-by-period
periods = fieldnames(codeMap);
for p = 1:numel(periods)
    fld   = periods{p};
    codes = codeMap.(fld);
    for c = codes(:)'
        k = find(TYP == c);
        for j = 1:numel(k)
            a = POS(k(j));
            b = a + max(DUR(k(j)) - 1, 0);
            a = max(a,1); b = min(b,nSamples);
            if a <= b
                switch upper(fld)
                    case 'FIX',  S.Fk(a:b)  = c;
                    case 'CUE',  S.Ak(a:b)  = c;
                    case 'FEED', S.CFk(a:b) = c;
                    case {'HIT','MISS'}
                        S.Xk(a:b) = c;
                end
            end
        end
    end
end

% Trial counter Tk: each cue onset starts a new trial (until next cue)
cueMask = S.Ak ~= 0;
cueOn   = find(diff([0; cueMask]) > 0);
cueOn(end+1) = nSamples + 1;   % sentinel
for i = 1:numel(cueOn)-1
    S.Tk(cueOn(i):cueOn(i+1)-1) = i;
end
end
