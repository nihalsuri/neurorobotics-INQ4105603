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
