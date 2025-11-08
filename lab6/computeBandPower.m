function [muPow, betaPow] = computeBandPower(S, fs, bu_mu, au_mu, bu_beta, au_beta)
avgN = round(1.0 * fs);
avgKernel = ones(avgN,1)/avgN;
muPow = zeros(size(S)); betaPow = zeros(size(S));

for ch = 1:size(S,2)
    x_mu   = filtfilt(bu_mu, au_mu, S(:,ch));
    x_beta = filtfilt(bu_beta, au_beta, S(:,ch));
    muPow(:,ch)   = log10(filter(avgKernel,1,x_mu.^2) + eps);
    betaPow(:,ch) = log10(filter(avgKernel,1,x_beta.^2) + eps);
end
end
