function Sx = applySpatialFilter(S, mode, lapMaskFile)
switch lower(mode)
    case 'none'
        Sx = S;
    case 'car'
        Sx = S - mean(S,2);
    case 'lap'
        L = load(lapMaskFile);
        fn = fieldnames(L);
        mask = L.(fn{1}); % assume only one variable in file
        Sx = S * mask;
    otherwise
        error('Unknown spatial mode.');
end
end
