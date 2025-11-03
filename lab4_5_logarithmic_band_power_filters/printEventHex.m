function printEventHex(EVENT)
    if isempty(EVENT) || ~isfield(EVENT,'TYP') || isempty(EVENT.TYP)
        disp('No events found.');
        return;
    end
    u = unique(EVENT.TYP(:))';
    fprintf('Unique EVENT.TYP codes: ');
    for k = 1:numel(u)
        fprintf('0x%04X ', u(k));
    end
    fprintf('\nTotal events: %d\n', numel(EVENT.TYP));
end


