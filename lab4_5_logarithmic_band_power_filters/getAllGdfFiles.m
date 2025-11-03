%% 2) Find all .gdf files recursively under rootFolder
function fileList = getAllGdfFiles(rootFolder)
    d = dir(fullfile(rootFolder, '**', '*.gdf'));
    % Sort for reproducibility (folderpath + filename)
    [~, order] = sort(lower(fullfile({d.folder}, {d.name})));
    d = d(order);
    % Build full paths
    fileList = fullfile({d.folder}, {d.name});
    fileList = fileList(:);
end
