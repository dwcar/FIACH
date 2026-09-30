% Reproducible R-compatible FIACH run for the packaged overt_0209 data.

repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
setup_fiach; % Or: setup_fiach('SPMRoot', 'C:\path\to\spm12')
inputFolder = fullfile(repositoryRoot, 'test_data', 'overt_0209', 'input');
workFolder = fullfile(repositoryRoot, 'work', 'overt_0209_example');
if exist(workFolder, 'dir') ~= 7
    mkdir(workFolder);
end

listing = dir(fullfile(inputFolder, 'rfmri_overt_*.nii'));
[~, order] = sort({listing.name});
listing = listing(order);
for index = 1:numel(listing)
    destination = fullfile(workFolder, listing(index).name);
    if exist(destination, 'file') ~= 2
        copyfile(fullfile(listing(index).folder, listing(index).name), destination);
    end
end
copyfile(fullfile(inputFolder, 'rp_fmri_overt_0004.txt'), workFolder);
functional = arrayfun(@(item) fullfile(workFolder, item.name), listing, ...
    'UniformOutput', false);

t = fiach_bold_contrast(3, 26);
result = fiach(functional, t, 1.25, ...
    'RP', fullfile(workFolder, 'rp_fmri_overt_0004.txt'), ...
    'MaxGap', 1, 'Freq', 128, 'NMads', 1.96, ...
    'GMMMethod', 'R', 'Overwrite', true, 'MakePlots', true);
disp(result)
