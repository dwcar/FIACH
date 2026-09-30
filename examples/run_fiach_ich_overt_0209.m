% Run the packaged overt_0209 data with the optional ModifiedICH tSNR GMM.
% Run run_overt_0209_validation first when a paired R-compatible MATLAB
% output is also wanted for the GMM comparison viewer.

repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
setup_fiach; % Or: setup_fiach('SPMRoot', 'C:\path\to\spm12')
sourceFolder = fullfile(repositoryRoot, 'test_data', 'overt_0209', 'input');
inputFolder = fullfile(repositoryRoot, 'work', 'overt_0209_ich');
maskFile = fullfile(repositoryRoot, 'test_data', 'overt_0209', ...
    'reference_matlab_diagnostics', 'mask.nii');
if exist(inputFolder, 'dir') ~= 7
    mkdir(inputFolder);
end

listing = dir(fullfile(sourceFolder, 'rfmri_overt_*.nii'));
[~, order] = sort({listing.name});
listing = listing(order);
assert(numel(listing) == 206, 'Expected 206 packaged input volumes.');
for index = 1:numel(listing)
    destination = fullfile(inputFolder, listing(index).name);
    if exist(destination, 'file') ~= 2
        copyfile(fullfile(listing(index).folder, listing(index).name), destination);
    end
end
rpFile = fullfile(inputFolder, 'rp_fmri_overt_0004.txt');
copyfile(fullfile(sourceFolder, 'rp_fmri_overt_0004.txt'), rpFile);
functional = arrayfun(@(item) fullfile(inputFolder, item.name), listing, ...
    'UniformOutput', false);

B0 = 3;
TEms = 26;
TR = 1.25;
t = fiach_bold_contrast(B0, TEms);
result = fiach(functional, t, TR, ...
    'RP', rpFile, 'MaxGap', 1, 'Freq', 128, 'NMads', 1.96, ...
    'UseUserMask', true, 'Mask', maskFile, ...
    'GMMMethod', 'ModifiedICH', 'GMMPosteriorThreshold', 0.5, ...
    'Overwrite', true, 'MakePlots', true);
save(fullfile(inputFolder, 'fiach_ich_result.mat'), ...
    'result', 'B0', 'TEms', 'TR', 't', 'maskFile', 'rpFile');
disp(result)

rCompatibleFolder = fullfile(repositoryRoot, 'work', 'overt_0209_r_validation');
if exist(fullfile(rCompatibleFolder, 'noise_basis6.txt'), 'file') == 2
    comparison = fiach_compare_gmm_options('OpenViewer', usejava('desktop'));
end
