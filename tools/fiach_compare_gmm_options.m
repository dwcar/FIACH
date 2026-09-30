function result = fiach_compare_gmm_options(varargin)
%FIACH_COMPARE_GMM_OPTIONS Compare R-compatible and ModifiedICH FIACH runs.
%
%   RESULT = FIACH_COMPARE_GMM_OPTIONS() opens the existing movie viewer for
%   the two filtered series, compares the tSNR/noise-mask diagnostics, and
%   creates the separate three-panel PCA-regressor comparison.

    thisFolder = fileparts(mfilename('fullpath'));
    repositoryRoot = fileparts(thisFolder);
    defaultRFolder = fullfile(repositoryRoot, 'work', 'overt_0209_r_validation');
    defaultICHFolder = fullfile(repositoryRoot, 'work', 'overt_0209_ich');
    defaultRDiagnostics = [defaultRFolder '_fiach_diagnostics'];
    defaultICHDiagnostics = [defaultICHFolder '_fiach_diagnostics'];
    defaultOutputFolder = fullfile(repositoryRoot, 'work', ...
        'fiach_gmm_comparison_overt_0209');

    parser = inputParser;
    addParameter(parser, 'RCompatibleFolder', defaultRFolder, @isTextScalar);
    addParameter(parser, 'ModifiedICHFolder', defaultICHFolder, @isTextScalar);
    addParameter(parser, 'RCompatibleDiagnostics', defaultRDiagnostics, @isTextScalar);
    addParameter(parser, 'ModifiedICHDiagnostics', defaultICHDiagnostics, @isTextScalar);
    addParameter(parser, 'OutputFolder', defaultOutputFolder, @isTextScalar);
    addParameter(parser, 'OpenViewer', true, @(value) islogical(value) && isscalar(value));
    parse(parser, varargin{:});
    opts = parser.Results;
    fields = {'RCompatibleFolder', 'ModifiedICHFolder', 'RCompatibleDiagnostics', ...
        'ModifiedICHDiagnostics', 'OutputFolder'};
    for index = 1:numel(fields)
        opts.(fields{index}) = char(opts.(fields{index}));
    end

    if exist('spm_vol', 'file') ~= 2
        error('fiach_compare_gmm_options:SPMNotFound', 'SPM12 must be on the MATLAB path.');
    end
    if exist(opts.OutputFolder, 'dir') ~= 7
        mkdir(opts.OutputFolder);
    end

    rMaskFile = fullfile(opts.RCompatibleDiagnostics, 'mask.nii');
    rTsnrFile = fullfile(opts.RCompatibleDiagnostics, 'rtsnr.nii');
    ichTsnrFile = fullfile(opts.ModifiedICHDiagnostics, 'rtsnr.nii');
    rNoiseFile = fullfile(opts.RCompatibleDiagnostics, 'noise_mask.nii');
    ichNoiseFile = fullfile(opts.ModifiedICHDiagnostics, 'noise_mask.nii');

    brainMask = readLogicalVolume(rMaskFile);
    rTsnr = readNumericVolume(rTsnrFile);
    ichTsnr = readNumericVolume(ichTsnrFile);
    rNoise = readLogicalVolume(rNoiseFile);
    ichNoise = readLogicalVolume(ichNoiseFile);
    assert(isequal(size(rTsnr), size(ichTsnr), size(rNoise), size(ichNoise), size(brainMask)), ...
        'FIACH diagnostic image dimensions do not match.');

    validTsnr = brainMask & isfinite(rTsnr) & isfinite(ichTsnr);
    tsnrDifference = ichTsnr(validTsnr) - rTsnr(validTsnr);
    disagreement = xor(rNoise, ichNoise);
    intersectionCount = nnz(rNoise & ichNoise);
    unionCount = nnz(rNoise | ichNoise);

    diagnostics = struct();
    diagnostics.rNoiseVoxelCount = nnz(rNoise);
    diagnostics.ichNoiseVoxelCount = nnz(ichNoise);
    diagnostics.sharedNoiseVoxelCount = intersectionCount;
    diagnostics.rOnlyNoiseVoxelCount = nnz(rNoise & ~ichNoise);
    diagnostics.ichOnlyNoiseVoxelCount = nnz(ichNoise & ~rNoise);
    diagnostics.disagreementVoxelCount = nnz(disagreement);
    diagnostics.dice = 2 * intersectionCount / max(1, nnz(rNoise) + nnz(ichNoise));
    diagnostics.jaccard = intersectionCount / max(1, unionCount);
    diagnostics.tsnrMaximumAbsoluteDifference = max(abs(tsnrDifference));
    diagnostics.tsnrRmse = sqrt(mean(tsnrDifference .^ 2));

    sliceDisagreement = squeeze(sum(sum(disagreement, 1), 2));
    [~, displaySlice] = max(sliceDisagreement);
    noiseFigureFile = fullfile(opts.OutputFolder, 'gmm_noise_mask_comparison.png');
    noiseFigure = createNoiseFigure(rNoise, ichNoise, displaySlice, diagnostics, ...
        noiseFigureFile, opts.OpenViewer);

    movie = fiach_compare_outputs( ...
        'MatlabFolder', opts.RCompatibleFolder, ...
        'RFolder', opts.ModifiedICHFolder, ...
        'Mask', rMaskFile, ...
        'OutputFolder', opts.OutputFolder, ...
        'LeftLabel', 'R-compatible GMM', ...
        'RightLabel', 'ModifiedICH GMM', ...
        'OpenViewer', opts.OpenViewer);
    pca = fiach_compare_pca_regressors( ...
        'RCompatibleFile', fullfile(opts.RCompatibleFolder, 'noise_basis6.txt'), ...
        'ModifiedICHFile', fullfile(opts.ModifiedICHFolder, 'noise_basis6.txt'), ...
        'OutputFolder', opts.OutputFolder, ...
        'Visible', opts.OpenViewer);

    summaryFile = fullfile(opts.OutputFolder, 'gmm_diagnostic_summary.txt');
    writeSummary(summaryFile, diagnostics, displaySlice);
    result = struct('movie', movie, 'pca', pca, 'diagnostics', diagnostics, ...
        'noiseMaskFigure', noiseFigure, 'noiseMaskFigureFile', noiseFigureFile, ...
        'summaryFile', summaryFile);
    fprintf(['Noise-mask Dice %.6f; %d differing voxels; ' ...
        'tSNR-map RMSE %.6g.\n'], diagnostics.dice, ...
        diagnostics.disagreementVoxelCount, diagnostics.tsnrRmse);
end


function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end


function values = readNumericVolume(file)
    if exist(file, 'file') ~= 2
        error('fiach_compare_gmm_options:FileMissing', 'Diagnostic file not found: %s', file);
    end
    header = spm_vol(file);
    values = double(spm_read_vols(header));
end


function values = readLogicalVolume(file)
    numeric = readNumericVolume(file);
    values = isfinite(numeric) & numeric ~= 0;
end


function figureHandle = createNoiseFigure(rNoise, ichNoise, slice, diagnostics, file, visible)
    visibility = 'off';
    if visible
        visibility = 'on';
    end
    figureHandle = figure('Name', 'FIACH GMM noise-mask comparison', ...
        'NumberTitle', 'off', 'Color', 'w', 'Visible', visibility, ...
        'Position', [100 150 1350 500]);
    layout = tiledlayout(figureHandle, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    rAxes = nexttile(layout);
    imagesc(rAxes, rNoise(:, :, slice).');
    title(rAxes, sprintf('R-compatible noise mask (n = %d)', diagnostics.rNoiseVoxelCount));
    configureImageAxes(rAxes);
    ichAxes = nexttile(layout);
    imagesc(ichAxes, ichNoise(:, :, slice).');
    title(ichAxes, sprintf('ModifiedICH noise mask (n = %d)', diagnostics.ichNoiseVoxelCount));
    configureImageAxes(ichAxes);
    differenceAxes = nexttile(layout);
    difference = double(ichNoise(:, :, slice)) - double(rNoise(:, :, slice));
    imagesc(differenceAxes, difference.', [-1 1]);
    title(differenceAxes, sprintf('Difference: ICH - R (Dice %.4f)', diagnostics.dice));
    configureImageAxes(differenceAxes);
    colormap(rAxes, gray(2));
    colormap(ichAxes, gray(2));
    colormap(differenceAxes, [0.10 0.32 0.75; 1 1 1; 0.82 0.16 0.12]);
    colorbar(differenceAxes, 'Ticks', [-1 0 1], ...
        'TickLabels', {'R only', 'same', 'ICH only'});
    title(layout, sprintf('Expanded tSNR noise classifications, axial slice %d', slice));
    exportgraphics(figureHandle, file, 'Resolution', 180);
end


function configureImageAxes(axesHandle)
    axis(axesHandle, 'image');
    axesHandle.XTick = [];
    axesHandle.YTick = [];
end


function writeSummary(file, diagnostics, displaySlice)
    fid = fopen(file, 'w');
    if fid < 0
        error('fiach_compare_gmm_options:SummaryWriteFailed', 'Could not write %s.', file);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'FIACH R-compatible versus ModifiedICH GMM comparison\n');
    fprintf(fid, 'R-compatible expanded noise voxels: %d\n', diagnostics.rNoiseVoxelCount);
    fprintf(fid, 'ModifiedICH expanded noise voxels: %d\n', diagnostics.ichNoiseVoxelCount);
    fprintf(fid, 'Shared noise voxels: %d\n', diagnostics.sharedNoiseVoxelCount);
    fprintf(fid, 'R-compatible only: %d\n', diagnostics.rOnlyNoiseVoxelCount);
    fprintf(fid, 'ModifiedICH only: %d\n', diagnostics.ichOnlyNoiseVoxelCount);
    fprintf(fid, 'Disagreement voxels: %d\n', diagnostics.disagreementVoxelCount);
    fprintf(fid, 'Dice coefficient: %.12g\n', diagnostics.dice);
    fprintf(fid, 'Jaccard coefficient: %.12g\n', diagnostics.jaccard);
    fprintf(fid, 'tSNR maximum absolute difference: %.12g\n', ...
        diagnostics.tsnrMaximumAbsoluteDifference);
    fprintf(fid, 'tSNR RMSE: %.12g\n', diagnostics.tsnrRmse);
    fprintf(fid, 'Displayed axial slice: %d\n', displaySlice);
    clear cleaner;
end
