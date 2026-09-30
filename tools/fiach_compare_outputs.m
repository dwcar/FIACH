function results = fiach_compare_outputs(varargin)
%FIACH_COMPARE_OUTPUTS Compare MATLAB and R FIACH result series.
%
%   RESULTS = FIACH_COMPARE_OUTPUTS() compares the 206 matched volumes in
%   the FIACH_test folders used for the overt_0209 validation. It writes a
%   per-volume CSV and summary figure, then opens an interactive viewer.
%
%   The viewer shows MATLAB, R, and MATLAB-minus-R images side by side. Its
%   frame and slice sliders, play button, crosshair, linked voxel timecourse,
%   and per-volume error plot all use the same state. By default, images are
%   masked so the comparison measures FIACH-processed voxels rather than the
%   different out-of-mask output conventions.
%
%   Name-value options:
%     'MatlabFolder' folder containing MATLAB filt_*.nii outputs
%     'RFolder'      folder containing R filt_*.nii outputs
%     'Mask'         aligned 3-D comparison mask; default is MATLAB mask.nii
%     'OutputFolder' destination for CSV, summary text, and PNG
%     'OpenViewer'   open the interactive figure (true)
%     'FrameRate'    playback frames per second (4)
%     'LeftLabel'    label for MatlabFolder data ('MATLAB')
%     'RightLabel'   label for RFolder data ('R')
%
%   SPM12 must be on the MATLAB path.

    thisFolder = fileparts(mfilename('fullpath'));
    repositoryRoot = fileparts(thisFolder);
    defaultMatlabFolder = fullfile(repositoryRoot, 'work', 'overt_0209_r_validation');
    defaultRFolder = fullfile(repositoryRoot, 'test_data', 'overt_0209', 'reference_R');

    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser, 'MatlabFolder', defaultMatlabFolder, @fiach_compare_is_text_scalar);
    addParameter(parser, 'RFolder', defaultRFolder, @fiach_compare_is_text_scalar);
    addParameter(parser, 'Mask', '', @fiach_compare_is_text_scalar);
    addParameter(parser, 'OutputFolder', fullfile(repositoryRoot, 'work', ...
        'fiach_comparison_overt_0209'), ...
        @fiach_compare_is_text_scalar);
    addParameter(parser, 'OpenViewer', true, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'FrameRate', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'LeftLabel', 'MATLAB', @fiach_compare_is_text_scalar);
    addParameter(parser, 'RightLabel', 'R', @fiach_compare_is_text_scalar);
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.MatlabFolder = char(opts.MatlabFolder);
    opts.RFolder = char(opts.RFolder);
    opts.Mask = char(opts.Mask);
    opts.OutputFolder = char(opts.OutputFolder);
    opts.LeftLabel = char(opts.LeftLabel);
    opts.RightLabel = char(opts.RightLabel);

    fiach_compare_require_spm();
    if isempty(opts.Mask)
        opts.Mask = fullfile([opts.MatlabFolder '_fiach_diagnostics'], 'mask.nii');
    end

    [matlabFiles, rFiles, frameNumbers] = fiach_compare_match_files( ...
        opts.MatlabFolder, opts.RFolder);
    Vmatlab = spm_vol(char(matlabFiles));
    Vr = spm_vol(char(rFiles));
    if numel(Vmatlab) ~= numel(matlabFiles) || numel(Vr) ~= numel(rFiles)
        error('fiach_compare_outputs:Unexpected4DFile', ...
            'Each matched file must contain exactly one 3-D volume.');
    end
    if ~isequal(Vmatlab(1).dim(1:3), Vr(1).dim(1:3))
        error('fiach_compare_outputs:DimensionMismatch', ...
            'MATLAB and R image dimensions differ.');
    end

    [mask, Vmask] = fiach_compare_read_mask(opts.Mask, Vmatlab(1));
    affineMaximumDifference = max(abs(Vmatlab(1).mat(:) - Vr(1).mat(:)));
    affineTranslationDifference = norm(Vmatlab(1).mat(1:3, 4) - Vr(1).mat(1:3, 4));
    maskAffineDifference = max(abs(Vmatlab(1).mat(:) - Vmask.mat(:)));
    if maskAffineDifference > 1e-4
        error('fiach_compare_outputs:MaskAffineMismatch', ...
            'The comparison mask affine does not match the MATLAB result grid.');
    end

    nFrames = numel(frameNumbers);
    values = nan(nFrames, 10);
    maskSampleCount = 0;
    maskSquaredDifference = 0;
    maskAbsoluteDifference = 0;
    wholeSampleCount = 0;
    wholeSquaredDifference = 0;
    wholeAbsoluteDifference = 0;
    maskDifferenceAboveOne = 0;
    maskDifferenceAboveFive = 0;
    maskDifferenceAboveHundred = 0;
    maximumMaskDifference = -inf;
    maximumMaskFrame = NaN;
    maximumMaskVoxel = [NaN NaN NaN];
    outsideSampleCount = 0;
    rOutsideZeroCount = 0;
    matlabOutsideEqualsInputCount = 0;
    rawInputsAvailable = true;

    outsideMask = ~mask;
    for index = 1:nFrames
        matlabData = double(spm_read_vols(Vmatlab(index)));
        rData = double(spm_read_vols(Vr(index)));
        valid = isfinite(matlabData) & isfinite(rData);

        wholeDifference = matlabData(valid) - rData(valid);
        maskValid = valid & mask;
        maskedDifference = matlabData(maskValid) - rData(maskValid);
        [wholeMean, wholeMae, wholeRmse, wholeCorrelation] = ...
            fiach_compare_metrics(matlabData(valid), rData(valid));
        [maskMean, maskMae, maskRmse, maskCorrelation] = ...
            fiach_compare_metrics(matlabData(maskValid), rData(maskValid));
        absoluteMaskedDifference = abs(maskedDifference);
        [frameMaximum, frameMaximumIndex] = max(absoluteMaskedDifference);
        values(index, :) = [frameNumbers(index), wholeMean, wholeMae, wholeRmse, ...
            wholeCorrelation, maskMean, maskMae, maskRmse, maskCorrelation, frameMaximum];

        wholeSampleCount = wholeSampleCount + numel(wholeDifference);
        wholeSquaredDifference = wholeSquaredDifference + sum(wholeDifference .^ 2);
        wholeAbsoluteDifference = wholeAbsoluteDifference + sum(abs(wholeDifference));
        maskSampleCount = maskSampleCount + numel(maskedDifference);
        maskSquaredDifference = maskSquaredDifference + sum(maskedDifference .^ 2);
        maskAbsoluteDifference = maskAbsoluteDifference + sum(absoluteMaskedDifference);
        maskDifferenceAboveOne = maskDifferenceAboveOne + nnz(absoluteMaskedDifference > 1);
        maskDifferenceAboveFive = maskDifferenceAboveFive + nnz(absoluteMaskedDifference > 5);
        maskDifferenceAboveHundred = maskDifferenceAboveHundred + nnz(absoluteMaskedDifference > 100);
        if frameMaximum > maximumMaskDifference
            maximumMaskDifference = frameMaximum;
            maximumMaskFrame = frameNumbers(index);
            validMaskIndices = find(maskValid);
            maximumLinearIndex = validMaskIndices(frameMaximumIndex);
            [maximumMaskVoxel(1), maximumMaskVoxel(2), maximumMaskVoxel(3)] = ...
                ind2sub(size(mask), maximumLinearIndex);
        end

        outsideValid = valid & outsideMask;
        outsideSampleCount = outsideSampleCount + nnz(outsideValid);
        rOutsideZeroCount = rOutsideZeroCount + nnz(rData(outsideValid) == 0);

        rawFile = fiach_compare_raw_input_path(matlabFiles{index});
        if exist(rawFile, 'file') == 2
            Vraw = spm_vol(rawFile);
            rawData = double(spm_read_vols(Vraw));
            rawValid = outsideValid & isfinite(rawData);
            matlabOutsideEqualsInputCount = matlabOutsideEqualsInputCount + ...
                nnz(matlabData(rawValid) == rawData(rawValid));
        else
            rawInputsAvailable = false;
        end
    end

    metrics = array2table(values, 'VariableNames', { ...
        'Frame', 'WholeMeanDifference', 'WholeMAE', 'WholeRMSE', 'WholeCorrelation', ...
        'MaskMeanDifference', 'MaskMAE', 'MaskRMSE', 'MaskCorrelation', 'MaskMaximumAbsoluteDifference'});

    summary = struct();
    summary.nFrames = nFrames;
    summary.dimensions = Vmatlab(1).dim(1:3);
    summary.maskVoxelCount = nnz(mask);
    summary.maskFraction = nnz(mask) / numel(mask);
    summary.wholeRMSE = sqrt(wholeSquaredDifference / wholeSampleCount);
    summary.wholeMAE = wholeAbsoluteDifference / wholeSampleCount;
    summary.maskRMSE = sqrt(maskSquaredDifference / maskSampleCount);
    summary.maskMAE = maskAbsoluteDifference / maskSampleCount;
    summary.meanMaskCorrelation = mean(metrics.MaskCorrelation);
    summary.maskFractionAboveOne = maskDifferenceAboveOne / maskSampleCount;
    summary.maskFractionAboveFive = maskDifferenceAboveFive / maskSampleCount;
    summary.maskFractionAboveHundred = maskDifferenceAboveHundred / maskSampleCount;
    summary.maximumMaskDifference = maximumMaskDifference;
    summary.maximumMaskDifferenceFrame = maximumMaskFrame;
    summary.maximumMaskDifferenceVoxel = maximumMaskVoxel;
    summary.rOutsideMaskZeroFraction = rOutsideZeroCount / outsideSampleCount;
    if rawInputsAvailable
        summary.matlabOutsideMaskEqualsInputFraction = ...
            matlabOutsideEqualsInputCount / outsideSampleCount;
    else
        summary.matlabOutsideMaskEqualsInputFraction = NaN;
    end
    summary.affineMaximumElementDifference = affineMaximumDifference;
    summary.affineTranslationDifferenceMm = affineTranslationDifference;

    if exist(opts.OutputFolder, 'dir') ~= 7
        mkdir(opts.OutputFolder);
    end
    metricsFile = fullfile(opts.OutputFolder, 'fiach_volume_comparison.csv');
    summaryFile = fullfile(opts.OutputFolder, 'fiach_comparison_summary.txt');
    figureFile = fullfile(opts.OutputFolder, 'fiach_volume_comparison.png');
    writetable(metrics, metricsFile);
    fiach_compare_write_summary(summaryFile, summary, opts, matlabFiles, rFiles);
    fiach_compare_write_metric_figure(figureFile, metrics, opts.LeftLabel, opts.RightLabel);

    viewerFigure = [];
    if opts.OpenViewer
        viewerFigure = fiach_compare_viewer(Vmatlab, Vr, frameNumbers, mask, metrics, ...
            summary, opts.FrameRate, opts.LeftLabel, opts.RightLabel);
    end

    results = struct();
    results.matlabFiles = matlabFiles;
    results.rFiles = rFiles;
    results.maskFile = opts.Mask;
    results.metrics = metrics;
    results.summary = summary;
    results.metricsFile = metricsFile;
    results.summaryFile = summaryFile;
    results.figureFile = figureFile;
    results.viewerFigure = viewerFigure;
    results.leftLabel = opts.LeftLabel;
    results.rightLabel = opts.RightLabel;

    fprintf(['Compared %d matched FIACH volumes. In-mask RMSE %.6g; ' ...
        'whole-grid RMSE %.6g.\n'], nFrames, summary.maskRMSE, summary.wholeRMSE);
end


function tf = fiach_compare_is_text_scalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end


function fiach_compare_require_spm()
    required = {'spm_vol', 'spm_read_vols', 'spm_get_data'};
    for index = 1:numel(required)
        if exist(required{index}, 'file') ~= 2
            error('fiach_compare_outputs:SPMNotFound', ...
                'SPM12 must be on the MATLAB path. Missing: %s', required{index});
        end
    end
end


function [matlabFiles, rFiles, frames] = fiach_compare_match_files(matlabFolder, rFolder)
    if exist(matlabFolder, 'dir') ~= 7
        error('fiach_compare_outputs:MatlabFolderMissing', ...
            'MATLAB output folder not found: %s', matlabFolder);
    end
    if exist(rFolder, 'dir') ~= 7
        error('fiach_compare_outputs:RFolderMissing', ...
            'R output folder not found: %s', rFolder);
    end
    matlabListing = dir(fullfile(matlabFolder, 'filt_*.nii'));
    rListing = dir(fullfile(rFolder, 'filt_*.nii'));
    matlabNames = {matlabListing.name};
    rNames = {rListing.name};
    commonNames = intersect(matlabNames, rNames, 'stable');
    if isempty(commonNames)
        error('fiach_compare_outputs:NoMatchedFiles', ...
            'No matching filt_*.nii files were found.');
    end
    frames = nan(numel(commonNames), 1);
    keep = false(numel(commonNames), 1);
    for index = 1:numel(commonNames)
        token = regexp(commonNames{index}, '_(\d+)\.nii$', 'tokens', 'once');
        if ~isempty(token)
            frames(index) = str2double(token{1});
            keep(index) = true;
        end
    end
    commonNames = commonNames(keep);
    frames = frames(keep);
    [frames, order] = sort(frames);
    commonNames = commonNames(order);
    matlabFiles = cellfun(@(name) fullfile(matlabFolder, name), commonNames, ...
        'UniformOutput', false);
    rFiles = cellfun(@(name) fullfile(rFolder, name), commonNames, ...
        'UniformOutput', false);
end


function [mask, Vmask] = fiach_compare_read_mask(maskFile, Vreference)
    if exist(maskFile, 'file') ~= 2
        error('fiach_compare_outputs:MaskMissing', ...
            'Comparison mask not found: %s', maskFile);
    end
    Vmask = spm_vol(maskFile);
    if numel(Vmask) ~= 1 || ~isequal(Vmask.dim(1:3), Vreference.dim(1:3))
        error('fiach_compare_outputs:MaskDimensionMismatch', ...
            'Mask must be one 3-D image matching the result dimensions.');
    end
    mask = spm_read_vols(Vmask);
    mask = isfinite(mask) & mask ~= 0;
    if ~any(mask(:))
        error('fiach_compare_outputs:EmptyMask', 'The comparison mask is empty.');
    end
end


function [meanDifference, mae, rmse, correlation] = fiach_compare_metrics(left, right)
    difference = left - right;
    meanDifference = mean(difference);
    mae = mean(abs(difference));
    rmse = sqrt(mean(difference .^ 2));
    leftCentred = left - mean(left);
    rightCentred = right - mean(right);
    denominator = sqrt(sum(leftCentred .^ 2) * sum(rightCentred .^ 2));
    if denominator == 0
        correlation = NaN;
    else
        correlation = sum(leftCentred .* rightCentred) / denominator;
    end
end


function rawFile = fiach_compare_raw_input_path(filteredFile)
    [folder, name, extension] = fileparts(filteredFile);
    if startsWith(name, 'filt_')
        name = name(6:end);
    end
    rawFile = fullfile(folder, [name extension]);
end


function fiach_compare_write_summary(file, summary, opts, matlabFiles, rFiles)
    fid = fopen(file, 'w');
    if fid < 0
        error('fiach_compare_outputs:SummaryWriteFailed', 'Could not write %s.', file);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'FIACH %s versus %s comparison\n', opts.LeftLabel, opts.RightLabel);
    fprintf(fid, '%s folder: %s\n', opts.LeftLabel, opts.MatlabFolder);
    fprintf(fid, '%s folder: %s\n', opts.RightLabel, opts.RFolder);
    fprintf(fid, 'First pair: %s | %s\n', matlabFiles{1}, rFiles{1});
    fprintf(fid, 'Matched volumes: %d\n', summary.nFrames);
    fprintf(fid, 'Dimensions: %d x %d x %d\n', summary.dimensions);
    fprintf(fid, 'Mask voxels: %d (%.6f%%)\n', summary.maskVoxelCount, 100 * summary.maskFraction);
    fprintf(fid, 'Whole-grid RMSE: %.12g\n', summary.wholeRMSE);
    fprintf(fid, 'Whole-grid MAE: %.12g\n', summary.wholeMAE);
    fprintf(fid, 'In-mask RMSE: %.12g\n', summary.maskRMSE);
    fprintf(fid, 'In-mask MAE: %.12g\n', summary.maskMAE);
    fprintf(fid, 'Mean in-mask correlation: %.12g\n', summary.meanMaskCorrelation);
    fprintf(fid, 'In-mask |difference| > 1: %.8f%%\n', 100 * summary.maskFractionAboveOne);
    fprintf(fid, 'In-mask |difference| > 5: %.8f%%\n', 100 * summary.maskFractionAboveFive);
    fprintf(fid, 'In-mask |difference| > 100: %.8f%%\n', 100 * summary.maskFractionAboveHundred);
    fprintf(fid, 'Maximum in-mask |difference|: %.12g at frame %d, voxel [%d %d %d]\n', ...
        summary.maximumMaskDifference, summary.maximumMaskDifferenceFrame, ...
        summary.maximumMaskDifferenceVoxel);
    fprintf(fid, '%s outside-mask zero fraction: %.8f%%\n', ...
        opts.RightLabel, 100 * summary.rOutsideMaskZeroFraction);
    fprintf(fid, '%s outside-mask equals input fraction: %.8f%%\n', ...
        opts.LeftLabel, 100 * summary.matlabOutsideMaskEqualsInputFraction);
    fprintf(fid, 'Affine maximum element difference: %.12g\n', ...
        summary.affineMaximumElementDifference);
    fprintf(fid, 'Affine translation displacement: %.12g mm\n', ...
        summary.affineTranslationDifferenceMm);
    clear cleaner;
end


function fiach_compare_write_metric_figure(file, metrics, leftLabel, rightLabel)
    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1100 700]);
    cleaner = onCleanup(@() close(figureHandle));
    layout = tiledlayout(figureHandle, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    wholeAxes = nexttile(layout);
    plot(wholeAxes, metrics.Frame, metrics.WholeRMSE, 'LineWidth', 1.25);
    grid(wholeAxes, 'on');
    ylabel(wholeAxes, 'Whole-grid RMSE');
    title(wholeAxes, sprintf('%s versus %s FIACH difference by volume', leftLabel, rightLabel));
    maskAxes = nexttile(layout);
    plot(maskAxes, metrics.Frame, metrics.MaskRMSE, 'LineWidth', 1.25);
    hold(maskAxes, 'on');
    plot(maskAxes, metrics.Frame, metrics.MaskMAE, 'LineWidth', 1.25);
    grid(maskAxes, 'on');
    xlabel(maskAxes, 'Source volume number');
    ylabel(maskAxes, 'In-mask error');
    legend(maskAxes, {'RMSE', 'MAE'}, 'Location', 'best');
    exportgraphics(figureHandle, file, 'Resolution', 160);
    clear cleaner;
end


function figureHandle = fiach_compare_viewer(Vmatlab, Vr, frameNumbers, mask, metrics, summary, ...
        frameRate, leftLabel, rightLabel)
    dimensions = size(mask);
    sliceCounts = squeeze(sum(sum(mask, 1), 2));
    [~, initialSlice] = max(sliceCounts);
    [initialX, initialY] = find(mask(:, :, initialSlice));
    if isempty(initialX)
        initialVoxel = [round(dimensions(1) / 2), round(dimensions(2) / 2), initialSlice];
    else
        initialVoxel = [round(median(initialX)), round(median(initialY)), initialSlice];
    end

    figureHandle = figure('Name', sprintf('FIACH %s versus %s comparison', leftLabel, rightLabel), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [50 50 1450 880], ...
        'CloseRequestFcn', @closeViewer);
    matlabAxes = axes(figureHandle, 'Position', [0.035 0.53 0.29 0.40]);
    rAxes = axes(figureHandle, 'Position', [0.355 0.53 0.29 0.40]);
    differenceAxes = axes(figureHandle, 'Position', [0.675 0.53 0.29 0.40]);
    timecourseAxes = axes(figureHandle, 'Position', [0.055 0.19 0.57 0.25]);
    metricAxes = axes(figureHandle, 'Position', [0.68 0.19 0.28 0.25]);

    frameSlider = uicontrol(figureHandle, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.17 0.105 0.45 0.03], 'Min', 1, 'Max', numel(frameNumbers), ...
        'Value', 1, 'SliderStep', [1 / max(1, numel(frameNumbers) - 1), ...
        min(10 / max(1, numel(frameNumbers) - 1), 1)], 'Callback', @frameChanged);
    uicontrol(figureHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.055 0.102 0.105 0.035], 'String', 'Frame', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
    frameText = uicontrol(figureHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.63 0.102 0.09 0.035], 'String', '', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
    sliceSlider = uicontrol(figureHandle, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.17 0.062 0.45 0.03], 'Min', 1, 'Max', dimensions(3), ...
        'Value', initialSlice, 'SliderStep', [1 / max(1, dimensions(3) - 1), ...
        min(5 / max(1, dimensions(3) - 1), 1)], 'Callback', @sliceChanged);
    uicontrol(figureHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.055 0.059 0.105 0.035], 'String', 'Axial slice', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
    sliceText = uicontrol(figureHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.63 0.059 0.09 0.035], 'String', '', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
    uicontrol(figureHandle, 'Style', 'togglebutton', 'Units', 'normalized', ...
        'Position', [0.75 0.095 0.09 0.045], 'String', 'Play', 'Callback', @playChanged);
    uicontrol(figureHandle, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.85 0.095 0.13 0.045], 'String', 'Show outside mask', ...
        'BackgroundColor', 'w', 'Value', 0, 'Callback', @outsideChanged);
    statusText = uicontrol(figureHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.055 0.012 0.91 0.035], 'String', '', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

    plot(metricAxes, metrics.Frame, max(metrics.WholeRMSE, eps), '--', ...
        'LineWidth', 1.0, 'Color', [0.45 0.45 0.45]);
    hold(metricAxes, 'on');
    plot(metricAxes, metrics.Frame, max(metrics.MaskRMSE, eps), '-', ...
        'LineWidth', 1.2, 'Color', [0.1 0.35 0.75]);
    set(metricAxes, 'YScale', 'log');
    grid(metricAxes, 'on');
    xlabel(metricAxes, 'Source volume');
    ylabel(metricAxes, 'RMSE (log scale)');
    title(metricAxes, 'Volume-by-volume difference');
    legend(metricAxes, {'Whole grid', 'Inside FIACH mask'}, 'Location', 'best');
    metricFrameLine = xline(metricAxes, frameNumbers(1), '-', 'Color', [0.85 0.2 0.1]);

    playbackTimer = timer('ExecutionMode', 'fixedSpacing', 'BusyMode', 'drop', ...
        'Period', 1 / frameRate, 'TimerFcn', @advanceFrame);
    state = struct('frameIndex', 1, 'slice', initialSlice, 'voxel', initialVoxel, ...
        'showOutside', false, 'matlabImage', [], 'rImage', [], 'differenceImage', [], ...
        'crosshair', gobjects(3, 2), 'timecourseFrameLine', [], ...
        'cachedVoxel', [NaN NaN NaN], 'matlabTimecourse', [], 'rTimecourse', []);
    guidata(figureHandle, state);
    renderFrame(true);

    function renderFrame(updateTimecourse)
        if ~ishghandle(figureHandle)
            return;
        end
        current = guidata(figureHandle);
        matlabData = double(spm_read_vols(Vmatlab(current.frameIndex)));
        rData = double(spm_read_vols(Vr(current.frameIndex)));
        difference = matlabData - rData;
        displayMatlab = matlabData;
        displayR = rData;
        displayDifference = difference;
        if ~current.showOutside
            displayMatlab(~mask) = 0;
            displayR(~mask) = 0;
            displayDifference(~mask) = 0;
        end

        matlabSlice = displayMatlab(:, :, current.slice).';
        rSlice = displayR(:, :, current.slice).';
        differenceSlice = displayDifference(:, :, current.slice).';
        scaleValues = [matlabData(mask); rData(mask)];
        intensityLimits = fiach_compare_robust_limits(scaleValues, 0.01, 0.99);
        differenceLimit = fiach_compare_robust_abs_limit(difference(mask), 0.99);

        if isempty(current.matlabImage) || ~ishghandle(current.matlabImage)
            current.matlabImage = imagesc(matlabAxes, matlabSlice, intensityLimits);
            current.rImage = imagesc(rAxes, rSlice, intensityLimits);
            current.differenceImage = imagesc(differenceAxes, differenceSlice, ...
                [-differenceLimit differenceLimit]);
            set([current.matlabImage current.rImage current.differenceImage], ...
                'ButtonDownFcn', @imageClicked);
            axis(matlabAxes, 'image'); axis(rAxes, 'image'); axis(differenceAxes, 'image');
            set([matlabAxes rAxes differenceAxes], 'XTick', [], 'YTick', []);
            colormap(matlabAxes, gray(256));
            colormap(rAxes, gray(256));
            colormap(differenceAxes, fiach_compare_diverging_colormap(256));
            colorbar(differenceAxes);
        else
            set(current.matlabImage, 'CData', matlabSlice);
            set(current.rImage, 'CData', rSlice);
            set(current.differenceImage, 'CData', differenceSlice);
            clim(matlabAxes, intensityLimits);
            clim(rAxes, intensityLimits);
            clim(differenceAxes, [-differenceLimit differenceLimit]);
        end
        title(matlabAxes, sprintf('%s — frame %04d', leftLabel, frameNumbers(current.frameIndex)));
        title(rAxes, sprintf('%s — frame %04d', rightLabel, frameNumbers(current.frameIndex)));
        title(differenceAxes, sprintf('%s - %s', leftLabel, rightLabel));
        frameSlider.Value = current.frameIndex;
        sliceSlider.Value = current.slice;
        frameText.String = sprintf('%d / %d', current.frameIndex, numel(frameNumbers));
        sliceText.String = sprintf('%d / %d', current.slice, dimensions(3));
        metricFrameLine.Value = frameNumbers(current.frameIndex);
        guidata(figureHandle, current);
        updateCrosshairs();
        if updateTimecourse
            updateTimecoursePlot();
        elseif ~isempty(current.timecourseFrameLine) && ishghandle(current.timecourseFrameLine)
            current.timecourseFrameLine.Value = frameNumbers(current.frameIndex);
        end
        updateStatus(matlabData, rData);
        drawnow limitrate;
    end

    function updateCrosshairs()
        current = guidata(figureHandle);
        axesList = [matlabAxes rAxes differenceAxes];
        for axesIndex = 1:3
            if ~ishghandle(current.crosshair(axesIndex, 1))
                hold(axesList(axesIndex), 'on');
                current.crosshair(axesIndex, 1) = line(axesList(axesIndex), ...
                    [current.voxel(1) current.voxel(1)], [0.5 dimensions(2) + 0.5], ...
                    'Color', [1 0.75 0], 'LineWidth', 1.0, 'HitTest', 'off');
                current.crosshair(axesIndex, 2) = line(axesList(axesIndex), ...
                    [0.5 dimensions(1) + 0.5], [current.voxel(2) current.voxel(2)], ...
                    'Color', [1 0.75 0], 'LineWidth', 1.0, 'HitTest', 'off');
            else
                set(current.crosshair(axesIndex, 1), 'XData', ...
                    [current.voxel(1) current.voxel(1)]);
                set(current.crosshair(axesIndex, 2), 'YData', ...
                    [current.voxel(2) current.voxel(2)]);
            end
        end
        guidata(figureHandle, current);
    end

    function updateTimecoursePlot()
        current = guidata(figureHandle);
        if ~isequal(current.cachedVoxel, current.voxel)
            coordinates = current.voxel(:);
            current.matlabTimecourse = spm_get_data(Vmatlab, coordinates);
            current.rTimecourse = spm_get_data(Vr, coordinates);
            current.cachedVoxel = current.voxel;
        end
        cla(timecourseAxes);
        plot(timecourseAxes, frameNumbers, current.matlabTimecourse, '-', ...
            'LineWidth', 1.1, 'Color', [0.1 0.35 0.75]);
        hold(timecourseAxes, 'on');
        plot(timecourseAxes, frameNumbers, current.rTimecourse, '--', ...
            'LineWidth', 1.1, 'Color', [0.85 0.25 0.15]);
        current.timecourseFrameLine = xline(timecourseAxes, ...
            frameNumbers(current.frameIndex), '-', 'Color', [0.2 0.2 0.2]);
        grid(timecourseAxes, 'on');
        xlabel(timecourseAxes, 'Source volume');
        ylabel(timecourseAxes, 'Signal intensity');
        title(timecourseAxes, sprintf('Voxel timecourse at [%d %d %d]', current.voxel));
        legend(timecourseAxes, {leftLabel, rightLabel}, 'Location', 'best');
        guidata(figureHandle, current);
    end

    function updateStatus(matlabData, rData)
        current = guidata(figureHandle);
        linearIndex = sub2ind(dimensions, current.voxel(1), current.voxel(2), current.voxel(3));
        insideText = 'outside mask';
        if mask(linearIndex)
            insideText = 'inside mask';
        end
        statusText.String = sprintf([ ...
            'Voxel [%d %d %d], %s: %s %.6g, %s %.6g, difference %.6g. ' ...
            'Array-index view; header translation differs by %.3f mm.'], ...
            current.voxel, insideText, leftLabel, matlabData(linearIndex), ...
            rightLabel, rData(linearIndex), ...
            matlabData(linearIndex) - rData(linearIndex), summary.affineTranslationDifferenceMm);
    end

    function imageClicked(source, ~)
        current = guidata(figureHandle);
        clickedAxes = ancestor(source, 'axes');
        point = clickedAxes.CurrentPoint;
        x = min(max(round(point(1, 1)), 1), dimensions(1));
        y = min(max(round(point(1, 2)), 1), dimensions(2));
        current.voxel = [x y current.slice];
        guidata(figureHandle, current);
        updateCrosshairs();
        updateTimecoursePlot();
        renderFrame(false);
    end

    function frameChanged(source, ~)
        current = guidata(figureHandle);
        current.frameIndex = min(max(round(source.Value), 1), numel(frameNumbers));
        guidata(figureHandle, current);
        renderFrame(false);
    end

    function sliceChanged(source, ~)
        current = guidata(figureHandle);
        current.slice = min(max(round(source.Value), 1), dimensions(3));
        current.voxel(3) = current.slice;
        current.cachedVoxel = [NaN NaN NaN];
        guidata(figureHandle, current);
        renderFrame(true);
    end

    function outsideChanged(source, ~)
        current = guidata(figureHandle);
        current.showOutside = logical(source.Value);
        guidata(figureHandle, current);
        renderFrame(false);
    end

    function playChanged(source, ~)
        if logical(source.Value)
            source.String = 'Pause';
            start(playbackTimer);
        else
            source.String = 'Play';
            if strcmp(playbackTimer.Running, 'on')
                stop(playbackTimer);
            end
        end
    end

    function advanceFrame(~, ~)
        if ~ishghandle(figureHandle)
            return;
        end
        current = guidata(figureHandle);
        current.frameIndex = mod(current.frameIndex, numel(frameNumbers)) + 1;
        guidata(figureHandle, current);
        renderFrame(false);
    end

    function closeViewer(~, ~)
        if strcmp(playbackTimer.Running, 'on')
            stop(playbackTimer);
        end
        delete(playbackTimer);
        delete(figureHandle);
    end
end


function limits = fiach_compare_robust_limits(values, lowerProbability, upperProbability)
    values = sort(values(isfinite(values)));
    if isempty(values)
        limits = [0 1];
        return;
    end
    lowerIndex = max(1, min(numel(values), round(1 + (numel(values) - 1) * lowerProbability)));
    upperIndex = max(1, min(numel(values), round(1 + (numel(values) - 1) * upperProbability)));
    limits = [values(lowerIndex) values(upperIndex)];
    if limits(1) == limits(2)
        limits = limits + [-0.5 0.5];
    end
end


function limit = fiach_compare_robust_abs_limit(values, probability)
    values = sort(abs(values(isfinite(values))));
    if isempty(values)
        limit = 1;
        return;
    end
    index = max(1, min(numel(values), round(1 + (numel(values) - 1) * probability)));
    limit = max(values(index), 1);
end


function map = fiach_compare_diverging_colormap(n)
    half = ceil(n / 2);
    blue = [0.10 0.32 0.75];
    white = [1 1 1];
    red = [0.82 0.16 0.12];
    lower = [linspace(blue(1), white(1), half).', ...
        linspace(blue(2), white(2), half).', linspace(blue(3), white(3), half).'];
    upperCount = n - half;
    upper = [linspace(white(1), red(1), upperCount + 1).', ...
        linspace(white(2), red(2), upperCount + 1).', ...
        linspace(white(3), red(3), upperCount + 1).'];
    map = [lower; upper(2:end, :)];
end
