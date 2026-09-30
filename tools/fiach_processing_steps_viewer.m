function result = fiach_processing_steps_viewer(varargin)
%FIACH_PROCESSING_STEPS_VIEWER Inspect original, FIACH, and nuisance-GLM data.
%
%   RESULT = FIACH_PROCESSING_STEPS_VIEWER() uses the completed overt_0209
%   ModifiedICH run. It projects each in-mask voxel's interpolated FIACH
%   timecourse away from the motion-plus-FIACH nuisance subspace, writes
%   float32 glm_filt_*.nii volumes, and opens a synchronized three-panel
%   movie viewer.
%
%   The nuisance model contains the six realignment parameters and the last
%   six columns of noise_basis6.txt (the FIACH PCA regressors). Regressors
%   are mean-centred before an SVD obtains an orthonormal basis. Therefore:
%
%       cleaned = filtered - Q * (Q' * filtered)
%
%   This is the voxelwise residual-forming projection for the nuisance
%   model while retaining each voxel's temporal mean.
%
%   Name-value options:
%     'RunFolder'       input/output FIACH run directory
%     'Diagnostics'     directory containing mask.nii
%     'MotionFile'      six-column realignment parameter file
%     'RegressorFile'   FIACH noise_basis6.txt
%     'OutputFolder'    destination for GLM volumes and summaries
%     'OverwriteGLM'    replace existing GLM volumes (false)
%     'OpenViewer'      open the synchronized viewer (true)
%     'FrameRate'       movie playback frames per second (4)

    thisFolder = fileparts(mfilename('fullpath'));
    repositoryRoot = fileparts(thisFolder);
    defaultRunFolder = fullfile(repositoryRoot, 'work', 'overt_0209_ich');
    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser, 'RunFolder', defaultRunFolder, @isTextScalar);
    addParameter(parser, 'Diagnostics', '', @isTextScalar);
    addParameter(parser, 'MotionFile', '', @isTextScalar);
    addParameter(parser, 'RegressorFile', '', @isTextScalar);
    addParameter(parser, 'OutputFolder', '', @isTextScalar);
    addParameter(parser, 'OverwriteGLM', false, @(value) islogical(value) && isscalar(value));
    addParameter(parser, 'OpenViewer', true, @(value) islogical(value) && isscalar(value));
    addParameter(parser, 'FrameRate', 4, @(value) isnumeric(value) && isscalar(value) && ...
        isfinite(value) && value > 0);
    parse(parser, varargin{:});
    opts = parser.Results;
    textFields = {'RunFolder', 'Diagnostics', 'MotionFile', 'RegressorFile', 'OutputFolder'};
    for index = 1:numel(textFields)
        opts.(textFields{index}) = char(opts.(textFields{index}));
    end
    if isempty(opts.Diagnostics)
        opts.Diagnostics = [opts.RunFolder '_fiach_diagnostics'];
    end
    if isempty(opts.MotionFile)
        motionListing = dir(fullfile(opts.RunFolder, 'rp_*.txt'));
        if numel(motionListing) ~= 1
            error('fiach_processing_steps_viewer:MotionFileAmbiguous', ...
                'Specify MotionFile; expected one rp_*.txt file in the run folder.');
        end
        opts.MotionFile = fullfile(motionListing.folder, motionListing.name);
    end
    if isempty(opts.RegressorFile)
        opts.RegressorFile = fullfile(opts.RunFolder, 'noise_basis6.txt');
    end
    if isempty(opts.OutputFolder)
        opts.OutputFolder = fullfile(opts.RunFolder, 'glm_cleaned');
    end

    requireSpm();
    [originalFiles, filteredFiles, frameNumbers] = matchSeries(opts.RunFolder);
    nFrames = numel(frameNumbers);
    Voriginal = spm_vol(char(originalFiles));
    Vfiltered = spm_vol(char(filteredFiles));
    validateVolumeSeries(Voriginal, Vfiltered, nFrames);
    [mask, Vmask] = readMask(fullfile(opts.Diagnostics, 'mask.nii'), Vfiltered(1));

    [design, designInfo] = buildNuisanceDesign(opts.MotionFile, opts.RegressorFile, nFrames);
    if exist(opts.OutputFolder, 'dir') ~= 7
        mkdir(opts.OutputFolder);
    end
    glmFiles = makeGlmFilenames(filteredFiles, opts.OutputFolder);
    outputsExist = all(cellfun(@(file) exist(file, 'file') == 2, glmFiles));
    if opts.OverwriteGLM || ~outputsExist
        glmFiles = writeGlmSeries(Vfiltered, glmFiles, mask, design.basis);
    end
    Vglm = spm_vol(char(glmFiles));
    validateVolumeSeries(Vfiltered, Vglm, nFrames);

    designFile = fullfile(opts.OutputFolder, 'glm_nuisance_design.csv');
    designVariables = ["Frame", compose('Motion_%d', 1:designInfo.nMotion), ...
        compose('FIACH_PC%d', 1:designInfo.nFiach)];
    designTable = array2table([frameNumbers(:), design.raw], ...
        'VariableNames', cellstr(designVariables));
    writetable(designTable, designFile);

    [metrics, summary] = calculateMetrics(Voriginal, Vfiltered, Vglm, mask, frameNumbers);
    metricsFile = fullfile(opts.OutputFolder, 'processing_step_metrics.csv');
    summaryFile = fullfile(opts.OutputFolder, 'processing_step_summary.txt');
    figureFile = fullfile(opts.OutputFolder, 'processing_step_volume_effects.png');
    writetable(metrics, metricsFile);
    writeSummary(summaryFile, summary, designInfo, opts);
    writeMetricFigure(figureFile, metrics);

    viewerFigure = [];
    if opts.OpenViewer
        viewerFigure = openViewer(Voriginal, Vfiltered, Vglm, frameNumbers, mask, ...
            metrics, opts.FrameRate);
    end

    result = struct();
    result.originalFiles = originalFiles;
    result.filteredFiles = filteredFiles;
    result.glmFiles = glmFiles;
    result.maskFile = Vmask.fname;
    result.design = design;
    result.designInfo = designInfo;
    result.metrics = metrics;
    result.summary = summary;
    result.designFile = designFile;
    result.metricsFile = metricsFile;
    result.summaryFile = summaryFile;
    result.figureFile = figureFile;
    result.viewerFigure = viewerFigure;

    fprintf(['FIACH processing viewer ready: %d volumes, %d nuisance dimensions; ' ...
        'mean in-mask GLM change RMSE %.6g.\n'], nFrames, designInfo.rank, ...
        mean(metrics.FilteredToGlmRMSE));
end


function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end


function requireSpm()
    required = {'spm_vol', 'spm_read_vols', 'spm_write_vol', 'spm_get_data'};
    for index = 1:numel(required)
        if exist(required{index}, 'file') ~= 2
            error('fiach_processing_steps_viewer:SPMNotFound', ...
                'SPM12 must be on the MATLAB path. Missing: %s', required{index});
        end
    end
end


function [originalFiles, filteredFiles, frames] = matchSeries(folder)
    if exist(folder, 'dir') ~= 7
        error('fiach_processing_steps_viewer:RunFolderMissing', ...
            'FIACH run folder not found: %s', folder);
    end
    listing = dir(fullfile(folder, 'filt_*.nii'));
    if isempty(listing)
        error('fiach_processing_steps_viewer:NoFilteredFiles', ...
            'No filt_*.nii files were found in %s.', folder);
    end
    names = {listing.name};
    frames = nan(numel(names), 1);
    keep = false(numel(names), 1);
    for index = 1:numel(names)
        token = regexp(names{index}, '_(\d+)\.nii$', 'tokens', 'once');
        if ~isempty(token)
            frames(index) = str2double(token{1});
            keep(index) = true;
        end
    end
    names = names(keep);
    frames = frames(keep);
    [frames, order] = sort(frames);
    names = names(order);
    filteredFiles = cellfun(@(name) fullfile(folder, name), names, 'UniformOutput', false);
    originalNames = cellfun(@(name) name(6:end), names, 'UniformOutput', false);
    originalFiles = cellfun(@(name) fullfile(folder, name), originalNames, ...
        'UniformOutput', false);
    missing = originalFiles(cellfun(@(file) exist(file, 'file') ~= 2, originalFiles));
    if ~isempty(missing)
        error('fiach_processing_steps_viewer:OriginalFileMissing', ...
            'Original volume corresponding to a filtered image is missing: %s', missing{1});
    end
end


function validateVolumeSeries(left, right, expectedFrames)
    if numel(left) ~= expectedFrames || numel(right) ~= expectedFrames
        error('fiach_processing_steps_viewer:Unexpected4DFile', ...
            'Each series must contain one 3-D file per frame.');
    end
    if ~isequal(left(1).dim(1:3), right(1).dim(1:3))
        error('fiach_processing_steps_viewer:DimensionMismatch', ...
            'Image series dimensions do not match.');
    end
    if max(abs(left(1).mat(:) - right(1).mat(:))) > 1e-4
        error('fiach_processing_steps_viewer:AffineMismatch', ...
            'Image series affines do not match.');
    end
end


function [mask, Vmask] = readMask(file, reference)
    if exist(file, 'file') ~= 2
        error('fiach_processing_steps_viewer:MaskMissing', 'Mask not found: %s', file);
    end
    Vmask = spm_vol(file);
    if numel(Vmask) ~= 1 || ~isequal(Vmask.dim(1:3), reference.dim(1:3)) || ...
            max(abs(Vmask.mat(:) - reference.mat(:))) > 1e-4
        error('fiach_processing_steps_viewer:MaskMismatch', ...
            'The mask must match the filtered data grid and affine.');
    end
    values = spm_read_vols(Vmask);
    mask = isfinite(values) & values ~= 0;
    if ~any(mask(:))
        error('fiach_processing_steps_viewer:EmptyMask', 'The FIACH mask is empty.');
    end
end


function [design, info] = buildNuisanceDesign(motionFile, regressorFile, nFrames)
    if exist(motionFile, 'file') ~= 2 || exist(regressorFile, 'file') ~= 2
        error('fiach_processing_steps_viewer:RegressorFileMissing', ...
            'The motion and FIACH regressor files are required.');
    end
    motion = readmatrix(motionFile, 'FileType', 'text');
    regressors = readmatrix(regressorFile, 'FileType', 'text');
    motion = motion(:, all(isfinite(motion), 1));
    regressors = regressors(:, all(isfinite(regressors), 1));
    if size(motion, 1) ~= nFrames || size(regressors, 1) ~= nFrames
        error('fiach_processing_steps_viewer:RegressorLengthMismatch', ...
            'Motion/FIACH regressors must contain one row per volume.');
    end
    if size(motion, 2) < 6 || size(regressors, 2) < 6
        error('fiach_processing_steps_viewer:RegressorColumnMismatch', ...
            'Expected at least six motion and six FIACH columns.');
    end
    motion = motion(:, 1:6);
    fiach = regressors(:, end-5:end);
    raw = [motion, fiach];
    centred = raw - mean(raw, 1);
    columnNorm = sqrt(sum(centred .^ 2, 1));
    usable = columnNorm > max(columnNorm) * eps(class(columnNorm)) * nFrames;
    scaled = centred(:, usable) ./ columnNorm(usable);
    [left, singularValues, ~] = svd(scaled, 'econ');
    spectrum = diag(singularValues);
    tolerance = max(size(scaled)) * eps(max(spectrum));
    nuisanceRank = nnz(spectrum > tolerance);
    basis = left(:, 1:nuisanceRank);
    design = struct('raw', raw, 'centredScaled', scaled, 'basis', basis);
    info = struct('nMotion', size(motion, 2), 'nFiach', size(fiach, 2), ...
        'rank', nuisanceRank, 'degreesOfFreedom', nFrames - nuisanceRank - 1, ...
        'preservesTemporalMean', true);
end


function files = makeGlmFilenames(filteredFiles, outputFolder)
    files = cell(size(filteredFiles));
    for index = 1:numel(filteredFiles)
        [~, name, extension] = fileparts(filteredFiles{index});
        files{index} = fullfile(outputFolder, ['glm_' name extension]);
    end
end


function outputFiles = writeGlmSeries(Vfiltered, outputFiles, mask, basis)
    nFrames = numel(Vfiltered);
    maskIndices = find(mask);
    nVoxels = numel(maskIndices);
    filteredMatrix = zeros(nFrames, nVoxels, 'single');
    for frame = 1:nFrames
        volume = double(spm_read_vols(Vfiltered(frame)));
        filteredMatrix(frame, :) = single(volume(maskIndices));
    end

    cleanedMatrix = zeros(size(filteredMatrix), 'single');
    chunkSize = 10000;
    for first = 1:chunkSize:nVoxels
        last = min(first + chunkSize - 1, nVoxels);
        values = double(filteredMatrix(:, first:last));
        cleanedMatrix(:, first:last) = single(values - basis * (basis.' * values));
    end

    for frame = 1:nFrames
        volume = double(spm_read_vols(Vfiltered(frame)));
        volume(maskIndices) = double(cleanedMatrix(frame, :));
        Vout = Vfiltered(frame);
        Vout.fname = outputFiles{frame};
        Vout.dt = [16 0]; % float32 retains the continuous GLM projection.
        Vout.pinfo = [1; 0; 0];
        Vout.descrip = 'FIACH interpolated + motion/FIACH nuisance GLM';
        spm_write_vol(Vout, volume);
    end
    fprintf('Wrote %d nuisance-GLM volumes to %s.\n', nFrames, fileparts(outputFiles{1}));
end


function [metrics, summary] = calculateMetrics(Voriginal, Vfiltered, Vglm, mask, frames)
    nFrames = numel(frames);
    values = zeros(nFrames, 4);
    for frame = 1:nFrames
        original = double(spm_read_vols(Voriginal(frame)));
        filtered = double(spm_read_vols(Vfiltered(frame)));
        glm = double(spm_read_vols(Vglm(frame)));
        originalDifference = filtered(mask) - original(mask);
        glmDifference = glm(mask) - filtered(mask);
        values(frame, :) = [sqrt(mean(originalDifference .^ 2)), ...
            mean(abs(originalDifference)), sqrt(mean(glmDifference .^ 2)), ...
            mean(abs(glmDifference))];
    end
    metrics = array2table([frames(:), values], 'VariableNames', ...
        {'Frame', 'OriginalToFilteredRMSE', 'OriginalToFilteredMAE', ...
        'FilteredToGlmRMSE', 'FilteredToGlmMAE'});
    summary = struct('nFrames', nFrames, 'maskVoxelCount', nnz(mask), ...
        'meanOriginalToFilteredRMSE', mean(metrics.OriginalToFilteredRMSE), ...
        'meanFilteredToGlmRMSE', mean(metrics.FilteredToGlmRMSE), ...
        'maximumFilteredToGlmRMSE', max(metrics.FilteredToGlmRMSE));
end


function writeSummary(file, summary, designInfo, opts)
    fid = fopen(file, 'w');
    if fid < 0
        error('fiach_processing_steps_viewer:SummaryWriteFailed', ...
            'Could not write %s.', file);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'FIACH processing-stage visualization\n');
    fprintf(fid, 'Run folder: %s\n', opts.RunFolder);
    fprintf(fid, 'Volumes: %d\n', summary.nFrames);
    fprintf(fid, 'Mask voxels: %d\n', summary.maskVoxelCount);
    fprintf(fid, 'Motion regressors: %d\n', designInfo.nMotion);
    fprintf(fid, 'FIACH PCA regressors: %d\n', designInfo.nFiach);
    fprintf(fid, 'Nuisance design rank: %d\n', designInfo.rank);
    fprintf(fid, 'Residual degrees of freedom: %d\n', designInfo.degreesOfFreedom);
    fprintf(fid, 'Temporal mean retained: %d\n', designInfo.preservesTemporalMean);
    fprintf(fid, 'Mean original-to-filtered in-mask RMSE: %.12g\n', ...
        summary.meanOriginalToFilteredRMSE);
    fprintf(fid, 'Mean filtered-to-GLM in-mask RMSE: %.12g\n', ...
        summary.meanFilteredToGlmRMSE);
    fprintf(fid, 'Maximum filtered-to-GLM in-mask RMSE: %.12g\n', ...
        summary.maximumFilteredToGlmRMSE);
    clear cleaner;
end


function writeMetricFigure(file, metrics)
    figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 600]);
    cleaner = onCleanup(@() close(figureHandle));
    plot(metrics.Frame, metrics.OriginalToFilteredRMSE, '-', 'LineWidth', 1.2);
    hold on;
    plot(metrics.Frame, metrics.FilteredToGlmRMSE, '-', 'LineWidth', 1.2);
    grid on;
    xlabel('Source volume');
    ylabel('In-mask RMSE');
    title('Magnitude of sequential FIACH processing effects');
    legend({'Original to filtered/interpolated FIACH', ...
        'Filtered/interpolated FIACH to nuisance GLM'}, ...
        'Location', 'best');
    exportgraphics(figureHandle, file, 'Resolution', 170);
    clear cleaner;
end


function figureHandle = openViewer(Voriginal, Vfiltered, Vglm, frameNumbers, mask, metrics, frameRate)
    dimensions = size(mask);
    sliceCounts = squeeze(sum(sum(mask, 1), 2));
    [~, initialSlice] = max(sliceCounts);
    [initialX, initialY] = find(mask(:, :, initialSlice));
    initialVoxel = [round(median(initialX)), round(median(initialY)), initialSlice];

    figureHandle = figure('Name', 'FIACH processing stages', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [40 40 1500 900], 'CloseRequestFcn', @closeViewer);
    originalAxes = axes(figureHandle, 'Position', [0.025 0.54 0.30 0.39]);
    filteredAxes = axes(figureHandle, 'Position', [0.35 0.54 0.30 0.39]);
    glmAxes = axes(figureHandle, 'Position', [0.675 0.54 0.30 0.39]);
    timecourseAxes = axes(figureHandle, 'Position', [0.05 0.19 0.58 0.25]);
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
        'Position', [0.04 0.012 0.93 0.035], 'String', '', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

    plot(metricAxes, metrics.Frame, max(metrics.OriginalToFilteredRMSE, eps), '-', ...
        'LineWidth', 1.1, 'Color', [0.85 0.25 0.15]);
    hold(metricAxes, 'on');
    plot(metricAxes, metrics.Frame, max(metrics.FilteredToGlmRMSE, eps), '-', ...
        'LineWidth', 1.1, 'Color', [0.15 0.55 0.25]);
    set(metricAxes, 'YScale', 'log');
    grid(metricAxes, 'on');
    xlabel(metricAxes, 'Source volume');
    ylabel(metricAxes, 'In-mask RMSE (log)');
    title(metricAxes, 'Processing effect by volume');
    legend(metricAxes, {'Original → FIACH', 'FIACH → nuisance GLM'}, 'Location', 'best');
    metricFrameLine = xline(metricAxes, frameNumbers(1), '-', 'Color', [0.2 0.2 0.2]);

    playbackTimer = timer('ExecutionMode', 'fixedSpacing', 'BusyMode', 'drop', ...
        'Period', 1 / frameRate, 'TimerFcn', @advanceFrame);
    state = struct('frameIndex', 1, 'slice', initialSlice, 'voxel', initialVoxel, ...
        'showOutside', false, 'images', gobjects(1, 3), 'crosshair', gobjects(3, 2), ...
        'timecourseFrameLine', [], 'cachedVoxel', [NaN NaN NaN], ...
        'originalTimecourse', [], 'filteredTimecourse', [], 'glmTimecourse', []);
    guidata(figureHandle, state);
    renderFrame(true);

    function renderFrame(updateTimecourse)
        if ~ishghandle(figureHandle)
            return;
        end
        current = guidata(figureHandle);
        original = double(spm_read_vols(Voriginal(current.frameIndex)));
        filtered = double(spm_read_vols(Vfiltered(current.frameIndex)));
        glm = double(spm_read_vols(Vglm(current.frameIndex)));
        displayVolumes = {original, filtered, glm};
        if ~current.showOutside
            for panel = 1:3
                displayVolumes{panel}(~mask) = 0;
            end
        end
        scaleValues = [original(mask); filtered(mask); glm(mask)];
        limits = robustLimits(scaleValues, 0.01, 0.99);
        axesList = [originalAxes filteredAxes glmAxes];
        titles = {'Original realigned', 'FIACH filtered/interpolated', ...
            'FIACH + nuisance GLM'};
        for panel = 1:3
            slice = displayVolumes{panel}(:, :, current.slice).';
            if ~ishghandle(current.images(panel))
                current.images(panel) = imagesc(axesList(panel), slice, limits);
                current.images(panel).ButtonDownFcn = @imageClicked;
                axis(axesList(panel), 'image');
                axesList(panel).XTick = [];
                axesList(panel).YTick = [];
                colormap(axesList(panel), gray(256));
            else
                current.images(panel).CData = slice;
                clim(axesList(panel), limits);
            end
            title(axesList(panel), sprintf('%s — frame %04d', ...
                titles{panel}, frameNumbers(current.frameIndex)));
        end
        frameSlider.Value = current.frameIndex;
        sliceSlider.Value = current.slice;
        frameText.String = sprintf('%d / %d', current.frameIndex, numel(frameNumbers));
        sliceText.String = sprintf('%d / %d', current.slice, dimensions(3));
        metricFrameLine.Value = frameNumbers(current.frameIndex);
        guidata(figureHandle, current);
        updateCrosshairs();
        if updateTimecourse
            updateTimecoursePlot();
        else
            current = guidata(figureHandle);
            if ~isempty(current.timecourseFrameLine) && ishghandle(current.timecourseFrameLine)
                current.timecourseFrameLine.Value = frameNumbers(current.frameIndex);
            end
        end
        updateStatus(original, filtered, glm);
        drawnow limitrate;
    end

    function updateCrosshairs()
        current = guidata(figureHandle);
        axesList = [originalAxes filteredAxes glmAxes];
        for panel = 1:3
            if ~ishghandle(current.crosshair(panel, 1))
                hold(axesList(panel), 'on');
                current.crosshair(panel, 1) = line(axesList(panel), ...
                    [current.voxel(1) current.voxel(1)], [0.5 dimensions(2) + 0.5], ...
                    'Color', [1 0.75 0], 'LineWidth', 1, 'HitTest', 'off');
                current.crosshair(panel, 2) = line(axesList(panel), ...
                    [0.5 dimensions(1) + 0.5], [current.voxel(2) current.voxel(2)], ...
                    'Color', [1 0.75 0], 'LineWidth', 1, 'HitTest', 'off');
            else
                current.crosshair(panel, 1).XData = [current.voxel(1) current.voxel(1)];
                current.crosshair(panel, 2).YData = [current.voxel(2) current.voxel(2)];
            end
        end
        guidata(figureHandle, current);
    end

    function updateTimecoursePlot()
        current = guidata(figureHandle);
        if ~isequal(current.cachedVoxel, current.voxel)
            coordinate = current.voxel(:);
            current.originalTimecourse = spm_get_data(Voriginal, coordinate);
            current.filteredTimecourse = spm_get_data(Vfiltered, coordinate);
            current.glmTimecourse = spm_get_data(Vglm, coordinate);
            current.cachedVoxel = current.voxel;
        end
        cla(timecourseAxes);
        plot(timecourseAxes, frameNumbers, current.originalTimecourse, '-', ...
            'LineWidth', 1, 'Color', [0.15 0.35 0.75]);
        hold(timecourseAxes, 'on');
        plot(timecourseAxes, frameNumbers, current.filteredTimecourse, '-', ...
            'LineWidth', 1.1, 'Color', [0.85 0.25 0.15]);
        plot(timecourseAxes, frameNumbers, current.glmTimecourse, '-', ...
            'LineWidth', 1.1, 'Color', [0.15 0.55 0.25]);
        current.timecourseFrameLine = xline(timecourseAxes, ...
            frameNumbers(current.frameIndex), '-', 'Color', [0.2 0.2 0.2]);
        grid(timecourseAxes, 'on');
        xlabel(timecourseAxes, 'Source volume');
        ylabel(timecourseAxes, 'Signal intensity');
        title(timecourseAxes, sprintf('Voxel timecourse at [%d %d %d]', current.voxel));
        legend(timecourseAxes, {'Original', 'FIACH filtered/interpolated', ...
            'FIACH + nuisance GLM'}, ...
            'Location', 'best');
        guidata(figureHandle, current);
    end

    function updateStatus(original, filtered, glm)
        current = guidata(figureHandle);
        linear = sub2ind(dimensions, current.voxel(1), current.voxel(2), current.voxel(3));
        location = 'outside mask';
        if mask(linear)
            location = 'inside mask';
        end
        statusText.String = sprintf([ ...
            'Voxel [%d %d %d], %s: original %.6g | FIACH %.6g | GLM %.6g | ' ...
            'interpolation change %.6g | nuisance change %.6g'], current.voxel, location, ...
            original(linear), filtered(linear), glm(linear), ...
            filtered(linear) - original(linear), glm(linear) - filtered(linear));
    end

    function imageClicked(source, ~)
        current = guidata(figureHandle);
        clickedAxes = ancestor(source, 'axes');
        point = clickedAxes.CurrentPoint;
        current.voxel = [min(max(round(point(1, 1)), 1), dimensions(1)), ...
            min(max(round(point(1, 2)), 1), dimensions(2)), current.slice];
        guidata(figureHandle, current);
        renderFrame(true);
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


function limits = robustLimits(values, lowerProbability, upperProbability)
    values = sort(values(isfinite(values)));
    lowerIndex = max(1, round(1 + (numel(values) - 1) * lowerProbability));
    upperIndex = min(numel(values), round(1 + (numel(values) - 1) * upperProbability));
    limits = [values(lowerIndex), values(upperIndex)];
    if limits(1) == limits(2)
        limits = limits + [-0.5 0.5];
    end
end
