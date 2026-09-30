function result = fiach_compare_pca_regressors(varargin)
%FIACH_COMPARE_PCA_REGRESSORS Compare R-GMM and ModifiedICH FIACH PCA scores.
%
%   RESULT = FIACH_COMPARE_PCA_REGRESSORS() reads the overt_0209 outputs,
%   aligns ModifiedICH PCA components to the R-compatible components by
%   maximum absolute temporal correlation, corrects arbitrary PCA signs,
%   and writes a three-panel image plus aligned numeric outputs.

    thisFolder = fileparts(mfilename('fullpath'));
    repositoryRoot = fileparts(thisFolder);
    defaultRFile = fullfile(repositoryRoot, 'work', 'overt_0209_r_validation', ...
        'noise_basis6.txt');
    defaultICHFile = fullfile(repositoryRoot, 'work', 'overt_0209_ich', ...
        'noise_basis6.txt');
    defaultOutputFolder = fullfile(repositoryRoot, 'work', ...
        'fiach_gmm_comparison_overt_0209');

    parser = inputParser;
    addParameter(parser, 'RCompatibleFile', defaultRFile, @isTextScalar);
    addParameter(parser, 'ModifiedICHFile', defaultICHFile, @isTextScalar);
    addParameter(parser, 'OutputFolder', defaultOutputFolder, @isTextScalar);
    addParameter(parser, 'Visible', true, @(value) islogical(value) && isscalar(value));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.RCompatibleFile = char(opts.RCompatibleFile);
    opts.ModifiedICHFile = char(opts.ModifiedICHFile);
    opts.OutputFolder = char(opts.OutputFolder);

    rScores = readPcaScores(opts.RCompatibleFile);
    ichScores = readPcaScores(opts.ModifiedICHFile);
    if size(rScores, 1) ~= size(ichScores, 1)
        error('fiach_compare_pca_regressors:LengthMismatch', ...
            'The regressor files contain different numbers of time points.');
    end

    nComponents = min(size(rScores, 2), size(ichScores, 2));
    rScores = rScores(:, 1:nComponents);
    ichScores = ichScores(:, 1:nComponents);
    correlationMatrix = componentCorrelations(rScores, ichScores);
    [assignment, signs] = alignComponents(correlationMatrix);
    alignedICH = ichScores(:, assignment) .* signs;
    difference = alignedICH - rScores;

    componentCorrelation = zeros(nComponents, 1);
    componentRmse = zeros(nComponents, 1);
    for component = 1:nComponents
        componentCorrelation(component) = correlationMatrix(component, assignment(component)) * signs(component);
        componentRmse(component) = sqrt(mean(difference(:, component) .^ 2));
    end
    componentCorrelation = abs(componentCorrelation);

    if exist(opts.OutputFolder, 'dir') ~= 7
        mkdir(opts.OutputFolder);
    end
    figureFile = fullfile(opts.OutputFolder, 'pca_regressor_comparison.png');
    alignedFile = fullfile(opts.OutputFolder, 'pca_regressors_aligned.csv');
    summaryFile = fullfile(opts.OutputFolder, 'pca_regressor_component_summary.csv');

    variableNames = [{'Frame'}, cellstr(compose('R_PC%d', 1:nComponents)), ...
        cellstr(compose('ICH_PC%d_aligned', 1:nComponents)), ...
        cellstr(compose('Difference_PC%d', 1:nComponents))];
    alignedTable = array2table([(1:size(rScores, 1)).', rScores, alignedICH, difference], ...
        'VariableNames', variableNames);
    writetable(alignedTable, alignedFile);
    summaryTable = table((1:nComponents).', assignment(:), signs(:), ...
        componentCorrelation, componentRmse, ...
        'VariableNames', {'RComponent', 'MatchedICHComponent', 'AppliedSign', 'Correlation', 'RMSE'});
    writetable(summaryTable, summaryFile);

    visibility = 'off';
    if opts.Visible
        visibility = 'on';
    end
    figureHandle = figure('Name', 'FIACH PCA regressors: R GMM versus ModifiedICH', ...
        'NumberTitle', 'off', 'Color', 'w', 'Visible', visibility, ...
        'Position', [80 120 1450 620]);
    layout = tiledlayout(figureHandle, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    commonLimit = max(abs([rScores(:); alignedICH(:)]));
    commonLimit = max(commonLimit, eps);
    differenceLimit = max(abs(difference(:)));
    differenceLimit = max(differenceLimit, eps);

    rAxes = nexttile(layout);
    imagesc(rAxes, rScores, [-commonLimit commonLimit]);
    title(rAxes, 'R-compatible GMM PCA regressors');
    configureAxes(rAxes, nComponents);

    ichAxes = nexttile(layout);
    imagesc(ichAxes, alignedICH, [-commonLimit commonLimit]);
    title(ichAxes, 'ModifiedICH PCA regressors (aligned)');
    configureAxes(ichAxes, nComponents);

    differenceAxes = nexttile(layout);
    imagesc(differenceAxes, difference, [-differenceLimit differenceLimit]);
    title(differenceAxes, sprintf('ModifiedICH - R-compatible (mean RMSE %.4g)', ...
        mean(componentRmse)));
    configureAxes(differenceAxes, nComponents);

    map = divergingMap(256);
    colormap(rAxes, map);
    colormap(ichAxes, map);
    colormap(differenceAxes, map);
    colorbar(rAxes);
    colorbar(ichAxes);
    colorbar(differenceAxes);
    title(layout, ['PCA components matched by maximum absolute correlation; ' ...
        'arbitrary component signs corrected']);
    exportgraphics(figureHandle, figureFile, 'Resolution', 180);

    result = struct();
    result.rCompatibleScores = rScores;
    result.modifiedICHScoresAligned = alignedICH;
    result.difference = difference;
    result.assignment = assignment;
    result.signs = signs;
    result.componentCorrelation = componentCorrelation;
    result.componentRmse = componentRmse;
    result.figure = figureHandle;
    result.figureFile = figureFile;
    result.alignedFile = alignedFile;
    result.summaryFile = summaryFile;

    fprintf('Mean aligned PCA correlation %.6f; mean component RMSE %.6g.\n', ...
        mean(componentCorrelation), mean(componentRmse));
end


function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end


function scores = readPcaScores(file)
    if exist(file, 'file') ~= 2
        error('fiach_compare_pca_regressors:FileMissing', 'File not found: %s', file);
    end
    values = readmatrix(file, 'FileType', 'text');
    values = values(:, all(isfinite(values), 1));
    if size(values, 2) < 6
        error('fiach_compare_pca_regressors:TooFewColumns', ...
            'Expected at least six numeric columns in %s.', file);
    end
    % With an RP file, FIACH writes six motion columns followed by six PCs.
    scores = values(:, end-5:end);
end


function correlations = componentCorrelations(left, right)
    nComponents = size(left, 2);
    correlations = zeros(nComponents);
    for leftIndex = 1:nComponents
        for rightIndex = 1:nComponents
            pair = corrcoef(left(:, leftIndex), right(:, rightIndex));
            correlations(leftIndex, rightIndex) = pair(1, 2);
        end
    end
end


function [assignment, signs] = alignComponents(correlations)
    nComponents = size(correlations, 1);
    candidates = perms(1:nComponents);
    score = zeros(size(candidates, 1), 1);
    rows = 1:nComponents;
    for index = 1:size(candidates, 1)
        linear = sub2ind(size(correlations), rows, candidates(index, :));
        score(index) = sum(abs(correlations(linear)));
    end
    [~, best] = max(score);
    assignment = candidates(best, :);
    linear = sub2ind(size(correlations), rows, assignment);
    signs = sign(correlations(linear));
    signs(signs == 0) = 1;
end


function configureAxes(axesHandle, nComponents)
    xlabel(axesHandle, 'PCA component');
    ylabel(axesHandle, 'Frame');
    axesHandle.XTick = 1:nComponents;
    axesHandle.YDir = 'normal';
    axis(axesHandle, 'tight');
end


function map = divergingMap(n)
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
