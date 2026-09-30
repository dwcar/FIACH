function report = run_overt_0209_validation(varargin)
%RUN_OVERT_0209_VALIDATION Run and compare the packaged R-compatible FIACH.

    repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
    parser = inputParser;
    addParameter(parser, 'SPMRoot', '', @isTextScalar);
    addParameter(parser, 'WorkFolder', fullfile(repositoryRoot, 'work', ...
        'overt_0209_r_validation'), @isTextScalar);
    addParameter(parser, 'OpenViewer', false, @(value) islogical(value) && isscalar(value));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.SPMRoot = char(opts.SPMRoot);
    opts.WorkFolder = char(opts.WorkFolder);
    setup_fiach('SPMRoot', opts.SPMRoot);

    inputFolder = fullfile(repositoryRoot, 'test_data', 'overt_0209', 'input');
    referenceFolder = fullfile(repositoryRoot, 'test_data', 'overt_0209', 'reference_R');
    if exist(opts.WorkFolder, 'dir') ~= 7
        mkdir(opts.WorkFolder);
    end
    inputListing = dir(fullfile(inputFolder, 'rfmri_overt_*.nii'));
    assert(numel(inputListing) == 206, ...
        'run_overt_0209_validation:InputCount', 'Expected 206 packaged input volumes.');
    for index = 1:numel(inputListing)
        destination = fullfile(opts.WorkFolder, inputListing(index).name);
        if exist(destination, 'file') ~= 2
            copyfile(fullfile(inputListing(index).folder, inputListing(index).name), destination);
        end
    end
    rpSource = fullfile(inputFolder, 'rp_fmri_overt_0004.txt');
    rpFile = fullfile(opts.WorkFolder, 'rp_fmri_overt_0004.txt');
    if exist(rpFile, 'file') ~= 2
        copyfile(rpSource, rpFile);
    end

    workListing = dir(fullfile(opts.WorkFolder, 'rfmri_overt_*.nii'));
    [~, order] = sort({workListing.name});
    workListing = workListing(order);
    functional = arrayfun(@(item) fullfile(item.folder, item.name), workListing, ...
        'UniformOutput', false);

    B0 = 3;
    TEms = 26;
    TR = 1.25;
    t = fiach_bold_contrast(B0, TEms);
    fiachResult = fiach(functional, t, TR, ...
        'RP', rpFile, 'MaxGap', 1, 'Freq', 128, 'NMads', 1.96, ...
        'GMMMethod', 'R', 'Overwrite', true, 'MakePlots', false);

    comparisonFolder = fullfile(opts.WorkFolder, 'comparison_to_R');
    comparison = fiach_compare_outputs( ...
        'MatlabFolder', opts.WorkFolder, ...
        'RFolder', referenceFolder, ...
        'Mask', fiachResult.diagnosticFiles.mask, ...
        'OutputFolder', comparisonFolder, ...
        'LeftLabel', 'MATLAB R-compatible', ...
        'RightLabel', 'R reference', ...
        'OpenViewer', opts.OpenViewer);

    assert(strcmp(fiachResult.gmmMethod, 'R'), ...
        'run_overt_0209_validation:WrongGMM', 'Validation did not use the R-compatible GMM.');
    assert(comparison.summary.nFrames == 206, ...
        'run_overt_0209_validation:FrameCount', 'Comparison did not contain 206 frames.');
    assert(comparison.summary.meanMaskCorrelation > 0.9999, ...
        'run_overt_0209_validation:Correlation', ...
        'Mean in-mask MATLAB/R correlation was below 0.9999.');
    assert(comparison.summary.maskMAE < 2, ...
        'run_overt_0209_validation:MAE', ...
        'Pooled in-mask MATLAB/R mean absolute error was at least 2.');

    report = struct('fiach', fiachResult, 'comparison', comparison, ...
        'workFolder', opts.WorkFolder, 'passed', true);
    fprintf(['overt_0209 validation passed: mean in-mask correlation %.9f; ' ...
        'pooled in-mask MAE %.6g.\n'], comparison.summary.meanMaskCorrelation, ...
        comparison.summary.maskMAE);
end


function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end
