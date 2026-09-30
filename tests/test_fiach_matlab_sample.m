function result = test_fiach_matlab_sample(varargin)
%TEST_FIACH_MATLAB_SAMPLE Run fiach.m on the upstream FIACH toy dataset.
%
%   RESULT = TEST_FIACH_MATLAB_SAMPLE() uses FIACH's motion_ex.nii.gz
%   example copied into fiach_matlab_sample_validation. It uses the same
%   B0, TE, and TR as the R package example (1.5 T, 30 ms, and 2.16 s).
%   Existing files in that dedicated validation folder may be replaced.

    thisFolder = fileparts(mfilename('fullpath'));
    defaultSpm = '';
    defaultSampleDirectory = fullfile(thisFolder, 'fiach_matlab_sample_validation');

    parser = inputParser;
    addParameter(parser, 'SPMRoot', defaultSpm, @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'SampleDirectory', defaultSampleDirectory, @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.SPMRoot = char(opts.SPMRoot);
    opts.SampleDirectory = char(opts.SampleDirectory);

    setup_fiach('SPMRoot', opts.SPMRoot);

    compressed = fullfile(opts.SampleDirectory, 'motion_ex.nii.gz');
    functional = fullfile(opts.SampleDirectory, 'motion_ex.nii');
    assert(exist(compressed, 'file') == 2 || exist(functional, 'file') == 2, ...
        'test_fiach_matlab_sample:SampleNotFound', 'Could not find motion_ex.nii(.gz) in %s.', opts.SampleDirectory);
    if exist(functional, 'file') ~= 2
        gunzip(compressed, opts.SampleDirectory);
    end

    t = fiach_bold_contrast(1.5, 30);
    result = fiach(functional, t, 2.16, 'Overwrite', true, 'MakePlots', false);

    assert(exist(result.functionalFiles{1}, 'file') == 2, ...
        'test_fiach_matlab_sample:OutputMissing', 'Corrected functional output was not written.');
    assert(exist(result.diagnosticFiles.mask, 'file') == 2, ...
        'test_fiach_matlab_sample:MaskMissing', 'Mask output was not written.');
    assert(exist(result.diagnosticFiles.noiseMask, 'file') == 2, ...
        'test_fiach_matlab_sample:NoiseMaskMissing', 'Noise-mask output was not written.');
    assert(exist(result.noiseBasisFile, 'file') == 2, ...
        'test_fiach_matlab_sample:RegressorMissing', 'Noise-regressor output was not written.');
    assert(numel(spm_vol(result.functionalFiles{1})) == numel(spm_vol(functional)), ...
        'test_fiach_matlab_sample:VolumeMismatch', 'Output has a different number of volumes than input.');
    assert(strcmp(result.gmmMethod, 'R') && ~result.useUserMask, ...
        'test_fiach_matlab_sample:DefaultOptions', 'The sample test should use the R GMM and automatic mask by default.');

    fprintf('FIACH MATLAB sample validation passed. Output: %s\n', result.functionalFiles{1});
end
