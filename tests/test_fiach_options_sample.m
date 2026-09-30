function results = test_fiach_options_sample(varargin)
%TEST_FIACH_OPTIONS_SAMPLE Exercise the new user-mask and GMM options.
%
%   RESULTS = TEST_FIACH_OPTIONS_SAMPLE() first runs the R-compatible
%   validation job, then uses its diagnostic mask as an explicit input mask.
%   If Statistics and Machine Learning Toolbox is available, it also runs
%   the ModifiedICH GMM selector. All output is contained in the existing
%   dedicated sample-validation directory.

    parser = inputParser;
    addParameter(parser, 'RunModifiedICH', true, @(x) islogical(x) && isscalar(x));
    parse(parser, varargin{:});
    opts = parser.Results;

    results = struct();
    results.automatic = test_fiach_matlab_sample();
    userMask = results.automatic.diagnosticFiles.mask;
    functional = results.automatic.input{1};
    t = fiach_bold_contrast(1.5, 30);

    results.userMask = fiach(functional, t, 2.16, ...
        'Mask', userMask, 'UseUserMask', true, 'GMMMethod', 'R', ...
        'Overwrite', true, 'MakePlots', false);
    assert(results.userMask.useUserMask, ...
        'test_fiach_options_sample:UserMaskIgnored', 'The explicit user mask was not used.');
    assert(strcmp(results.userMask.maskSource, userMask), ...
        'test_fiach_options_sample:MaskSource', 'The reported mask path differs from the requested mask.');
    assert(strcmp(results.userMask.gmmMethod, 'R'), ...
        'test_fiach_options_sample:RMethod', 'The requested R GMM was not used.');

    results.modifiedICH = [];
    if opts.RunModifiedICH && exist('fitgmdist', 'file') == 2 && exist('statset', 'file') == 2
        results.modifiedICH = fiach(functional, t, 2.16, ...
            'Mask', userMask, 'GMMMethod', 'ModifiedICH', ...
            'GMMPosteriorThreshold', .5, 'Overwrite', true, 'MakePlots', false);
        assert(strcmp(results.modifiedICH.gmmMethod, 'ModifiedICH'), ...
            'test_fiach_options_sample:ICHMethod', 'The requested ModifiedICH GMM was not used.');
    elseif opts.RunModifiedICH
        fprintf(['ModifiedICH validation skipped because the Statistics and Machine ' ...
            'Learning Toolbox is not available.\n']);
    end

    fprintf('FIACH option validation passed. User mask: %s\n', userMask);
end
