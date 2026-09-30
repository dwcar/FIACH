function paths = setup_fiach(varargin)
%SETUP_FIACH Add the packaged MATLAB FIACH implementation and SPM12 to path.

    repositoryRoot = fileparts(mfilename('fullpath'));
    parser = inputParser;
    addParameter(parser, 'SPMRoot', '', @(value) ischar(value) || ...
        (isstring(value) && isscalar(value)));
    parse(parser, varargin{:});
    spmRoot = char(parser.Results.SPMRoot);

    addpath(fullfile(repositoryRoot, 'src'));
    addpath(fullfile(repositoryRoot, 'tools'));
    addpath(fullfile(repositoryRoot, 'tests'));
    addpath(fullfile(repositoryRoot, 'examples'));
    if ~isempty(spmRoot)
        if exist(spmRoot, 'dir') ~= 7
            error('setup_fiach:SPMNotFound', 'SPM directory not found: %s', spmRoot);
        end
        addpath(spmRoot);
    end
    if exist('spm', 'file') ~= 2
        error('setup_fiach:SPMNotOnPath', ...
            'SPM12 was not found. Call setup_fiach(''SPMRoot'', <spm12-folder>).');
    end
    spm('defaults', 'fmri');

    paths = struct('repositoryRoot', repositoryRoot, ...
        'source', fullfile(repositoryRoot, 'src'), ...
        'tools', fullfile(repositoryRoot, 'tools'), ...
        'tests', fullfile(repositoryRoot, 'tests'), ...
        'testData', fullfile(repositoryRoot, 'test_data'));
    fprintf('FIACH MATLAB package ready. Repository: %s\n', repositoryRoot);
end
