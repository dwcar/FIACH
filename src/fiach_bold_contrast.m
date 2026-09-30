function contrast = fiach_bold_contrast(B0, TEms, varargin)
%FIACH_BOLD_CONTRAST MATLAB port of FIACH::boldContrast.
%
%   T = FIACH_BOLD_CONTRAST(B0, TEMS) returns the expected BOLD signal
%   difference at the supplied echo time (milliseconds). This is the value
%   used as input T in FIAch's deterministic bad-data threshold:
%
%       results = fiach('realigned_bold.nii', ...
%           fiach_bold_contrast(3, 26), 1.25);
%
%   Name-value options reproduce R/boldContrast.R:
%     'Random'  random-vessel model (true)
%     'Alpha'   CBV/CBF exponent (.38)
%     'Hct'     haematocrit (.4)
%     'CBFBase' baseline CBF (55)
%     'Chi'     susceptibility difference (1.8e-7)
%     'EBase'   baseline oxygen extraction (.4)
%     'W0'      proton gyromagnetic ratio (2.675e8 rad/s/T)
%     'EAct'    activated oxygen extraction (.1)
%     'CBFAct'  activated CBF (2*CBFBase)
%
%   This mirrors the upstream helper's numerical model; it does not produce
%   a figure. TEms must be an integer from 0 to 250 because the R code
%   samples its signal model on that one-millisecond grid.
%
% DC to check at some point accuracy of short TRs - is T1 accounted for?

    validateattributes(B0, {'numeric'}, {'scalar', 'real', 'finite', '>', 0}, mfilename, 'B0', 1);
    validateattributes(TEms, {'numeric'}, {'scalar', 'integer', '>=', 0, '<=', 250}, mfilename, 'TEms', 2);

    parser = inputParser;
    addParameter(parser, 'Random', true, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'Alpha', .38, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'Hct', .4, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'CBFBase', 55, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'Chi', 1.8e-7, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'EBase', .4, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'W0', 2.675e8, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'EAct', .1, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'CBFAct', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(parser, varargin{:});
    opts = parser.Results;
    if isempty(opts.CBFAct)
        opts.CBFAct = 2 * opts.CBFBase;
    end

    cbvBase = .8 * opts.CBFBase ^ opts.Alpha;
    volumeBase = cbvBase / 100;
    cbvAct = .8 * opts.CBFAct ^ opts.Alpha;
    volumeAct = cbvAct / 100;
    r2Intrinsic = 1.74 * B0 + 7.77;
    chiSI = opts.Chi * 4 * pi;
    w0SI = opts.W0 / (2 * pi);

    if opts.Random
        geometry = 4 * pi / 3;
    else
        geometry = 2 * pi;
    end
    frequencyBase = w0SI * B0 * chiSI * opts.Hct * opts.EBase * geometry;
    r2StarBase = volumeBase * frequencyBase + r2Intrinsic;
    frequencyAct = w0SI * B0 * chiSI * opts.Hct * opts.EAct * geometry;
    r2StarAct = r2Intrinsic + volumeAct * frequencyAct;

    t2StarBaseMs = 1000 / r2StarBase;
    t2StarActMs = 1000 / r2StarAct;
    timeMs = 0:250;
    baseline = 100 * exp(-timeMs / t2StarBaseMs);
    activated = 100 * exp(-timeMs / t2StarActMs);
    signalDifference = activated - baseline;
    contrast = signalDifference(TEms + 1);
end
