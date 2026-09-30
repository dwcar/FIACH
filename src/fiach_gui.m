function fig = fiach_gui()
%FIACH_GUI Interactive front end for the MATLAB/SPM FIACH implementation.
%
%   FIACH_GUI opens a single-subject FIACH job window modelled on the
%   original R/Shiny interface. It collects the functional image(s), B0,
%   TE, TR, FIACH parameters, optional motion parameters, optional aligned
%   user brain mask, and the desired GMM segmentation method. The Save
%   batch script button writes the current job as reproducible MATLAB code.
%
%   A single functional 4-D NIfTI file represents one run. Multiple files
%   may instead be selected when each represents one 3-D time point.

    if exist('uifigure', 'file') ~= 2
        error('fiach_gui:ModernMATLABRequired', ...
            'fiach_gui requires MATLAB support for uifigure (R2016a or later).');
    end

    fig = uifigure('Name', 'FIACH for MATLAB/SPM', ...
        'Position', [100 80 930 680], 'Resize', 'on');

    % ---- Functional inputs -------------------------------------------------
    uilabel(fig, 'Position', [25 620 405 22], 'Text', ...
        'Functional image(s): one 4-D file or sequential 3-D files', ...
        'FontWeight', 'bold');
    functionalList = uilistbox(fig, 'Position', [25 435 405 180], ...
        'Items', {}, 'Multiselect', 'on');
    uibutton(fig, 'push', 'Text', 'Add functional image(s)...', ...
        'Position', [25 398 195 28], 'ButtonPushedFcn', @onAddFunctional);
    uibutton(fig, 'push', 'Text', 'Remove selected', ...
        'Position', [235 398 195 28], 'ButtonPushedFcn', @onRemoveFunctional);

    uilabel(fig, 'Position', [25 363 405 20], 'Text', ...
        'Realignment parameters (optional: six-column RP text file)');
    rpField = uieditfield(fig, 'text', 'Position', [25 333 320 27]);
    uibutton(fig, 'push', 'Text', 'Browse...', 'Position', [355 333 75 27], ...
        'ButtonPushedFcn', @onBrowseRP);

    userMaskCheck = uicheckbox(fig, 'Position', [25 296 405 24], ...
        'Text', 'Use a supplied brain mask (otherwise use FIACH automatic mask)', ...
        'Value', false, 'ValueChangedFcn', @onUserMaskChanged);
    maskField = uieditfield(fig, 'text', 'Position', [25 263 320 27], ...
        'Enable', 'off');
    maskBrowseButton = uibutton(fig, 'push', 'Text', 'Browse...', ...
        'Position', [355 263 75 27], 'Enable', 'off', ...
        'ButtonPushedFcn', @onBrowseMask);
    uilabel(fig, 'Position', [25 233 405 23], 'Text', ...
        'Mask must be one 3-D image already aligned to the functional grid.', ...
        'FontAngle', 'italic', 'FontSize', 10);

    % ---- The controls exposed by the original R GUI -----------------------
    settingsPanel = uipanel(fig, 'Title', 'FIACH settings', ...
        'Position', [470 230 435 415]);
    makeSettingLabel(settingsPanel, 'Magnetic field B0 (T)', [20 350 185 22]);
    b0Field = uieditfield(settingsPanel, 'numeric', 'Position', [235 350 155 25], 'Value', 3);
    makeSettingLabel(settingsPanel, 'Repetition time TR (s)', [20 310 185 22]);
    trField = uieditfield(settingsPanel, 'numeric', 'Position', [235 310 155 25], 'Value', 2);
    makeSettingLabel(settingsPanel, 'Echo time TE (ms)', [20 270 185 22]);
    teField = uieditfield(settingsPanel, 'numeric', 'Position', [235 270 155 25], 'Value', 30);
    makeSettingLabel(settingsPanel, 'No. MADS', [20 230 185 22]);
    nmadsField = uieditfield(settingsPanel, 'numeric', 'Position', [235 230 155 25], 'Value', 1.96);
    makeSettingLabel(settingsPanel, 'High-pass period (s)', [20 190 185 22]);
    frequencyField = uieditfield(settingsPanel, 'numeric', 'Position', [235 190 155 25], 'Value', 128);
    makeSettingLabel(settingsPanel, 'Maximum gap (volumes)', [20 150 185 22]);
    maxGapField = uieditfield(settingsPanel, 'numeric', 'Position', [235 150 155 25], 'Value', 1);
    makeSettingLabel(settingsPanel, 'GMM segmentation', [20 105 185 22]);
    gmmDropDown = uidropdown(settingsPanel, 'Position', [235 105 155 25], ...
        'Items', {'R-compatible', 'FIACH_modified_ICH'}, ...
        'ItemsData', {'R', 'ModifiedICH'}, 'Value', 'R', ...
        'ValueChangedFcn', @onGMMChanged);
    makeSettingLabel(settingsPanel, 'ICH posterior threshold', [20 65 185 22]);
    posteriorField = uieditfield(settingsPanel, 'numeric', ...
        'Position', [235 65 155 25], 'Value', .5, 'Enable', 'off');
    uilabel(settingsPanel, 'Position', [20 24 375 27], ...
        'Text', 'R-compatible is the default. ModifiedICH requires Statistics Toolbox.', ...
        'FontAngle', 'italic', 'FontSize', 10);

    overwriteCheck = uicheckbox(fig, 'Position', [470 195 205 24], ...
        'Text', 'Replace existing outputs', 'Value', false);
    plotCheck = uicheckbox(fig, 'Position', [685 195 220 24], ...
        'Text', 'Write GMM plot', 'Value', false);

    % ---- Status and actions ------------------------------------------------
    uilabel(fig, 'Position', [25 194 380 22], 'Text', 'Status', 'FontWeight', 'bold');
    statusArea = uitextarea(fig, 'Position', [25 55 880 135], 'Editable', 'off', ...
        'Value', {'Choose a functional image, set the acquisition parameters, then run FIACH.'});
    uibutton(fig, 'push', 'Text', 'Help', 'Position', [25 18 100 28], ...
        'ButtonPushedFcn', @onHelp);
    uibutton(fig, 'push', 'Text', 'Save batch script...', 'Position', [135 18 170 28], ...
        'ButtonPushedFcn', @onSaveBatch);
    runButton = uibutton(fig, 'push', 'Text', 'Run FIACH', ...
        'Position', [735 18 170 28], 'FontWeight', 'bold', ...
        'ButtonPushedFcn', @onRun);

    onUserMaskChanged([], []);
    onGMMChanged([], []);

    function makeSettingLabel(parent, text, position)
        uilabel(parent, 'Position', position, 'Text', text);
    end

    function onAddFunctional(~, ~)
        [names, folder] = uigetfile( ...
            {'*.nii;*.nii.gz;*.img;*.hdr', 'Functional images (*.nii, *.nii.gz, *.img, *.hdr)'; ...
             '*.*', 'All files'}, ...
            'Select a 4-D functional image or sequential 3-D images', ...
            'MultiSelect', 'on');
        if isequal(names, 0)
            return;
        end
        if ischar(names)
            names = {names};
        end
        newFiles = cellfun(@(name) fullfile(folder, name), names, 'UniformOutput', false);
        functionalList.Items = [functionalList.Items(:); newFiles(:)];
        functionalList.Value = newFiles{end};
        writeStatus(sprintf('Added %d functional image(s).', numel(newFiles)));
    end

    function onRemoveFunctional(~, ~)
        selected = functionalList.Value;
        if isempty(selected)
            return;
        end
        if ischar(selected)
            selected = {selected};
        end
        keep = ~ismember(functionalList.Items, selected);
        functionalList.Items = functionalList.Items(keep);
        if ~isempty(functionalList.Items)
            functionalList.Value = functionalList.Items{1};
        end
        writeStatus(sprintf('Removed %d functional image(s).', numel(selected)));
    end

    function onBrowseRP(~, ~)
        [name, folder] = uigetfile({'*.txt;*.tsv;*.csv', 'Motion parameter text files'; '*.*', 'All files'}, ...
            'Select realignment parameters');
        if ~isequal(name, 0)
            rpField.Value = fullfile(folder, name);
        end
    end

    function onUserMaskChanged(~, ~)
        enabled = logical(userMaskCheck.Value);
        maskField.Enable = onOff(enabled);
        maskBrowseButton.Enable = onOff(enabled);
    end

    function onBrowseMask(~, ~)
        [name, folder] = uigetfile({'*.nii;*.nii.gz;*.img;*.hdr', 'Brain-mask images'; '*.*', 'All files'}, ...
            'Select an aligned 3-D brain mask');
        if ~isequal(name, 0)
            maskField.Value = fullfile(folder, name);
        end
    end

    function onGMMChanged(~, ~)
        posteriorField.Enable = onOff(strcmp(gmmDropDown.Value, 'ModifiedICH'));
    end

    function onRun(~, ~)
        try
            [functional, t, tr, arguments] = currentJob();
        catch exception
            showError(exception);
            return;
        end
        runButton.Enable = 'off';
        fig.Pointer = 'watch';
        writeStatus('Running FIACH. This can take several minutes for a full fMRI run...');
        drawnow;
        cleanup = onCleanup(@() restoreRunControls(fig, runButton));
        try
            result = fiach(functional, t, tr, arguments{:});
            fig.UserData = result;
            writeStatus({ ...
                'FIACH completed successfully.', ...
                ['Corrected image: ' result.functionalFiles{1}], ...
                ['Diagnostics: ' result.diagnosticDirectory], ...
                sprintf('Mask: %s | GMM: %s', result.maskSource, result.gmmMethod)});
        catch exception
            writeStatus(['FIACH failed: ' exception.message]);
            showError(exception);
        end
    end

    function restoreRunControls(figureHandle, buttonHandle)
        if isvalid(figureHandle)
            figureHandle.Pointer = 'arrow';
            buttonHandle.Enable = 'on';
        end
    end

    function onSaveBatch(~, ~)
        try
            [functional, ~, tr, arguments, settings] = currentJob();
        catch exception
            showError(exception);
            return;
        end
        [name, folder] = uiputfile('fiach_batch.m', 'Save FIACH batch script');
        if isequal(name, 0)
            return;
        end
        scriptFile = fullfile(folder, name);
        try
            writeBatchScript(scriptFile, functional, tr, arguments, settings);
            writeStatus(['Saved reproducible FIACH batch script: ' scriptFile]);
        catch exception
            showError(exception);
        end
    end

    function [functional, t, tr, arguments, settings] = currentJob()
        functional = functionalList.Items;
        if isempty(functional)
            error('fiach_gui:MissingFunctional', 'Select at least one functional image.');
        end
        settings = struct();
        settings.B0 = validNumber(b0Field.Value, 'B0', 0, false);
        tr = validNumber(trField.Value, 'TR', 0, false);
        settings.TEms = validNumber(teField.Value, 'TE', 0, true);
        settings.NMads = validNumber(nmadsField.Value, 'No. MADS', 0, true);
        settings.Freq = validNumber(frequencyField.Value, 'High-pass period', 0, false);
        settings.MaxGap = validNumber(maxGapField.Value, 'Maximum gap', 1, true);
        if settings.MaxGap ~= round(settings.MaxGap)
            error('fiach_gui:InvalidMaxGap', 'Maximum gap must be a whole number of volumes.');
        end
        settings.GMMMethod = gmmDropDown.Value;
        settings.GMMPosteriorThreshold = validNumber(posteriorField.Value, ...
            'ICH posterior threshold', 0, false);
        if settings.GMMPosteriorThreshold >= 1
            error('fiach_gui:InvalidPosteriorThreshold', ...
                'ICH posterior threshold must be greater than 0 and less than 1.');
        end
        t = fiach_bold_contrast(settings.B0, settings.TEms);
        arguments = {'MaxGap', settings.MaxGap, 'Freq', settings.Freq, ...
            'NMads', settings.NMads, 'GMMMethod', settings.GMMMethod, ...
            'GMMPosteriorThreshold', settings.GMMPosteriorThreshold, ...
            'Overwrite', logical(overwriteCheck.Value), ...
            'MakePlots', logical(plotCheck.Value)};
        if userMaskCheck.Value
            if isempty(strtrim(maskField.Value))
                error('fiach_gui:MissingMask', ...
                    'Choose a mask file or untick "Use a supplied brain mask".');
            end
            arguments = [arguments, {'UseUserMask', true, 'Mask', maskField.Value}];
        else
            arguments = [arguments, {'UseUserMask', false}];
        end
        if ~isempty(strtrim(rpField.Value))
            arguments = [arguments, {'RP', rpField.Value}];
        end
    end

    function value = validNumber(value, label, minimum, allowEqual)
        if ~isscalar(value) || ~isfinite(value) || ...
                (allowEqual && value < minimum) || (~allowEqual && value <= minimum)
            comparison = 'greater than';
            if allowEqual
                comparison = 'greater than or equal to';
            end
            error('fiach_gui:InvalidNumber', '%s must be finite and %s %g.', ...
                label, comparison, minimum);
        end
    end

    function writeBatchScript(scriptFile, functional, tr, arguments, settings)
        fid = fopen(scriptFile, 'w');
        if fid < 0
            error('fiach_gui:WriteFailed', 'Could not write %s.', scriptFile);
        end
        cleaner = onCleanup(@() fclose(fid));
        savedAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
        fprintf(fid, '%% FIACH batch script saved by fiach_gui on %s\n', savedAt);
        fprintf(fid, '%% Add SPM12 and the folder containing fiach.m to the MATLAB path first.\n\n');
        fprintf(fid, 'functional = { ...\n');
        for i = 1:numel(functional)
            comma = ',';
            if i == numel(functional)
                comma = '';
            end
            fprintf(fid, '    %s%s ...\n', matlabQuote(functional{i}), comma);
        end
        fprintf(fid, '    };\n');
        fprintf(fid, 'B0 = %.15g;\nTEms = %.15g;\nTR = %.15g;\n', ...
            settings.B0, settings.TEms, tr);
        fprintf(fid, 't = fiach_bold_contrast(B0, TEms);\n');
        fprintf(fid, 'result = fiach(functional, t, TR, ...\n');
        for i = 1:2:numel(arguments)
            name = arguments{i};
            value = arguments{i+1};
            if ischar(value)
                valueText = matlabQuote(value);
            elseif islogical(value)
                valueText = lower(mat2str(value));
            else
                valueText = mat2str(value);
            end
            comma = ',';
            if i == numel(arguments) - 1
                comma = '';
            end
            fprintf(fid, '    %s, %s%s ...\n', matlabQuote(name), valueText, comma);
        end
        fprintf(fid, '    );\n');
        fprintf(fid, 'disp(result)\n');
    end

    function quoted = matlabQuote(value)
        quoted = ['''' strrep(char(value), '''', '''''') ''''];
    end

    function value = onOff(logicalValue)
        if logicalValue
            value = 'on';
        else
            value = 'off';
        end
    end

    function writeStatus(message)
        if ischar(message)
            message = {message};
        end
        statusArea.Value = message;
    end

    function showError(exception)
        if isvalid(fig)
            uialert(fig, exception.message, 'FIACH could not start or complete');
        end
    end

    function onHelp(~, ~)
        helpText = sprintf([ ...
            'Select one 4-D functional image, or multiple sequential 3-D images.\n\n' ...
            'B0, TE, and TR are used to calculate the FIACH deterministic BOLD threshold. ' ...
            'No. MADS, high-pass period, and maximum gap correspond to the controls in the original R GUI.\n\n' ...
            'If a user mask is selected it must already have the same dimensions and affine as the functional image. ' ...
            'If not selected, FIACH creates its own mask.\n\n' ...
            'R-compatible GMM is the default. FIACH_modified_ICH reproduces the histogram-initialised, ' ...
            'posterior-based GMM choice from run_fiach_DC_modified_ICH.m and requires the Statistics and Machine Learning Toolbox.']);
        uialert(fig, helpText, 'FIACH GUI help');
    end
end
