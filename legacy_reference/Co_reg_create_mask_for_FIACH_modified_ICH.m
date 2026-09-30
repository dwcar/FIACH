

%% ============================================================
%  Co_reg_create_mask_for_FIACH.m
%  Automatically detects subject ID (e.g. S04)
% ============================================================

clear; clc;
spm('defaults','fmri');
spm_jobman('initcfg');

fprintf('\n--- Running coregistration and mask creation for FIACH ---\n');

%% ------------------------------------------------------------------------
%  Detect subject ID (robust search up the directory tree)
% -------------------------------------------------------------------------
currDir = 'S:\ICH_DCNU_ESM_fMRI\To be analysed Prisma\RaHa\language\T1_FLASH'%pwd;
cd(currDir)
subjectID = '';

parts = strsplit(currDir, filesep);
for i = numel(parts):-1:1
    if ~isempty(regexp(parts{i}, '^S\d{2}$', 'once'))
        subjectID = parts{i};
        break;
    end
end

if isempty(subjectID)
    warning('⚠️  Could not detect subject folder (expected S##). Using "UnknownSubject".');
    subjectID = 'UnknownSubject';
end

fprintf('Detected subject ID: %s\n', subjectID);

%% ------------------------------------------------------------------------
%  File setup
% -------------------------------------------------------------------------
ref_epi = {'..\VG_OV_temp\meanfmri_overt_0004.nii'};  % EPI reference (target)
src_t1  = {'t1_mprage.nii'};             % skull-stripped T1 (source)
segs    = {                        % tissue probability maps
    'c1t1_mprage.nii'
    'c2t1_mprage.nii'
    'c3t1_mprage.nii'
};

%% ------------------------------------------------------------------------
%  Step 1: Coregister UNI → meanFunctional
% -------------------------------------------------------------------------
fprintf('Step 1: Coregister UNI → meanFunctional...\n');

clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.estimate.ref    = ref_epi;
matlabbatch{1}.spm.spatial.coreg.estimate.source = src_t1;
matlabbatch{1}.spm.spatial.coreg.estimate.other  = segs;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
spm_jobman('run',matlabbatch);
clear matlabbatch;

%% ------------------------------------------------------------------------
%  Step 2: Reslice the segmentations into EPI space
% -------------------------------------------------------------------------
fprintf('Step 2: Reslice c1/c2/c3 to EPI grid...\n');

clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.write.ref = ref_epi;
matlabbatch{1}.spm.spatial.coreg.write.source = segs;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch;

%% ------------------------------------------------------------------------
%  Step 3: Combine rc1/rc2/rc3 to create brain mask
% -------------------------------------------------------------------------
fprintf('Step 3: Combine rc1/rc2/rc3 to create brain mask...\n');

V1 = spm_vol('rc1t1_mprage.nii');
V2 = spm_vol('rc2t1_mprage.nii');
V3 = spm_vol('rc3t1_mprage.nii');
c1 = spm_read_vols(V1);
c2 = spm_read_vols(V2);
c3 = spm_read_vols(V3);

mask = (c1 + c2 + c3) > 0.3;   % liberal inclusion threshold
mask = imfill(mask,'holes');

Vmask = V1;
Vmask.fname = sprintf('rfBrainMask_%s.nii', subjectID);
Vmask.dt = [2 0];
spm_write_vol(Vmask, mask);

fprintf('\n✅ Mask saved as rfBrainMask_%s.nii\n', subjectID);

%% ------------------------------------------------------------------------
%  Optional QC overlay
% -------------------------------------------------------------------------
try
    fprintf('Launching SPM overlay for visual check...\n');
    spm_image('Display','meanFunctional.nii');
    spm_image('AddOverlay',sprintf('rfBrainMask_%s.nii', subjectID));
catch
    fprintf('SPM overlay skipped (headless mode or display issue)\n');
end

fprintf('\n--- Done. Check overlay alignment visually ---\n');









% %% ============================================================
% %  Co_reg_create_mask_for_FIACH.m
% %  Coregisters UNI → meanFunctional and builds EPI-space brain mask
% %  for use with FIACH preprocessing
% %  ============================================================
% 
% clear; clc;
% spm('defaults','fmri');
% spm_jobman('initcfg');
% 
% fprintf('\n--- Running coregistration and mask creation for FIACH ---\n');
% 
% %% ------------------------------------------------------------------------
% %  File setup
% % -------------------------------------------------------------------------
% ref_epi = {'meanFunctional.nii'};  % EPI reference (target)
% src_t1  = {'UNI.nii'};             % skull-stripped T1 (source)
% segs    = {                       % tissue probability maps
%     'c1UNI.nii'
%     'c2UNI.nii'
%     'c3UNI.nii'
% };
% 
% %% ------------------------------------------------------------------------
% %  Step 1: Estimate coregistration (UNI -> meanFunctional)
% % -------------------------------------------------------------------------
% fprintf('Step 1: Coregister UNI → meanFunctional...\n');
% 
% clear matlabbatch
% matlabbatch{1}.spm.spatial.coreg.estimate.ref    = ref_epi;
% matlabbatch{1}.spm.spatial.coreg.estimate.source = src_t1;
% matlabbatch{1}.spm.spatial.coreg.estimate.other  = segs;
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
%     [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
% spm_jobman('run',matlabbatch);
% clear matlabbatch;
% 
% %% ------------------------------------------------------------------------
% %  Step 2: Reslice the segmentations into EPI space
% % -------------------------------------------------------------------------
% fprintf('Step 2: Reslice c1/c2/c3 to EPI grid...\n');
% 
% clear matlabbatch
% matlabbatch{1}.spm.spatial.coreg.write.ref = ref_epi;
% matlabbatch{1}.spm.spatial.coreg.write.source = segs;
% matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
% matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
% matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% spm_jobman('run',matlabbatch);
% clear matlabbatch;
% 
% %% ------------------------------------------------------------------------
% %  Step 3: Combine resliced segmentations into brain mask
% % -------------------------------------------------------------------------
% fprintf('Step 3: Combine rc1/rc2/rc3 to create brain mask...\n');
% 
% V1 = spm_vol('rc1UNI.nii');
% V2 = spm_vol('rc2UNI.nii');
% V3 = spm_vol('rc3UNI.nii');
% c1 = spm_read_vols(V1);
% c2 = spm_read_vols(V2);
% c3 = spm_read_vols(V3);
% 
% % liberal threshold (include GM+WM+CSF > 0.3)
% mask = (c1 + c2 + c3) > 0.3;
% 
% % fill small holes
% mask = imfill(mask,'holes');
% 
% Vmask = V1;
% Vmask.fname = 'rfBrainMask_S03.nii';
% Vmask.dt = [2 0];
% spm_write_vol(Vmask, mask);
% 
% fprintf('\n✅ Mask saved as rfBrainMask_S03.nii\n');
% 
% %% ------------------------------------------------------------------------
% %  Optional quick QC overlay
% % -------------------------------------------------------------------------
% fprintf('Launching SPM overlay for visual check...\n');
% spm_image('Display','meanFunctional.nii');
% spm_image('AddOverlay','rfBrainMask_S04.nii');
% fprintf('\n--- Done. Check overlay alignment visually ---\n');
