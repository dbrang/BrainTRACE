%% Pre-Process Tumor MRIs for Electrode Registration

% This script takes MRIs and (optionally) tumor masks to generate the
% necessary files used with BrainTRACE.
% T1 MRIs are skull stripped
% Other MRIs are registered and resliced to the T1
% A vasculature mask is generated from the Post-Gad or Flair image

% This code requires:
% SPM
% Freesurfer
% ANTs


clear
clc
close all

%%
% Pre-Process Tumor MRIs for Electrode Registration

% Steps
% 1. Export images from horos. 
%    Export to dicom file
%    Click 'Options' Folder Tree --> Flat
%
% 2. Rename folders one at a time as the following (rename as you export):
%    T1_NonGad
%    T1_PostGad
%    Flair
%    T2 (only export this image if no Flair)
% 
% 3. Update the SubDir

%% 

FS_Sub = '0534SF';
SubDir = ['~/Desktop/grid_registration/',FS_Sub,'/Electrodes/']; 

% Import MRIs from the renamed folders
cd(SubDir)

% Run through MRIs
mri_array = {'Flair','T1_NonGad','T1_PostGad','T2'};

% Check which folders exist
mri_array_temp = {};
currfiles = dir;
for mri_x=1:length(currfiles)
    temp = strcmp(currfiles(mri_x).name,mri_array);
    temp = find(temp);
    if ~isempty(temp)
        mri_array_temp = [mri_array_temp,mri_array(temp)];
    end
end

% Issue warning for missing ones
for mri_x=1:length(mri_array)
    temp = strcmp(mri_array{mri_x},mri_array_temp);
    temp = find(temp);
    if isempty(temp)
        warning(['Missing Folder: ',mri_array{mri_x}])
    end
end

% Update files to work on
mri_array = mri_array_temp;

%% Import files
for mri_x=1:length(mri_array)
    mri_curr = mri_array{mri_x};
    mri_files = dir([mri_curr,'/*.dcm']);
    
    clear output_text
    output_text(1,1) = {['matlabbatch{1}.spm.util.import.dicom.data = {']};
    
    % Get list of files in SubDir
    for file_x = 1:length(mri_files)
        output_text(1+file_x,1) = {[' ''',mri_files(file_x).folder,'/',mri_files(file_x).name,' ''']};
    end
    output_text(end+1,1) = {'};'};
    
    output_text(end+1,1) = {'matlabbatch{1}.spm.util.import.dicom.root = ''flat'';'};
    output_text(end+1,1) = {['matlabbatch{1}.spm.util.import.dicom.outdir = {''',pwd,'/',mri_curr,'''};']};
    output_text(end+1,1) = {'matlabbatch{1}.spm.util.import.dicom.protfilter = ''.*'';'};
    output_text(end+1,1) = {'matlabbatch{1}.spm.util.import.dicom.convopts.format = ''nii'';'};
    output_text(end+1,1) = {'matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;'};
    output_text(end+1,1) = {'matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;'};
    writecell(output_text,'SPM_Job.m','Delimiter','space','FileType','text','QuoteStrings',false)
    
    % List of open inputs
    nrun = 1; % enter the number of runs here
    jobfile = {'SPM_Job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});
    
    % Remove temp job script
    delete('SPM_Job.m')
    
    % Rename file and move up SubDir
    mri_curr_name = dir([mri_curr,'/*.nii']);
    mri_curr_name = mri_curr_name.name;
    movefile([mri_curr,'/',mri_curr_name],[pwd,'/',mri_curr,'.nii'])
end



%% Export images from Horos and Align images
FS_Dir = '/Applications/freesurfer/7.1.1/';

MRI_array.T1_NonGad = 'T1_NonGad.nii';
MRI_array.T1_PostGad = 'T1_PostGad.nii';
MRI_array.Flair = 'Flair.nii';
MRI_array.TumorMask = ''; 

%% Run Recon-All Stage 1

% If masks are separate, enter both
if iscell(MRI_array.TumorMask)
    input_list = '';
    for i=1:length(MRI_array.TumorMask)
        input_list = [input_list,' --i ',MRI_array.TumorMask{i},' '];
    end
   system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_concat ',input_list,...
    ' --o Combined_Mask.nii.gz --sum']);
    MRI_array.TumorMask = 'Combined_Mask.nii.gz';
    
    system(['export FREESURFER_HOME=/Applications/freesurfer;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_binarize --i Combined_Mask.nii.gz --min .1 --o Combined_Mask.nii.gz']);
end

% Run Recon-All (stage 1 only). Should take around 15 minutes.
% This will generate brainmask.mgz 
system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'recon-all -subjid ',FS_Sub,' -i ',MRI_array.T1_NonGad,' -autorecon1']);

% Run Recon-All Clinical. Should take around ~2 hours.
system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'recon-all-clinical.sh ',MRI_array.T1_NonGad,' ',FS_Sub,'_Clinical 4']);

% Copy T1.mgz and brainmask.mgz and clinical folder
system(['cp ',FS_Dir,'subjects/',FS_Sub,'/mri/T1.mgz ',SubDir]);
system(['cp ',FS_Dir,'subjects/',FS_Sub,'/mri/brainmask.mgz ',SubDir]);
system(['cp -r ',FS_Dir,'subjects/',FS_Sub,'_Clinical ',SubDir]);

% Convert T1.mgz to T1.nii.gz
system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_convert T1.mgz T1.nii']);

% Rename image
MRI_array.T1_NonGad = 'T1.nii';

% Convert Mask and MRIs if necessary
% Need to be .nii for SPM (not nii.gz)
MRI_array_names = fieldnames(MRI_array);
for i=1:length(MRI_array_names)
    mri_x = MRI_array.([MRI_array_names{i}]);
    if contains(mri_x,'.gz')
        % Convert .nii.gz to .nii
        system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
            'mri_convert ',mri_x,' ',mri_x(1:end-3)]);
        
        % Update variable
        MRI_array.([MRI_array_names{i}]) = mri_x(1:end-3);
    end
end

%% Use this code after recon-all to create RAS_Transforms.mat
% Check if files exist, if not, load mri and create
% This file will be loaded from the same path as the surface

mri = MRIread('T1.mgz');
RAS.tkrvox2ras = mri.tkrvox2ras;
RAS.vox2ras = mri.vox2ras;
save('RAS_Transforms.mat','RAS');

%% REGISTER AND RESLICE POSTGAD OR FLAIR TO T1.MGZ IN ORDER TO CREATE VASCULATURE

output_text = {};
Orig = MRI_array.T1_NonGad;

% if PostGad exists, use that 
if ~strcmp(MRI_array.T1_PostGad,'')
    Reg = MRI_array.T1_PostGad;
    %Other_Images = {MRI_array.Flair,MRI_array.TumorMask};
    
    Other_Images = {};
    if ~strcmp(MRI_array.Flair,'')
        Other_Images = [Other_Images,MRI_array.Flair];
    end
     if ~strcmp(MRI_array.TumorMask,'')
        Other_Images = [Other_Images,MRI_array.TumorMask];
     end
elseif ~strcmp(MRI_array.Flair,'')
    % Otherwise use flair
    Reg = MRI_array.Flair;
    Other_Images = {MRI_array.TumorMask};
else
    error('Need PostGad or Flair for registration. Neither found.')
end

% Remove empty cells
Other_Images(find(strcmp(Other_Images,'')))=[];

output_text(1,1) = {['matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {''',Orig,',1''};']};
output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.coreg.estwrite.source = {''',Reg,',1''};']};
if length(Other_Images)>0
    output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.other = {'};
    for i=1:length(Other_Images)
        output_text(end+1,1) ={[' ''',Other_Images{i},',1''']};
    end
    output_text(end+1,1) = {'};'};
end
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = ''nmi'';'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;'};
output_text(end+1,1) = {'matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = ''r'';'};
writecell(output_text,'SPM_Job.m','Delimiter','space','FileType','text','QuoteStrings',false)

% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'SPM_Job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

% Remove temp job script
delete('SPM_Job.m')

% Update file names
if ~strcmp(MRI_array.T1_PostGad,''),MRI_array.T1_PostGad = ['r',MRI_array.T1_PostGad];end
if ~strcmp(MRI_array.Flair,''),MRI_array.Flair = ['r',MRI_array.Flair];end
if ~strcmp(MRI_array.TumorMask,''),MRI_array.TumorMask = ['r',MRI_array.TumorMask];end
    
%% Process vasculature extraction

% Add in SubDir info for spm location
spm_dir = [fileparts(which('spm')),'/'];

% Iterate through MRIs
MRI_array_names = {};
if ~strcmp(MRI_array.T1_NonGad,'')
    MRI_array_names = [MRI_array_names,'T1_NonGad'];
end
if ~strcmp(MRI_array.T1_PostGad,'')
    MRI_array_names = [MRI_array_names,'T1_PostGad'];
end
if ~strcmp(MRI_array.Flair,'')
    MRI_array_names = [MRI_array_names,'Flair'];
end
%MRI_array_names = {'T1_NonGad','T1_PostGad','Flair'};
for i=1:length(MRI_array_names)
    mri_x = MRI_array.([MRI_array_names{i}]);
    
    % skip if empty
    if ~strcmp(mri_x,'')
        
        % Create data structure
        output_text = {};
        output_text(1,1) = {['matlabbatch{1}.spm.spatial.preproc.channel.vols = {''',mri_x,',1''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {''',spm_dir,'/tpm/TPM.nii,1''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {''',spm_dir,'/tpm/TPM.nii,2''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {''',spm_dir,'/tpm/TPM.nii,3''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {''',spm_dir,'/tpm/TPM.nii,4''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {''',spm_dir,'/tpm/TPM.nii,5''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {''',spm_dir,'/tpm/TPM.nii,6''};']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.affreg = ''mni'';']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;']};
        output_text(end+1,1) = {['matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];']};
        writecell(output_text,'SPM_Job.m','Delimiter','space','FileType','text','QuoteStrings',false)
        
        % List of open inputs
        nrun = 1; % enter the number of runs here
        jobfile = {'SPM_Job.m'};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(0, nrun);
        spm('defaults', 'PET');
        spm_jobman('run', jobs, inputs{:});
        
        % Remove temp job script
        delete('SPM_Job.m')
        
    end
end

%% Binarize the tumor mask
% Since resliced, need to double check it's binary
if ~strcmp(MRI_array.TumorMask,'')
    system(['export FREESURFER_HOME=/Applications/freesurfer;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_binarize --i ',MRI_array.TumorMask,' --min .1 --o ',MRI_array.TumorMask]);

    % Convert brain from ANTs to mgz
    system(['export FREESURFER_HOME=/Applications/freesurfer;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_convert ',MRI_array.TumorMask,' ',MRI_array.TumorMask(1:end-4),'.mgz']);
end

img1 = MRIread('T1.nii');

if ~strcmp(MRI_array.TumorMask,'')
    img2 = MRIread(MRI_array.TumorMask);
    
    % Check if same dims
    if img1.nvoxels==img1.nvoxels && isempty(find(round(img1.vox2ras,1) == round(img2.vox2ras,1)==0))
    else
        error('DIMs do not match. Check if using registered image.')
    end

end

%% Run ANTs for brain extraction

% ANTs expects NA for missing images
if strcmp(MRI_array.TumorMask,'')
    TumorMask = 'NA';
else
    TumorMask = MRI_array.TumorMask;
end

% setup path
if contains(mfilename,'LiveEditor')
    ScriptDir = [fileparts(which('Tumor_Grid_Preprocessing')),'/'];
    disp('Check script SubDir is correct')
else
    ScriptDir = [fileparts(which(mfilename)),'/'];
end
system(['/Library/Frameworks/R.framework/Resources/bin/Rscript ',ScriptDir,'/BrainExtract.r ',SubDir,' ',TumorMask,' T1.nii']);

%% Make 2 brain images

% Create brain.mgz since stop FS before this point
system(['export FREESURFER_HOME=/Applications/freesurfer/7.1.1/;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_mask -T .01 T1.mgz brainmask.mgz brain_FS.mgz']);

% Convert brain from ANTs to mgz
system(['export FREESURFER_HOME=/Applications/freesurfer/7.1.1/;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
    'mri_convert Syn__subBrainOnly.nii.gz brain_ANTs.mgz']);

% Copy header information from the T1.mgz to these files
img1 = MRIread('brain_FS.mgz');
img2 = MRIread('brain_ANTs.mgz');
img1.vol = img2.vol;
MRIwrite(img1,'brain_ANTs.mgz');


%% Create Vascularization masks; convert to .mgz

% Check if postgad exists
if ~strcmp(MRI_array.T1_PostGad,'')
    % Grow brainmask and apply to postgad
    system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
        'mri_binarize --i brainmask.mgz --min .00001 --dilate 1 --o brainmask_dilated.mgz']);
    system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
        'mri_mask ',MRI_array.T1_PostGad,' brainmask_dilated.mgz PostGad_BrainMask_Vasc.mgz']);
    
    % Use c3PostGad
    curr_MRI = ['c3',MRI_array.T1_PostGad];
    system(['export FREESURFER_HOME=/Applications/freesurfer/7.1.1/;source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
        'mri_binarize --i ',curr_MRI,' --min .00001 --o PostGad_c3_Mask.mgz']);
    system(['export FREESURFER_HOME=',FS_Dir,';source $FREESURFER_HOME/SetUpFreeSurfer.sh;',...
        'mri_mask ',MRI_array.T1_PostGad,' PostGad_c3_Mask.mgz PostGad_c3_Vasc.mgz']);
end

% Can also edit these output images if necessary
% Sometimes dura will show brightly and needs to be deleted





