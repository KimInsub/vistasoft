function [scan,hdr]=Dicom2ThreeDAnalyze_SPM5_NIC_viaMPRAGE(scan,baseDir,dicomServer)
% scan=Dicom2ThreeDAnalyze_SPM5_NIC(Scan,baseDir,dicomServer);
% where scan has to be a struct with the fields
% Scan.DirName
% Scan.Slices
% Scan.Volumes
% Scan.Filenames
% Scan.Filename4d
%
% dicomServer indicates where the data have come from
% 1: NIC vis Hurricane
% 2: NIC via MPRAGE
% 3: China Basin via Hurricane
%
% the script will generate a 3danalyze file
% and put the name into the field Scan.Filename3d
% written 2006.01.12 by Mark Schira mark@ski.org
% Modified to use the SPM5 dicom maker rather than medcon. We still use fsl
% for the 4d file generation.
% The critical step is...
%  % convert dicoms to analyzer
%        hdr = spm_dicom_headers(files);
%        spm_dicom_convert(hdr,opts,root_dir,format);
% Rather than using medcon (great though medcon is!)
%  Possible Options and Values are
%
%       opts - options
%              'all'      - all DICOM files (default)
%              'mosaic'   - the mosaic images
%              'standard' - standard DICOM files
%              'spect'    - SIEMENS Spectroscopy DICOMs (position only)
%                           This will write out a mask image volume with 1's
%                           set at the position of spectroscopy voxel(s).
%              'raw'      - convert raw FIDs (not implemented)
%       root_dir - 'flat' - SPM5 standard, do not produce file tree
%                  With all other options, files will be sorted into
%                  directories according to their sequence/protocol names
%                  'date_time' - Place files under ./<StudyDate-StudyTime>
%                  'patid'         - Place files under ./<PatID>
%                  'patid_date'    - Place files under ./<PatID-StudyDate>
%                  'patname'          - Place files under ./<PatName>
%       format - output format
%                'img' Two file (hdr+img) NIfTI format (default)
%                'nii' Single file NIfTI format
%                      All images will contain a single 3D dataset, 4D images
%                      will not be created.
%% EDIT THESE VARIABLES
opts     = 'all';
root_dir = 'patid';
format   = 'img'; % These are 2-file nifti's rather than analyze. Beware!


if ~exist('baseDir','var')
    baseDir = pwd;
end

Directory=[baseDir,filesep,'RawDicom',filesep,scan.DirName,filesep];
Filename=scan.Filenames;
Slices=scan.Slices;
%Volumes=scan.Volumes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Stuff for fsl
% fslBase='/raid/MRI/toolbox/FSL/fsl';
% if (ispref('VISTA','fslBase'))
%     disp('Setting fslBase to the one specified in the VISTA matlab preferences:');
%     fslBase=getpref('VISTA','fslBase');
%     disp(fslBase);
% end
% 
% fslPath=fullfile(fslBase,'bin'); % This is where FSL lives - should also be able to get this from a Matlab pref
% 

 
%******************************************************************
%now do the Recon




switch dicomServer
    case 1 % Parnassus NIC Via MPRAGE
        % Do a 'dir' list of everything inside the directory 
        fileNameList=dir(Directory);
        fileNameList=fileNameList(3:end);
        
%         
%         fileBase=[Directory,scan.Filenames,'M-'];
% 
%         % Generate a list of filenames to pass into the spm world
% 
       for thisVol=1:Slices
           fileName{thisVol}=[Directory,filesep,fileNameList(thisVol).name];
         end 
case 1 % Parnassus NIC Via MPRAGE
%         fileBase=[Directory,scan.Filenames,'M-'];
% 
%         % Generate a list of filenames to pass into the spm world
% 
%         for thisVol=1:Slices;
% 
%             fileName{thisVol}=sprintf('%s%.4d-%.4d.dcm',fileBase,str2num(scan.DirName),thisVol);
%         end

    case 2 % Parnassus NIC Via Hurricane
        fileBase=[Directory,scan.Filenames];

        % Generate a list of filenames to pass into the spm world

        for thisVol=1:Slices;
            fileName{thisVol}=sprintf('%s%d.DCM',fileBase,thisVol)
        end
end


% Generate a list of filenames to pass into the spm world


disp('Generating headers');
hdr = spm_dicom_headers(strvcat(fileName));
disp('Converting dicoms');
spm_dicom_convert(hdr,opts,root_dir,format);
disp('Converted');
% Files come out of spm_convert in a directory composed of the
% patientID (here listed as SKERI' + the SeriesDescription (for
% example ep2d_pace_dynt_moco_128) + the scan number eg 0014

outBaseDir=sprintf([baseDir,filesep,strtrim(hdr{1}.PatientID),filesep,hdr{1}.SeriesDescription(~isspace(hdr{1}.SeriesDescription)),'_%0.4d'],hdr{1}.SeriesNumber);
 

if (strcmp(format,'img')) % using 2-part niftis.
    existFiles=dir([outBaseDir,filesep,'*.hdr']); % This does list in ascending order
    nFiles=length(existFiles);
	 if nFiles==0
		 % try dropping periods & hyphens
		 outBaseDir = strrep(outBaseDir,'.','');
		 [tmpBase,tmpDir] = fileparts(outBaseDir);
		 outBaseDir = fullfile(tmpBase,strrep(tmpDir,'-',''));
		 existFiles=dir([outBaseDir,filesep,'*.hdr']);
		 nFiles = length(existFiles);
	 end
    disp(nFiles);
    if (nFiles~=scan.Volumes)
        warning('Problem detected. Wrong number of img files compared to volumes. Taking only the first volume..');
	 end

	 outBaseDir = strrep(outBaseDir,' ','\ ');	% make sure this is after dir([outBaseDir...]) 
    FilenameList=[''];
    for thisFileName=1:1 % was nFiles
        FilenameList=[FilenameList,[' ',outBaseDir,filesep,existFiles(thisFileName).name]];
    end


    % Finally, call FSL to generate 3d analyze files.
%    OutVolStr=sprintf('%s_%.3d_3d.hdr',strtrim(hdr{1}.PatientID),hdr{1}.SeriesNumber);
   OutVolStr=sprintf('%s_%.3d_3d',strtrim(hdr{1}.PatientID),hdr{1}.SeriesNumber);
	OutVolPath=[baseDir,filesep,'Reconned',filesep,OutVolStr];
	if isunix		%  in case there're spaces in the path string
		OutVolPathSys=[strrep(baseDir,' ','\ '),filesep,'Reconned',filesep,OutVolStr];
	else
		OutVolPathSys = OutVolPath;
	end

end

cmd=['cp ',outBaseDir,filesep,'*.hdr ',OutVolPathSys,'.hdr'];
disp(cmd);
system(cmd);
cmd=['cp ',outBaseDir,filesep,'*.img ',OutVolPathSys,'.img'];
disp(cmd);
system(cmd);

scan.Filename3d=[OutVolPath];
display(['saved ',OutVolPath]);




