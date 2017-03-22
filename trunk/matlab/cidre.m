function model = cidre(source, varargin)
% Illumination correction for optical microscopy. Apply to a collection of 
% monochromatic images. Multi-channel images should be separated, and each 
% channel corrected separately.
% 
% Usage:  MODEL = cidre(source, ...)
%
% Input:  SOURCE can be a path to a folder containing images, a path with a
%                filter (eg "/images/*.tif"), or an RxCxN array containing 
%                the images to be corrected (R = height, C = width, N = 
%                number of images).
%
%                CIDRE also supports the following optional arguments,
%                which override default parameters.
%       
%         'destination' specifies a folder where corrected images and 
%                correction model are stored. If empty, CIDRE will return a
%                model but not store corrected images.
%       
%         'correction_mode' (default = 2) specifies the type of correction
%                to perform:  0 "Zero-light preserved" retains the original 
%                intensity range and zero-light level of the original 
%                images. 1 "Dynamic range corrected" retains the intensity 
%                range of the original images. 2 "Direct" subtracts the 
%                zero-light term and divides the illumination gain. 
%
%         'lambda_v' (default = 6.0) high values (eg. 9.5) increase the
%                spatial regularization strength, yielding a more smooth v
%                intensity gain surface.
%
%         'lambda_z' (default = 0.5) high values (eg. 3) increase the
%                zero-light regularization strength, enforcing a more
%                uniform z surface.
%
%         'z_limits' (default none) user may specify known limits for the z
%                z surface, using dark frame images for example. Specify
%                the limits as a 2-element array (e.g. 'z_limits',[98 105])
%
%         'bit_depth' specifies the bit depth of the images to be correct.
%                The bit depth is automatically detected, only use this
%                option if the detection fails.
%
%          'q_percent' (default = 0.25) specifies the proportion of the 
%                data used to compute the robust mean, Q.
%
%          'max_lbfgs_iterations' (default = 500) specifies the maximum
%                number of iterations allowed in the optimization.
%
% Output: MODEL  a structure containing the correction model used by CIDRE
%                to correct the source images
%        
%
% See also: cidreGui

% From the CIDRE project, an illumination correction method for optical
% microscopy (https://github.com/smithk/cidre).
% Copyright © 2015 Kevin Smith and Peter Horvath. Scientific Center for 
% Optical and Electron Microscopy (SCOPEM), Swiss Federal Institute of 
% Technology Zurich (ETH Zurich), Switzerland. All rights reserved.
%
% CIDRE is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License version 2 (or higher) 
% as published by the Free Software Foundation. See the license file in
% the root folder. This program is distributed WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% This software includes a folder "3rdparty" containing minFunc, a 3rd
% party software implementing L-BFGS. MinFunc is licensed under the
% Creative Commons, Attribute, Non-Commercial license. To use this software 
% for commercial purposes, you must replace minFunc with other software. 
% Matlab offers an alternative (slower) implementation in the function 
% fminlbfgs.

if nargin==0
	help(mfilename)
	return
end

% add necessary paths
if (~isdeployed)
    p = mfilename('fullpath');
    p = p(1:end-5);
    addpath([p '/3rdparty/']); 
    addpath([p '/gui/']);
    addpath([p '/io/']);
    addpath([p '/main/']);
end


% load names of inside folders
% userConfig=readStitchItINI;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%Find the raw directories we will descend into. 
% paramFile=getTiledAcquisitionParamFile;
% param=readMetaData2Stitchit(paramFile); 


% parse the input arguments, return a structure containing parameters
options = cdr_parseInputs(varargin);
options.cidrePath = p;

cdr_gui_toggle(options)

% load all the images 
% break source into a path, filter, and extension
% source = fullfile(userConfig.subdir.rawDataDir,source);
[pth filter ext] = fileparts(source);
pth = [pth '/**/'];

% store the source path in the options structure
options.folder_source = pth;  
options.source = source; 
% read all the files in d
d = rdir([options.folder_source filter ext]);
% create a new field for all file names of the channel in options.
options.ALLfilenames = cell(numel(d),1);

for i = 1:numel(d)
    options.ALLfilenames{i} = d(i).name;
end
% store the number of source images into the options structure
options.ALLnum_images_provided = numel(options.ALLfilenames);
    
    % check if results already exist. If at least one image does not exist,
    % do the calculation
%     count = 0;
%     for z=1:options.ALLnum_images_provided
%         FileName = options.ALLfilenames{z};
%         fold_name = FileName(length(options.folder_source)-2:end);
%         filename = sprintf('%s%s', options.folder_destination , fold_name);
%         if ~exist(filename) || options.overwrite
%             count = count+1;
%         else 
%             if mod(z,options.ALLnum_images_provided) == 0; fprintf('Some images are already corrected with CIDRE. \n'); end
%         end
%     end
    

% load all the images from the given channel into the memory
I = imread( options.ALLfilenames{1});
if numel(size(I)) == 3
    error('CIDRE:loadImages', 'Non-monochromatic image provided. CIDRE is designed for monochromatic images. Store each channel as a separate image and re-run CIDRE.'); 
end
options.image_size = size(I);
[R C] = determine_working_size(options.image_size, options.target_num_pixels);
options.working_size = [R C];

% preallocate memory for loaded images
S_all = zeros([options.working_size options.ALLnum_images_provided]);
NAMES = options.ALLfilenames;
fprintf(' Reading %d images from %s, chanel %s\n', options.ALLnum_images_provided, options.folder_source(1:end-3), filter(end-1:end));
parfor z = 1:options.ALLnum_images_provided
%             if mod(z,100) == 0; fprintf('.'); end  % progress to the command line
    try
        I = imread(NAMES{z});
    catch 
        disp 'The image is broken. Coping its neighbour.';
        % if an image is broken but exist
        SplitedName1 = strsplit (NAMES{z},'-');
        SplitedName2 = SplitedName1{4};
        SplitedName = strsplit (SplitedName2,'_');
        nUpperLayer = str2double(SplitedName(:,1))+nImagesPerLayer
        nLowerLayer = str2double(SplitedName(:,1))-nImagesPerLayer
        UpperLayerImage = strcat (SplitedName1(:,1),'-',SplitedName1(:,2),'-', SplitedName1(:,3),'-',int2str(nUpperLayer),'_', SplitedName(:,2))
        LowerLayerImage = strcat (SplitedName1(:,1),'-',SplitedName1(:,2),'-', SplitedName1(:,3),'-',int2str(nLowerLayer),'_', SplitedName(:,2))
        if exist (UpperLayerImage{1}, 'file')
               copyfile (UpperLayerImage{1}, NAMES{z});
        elseif exist (LowerLayerImage{1}, 'file')
               copyfile (LowerLayerImage{1}, NAMES{z});
        else 
               disp 'You have no image to fill the hole, haha';% copy black tile!!
        end
        I = imread(NAMES{z});
    end
%     I = double(I);
%         maxI = max(maxI, max(I(:)));
    Irescaled = imresize(I, [R C]);
    S_all(:,:,z) = double(Irescaled);

end


% check if there is enough contrast
  
% NotEnough = check_entropy(S_all, options);
% 
% if ~NotEnough
%     % separate optical sections to find background models
%     % read from Mosaic file number of layers
%     mosaicFile=getTiledAcquisitionParamFile;
%     param=readMetaData2Stitchit(mosaicFile);
%     numOptsec = param.mosaic.numOpticalPlanes;
%     for optical_section = 1:numOptsec
% 
%         % load the data, either from a folder or from a passed array
%         [S, options] = cdr_loadImages(source, options,optical_section, S_all);
% 
%         [model, options] = cdr_cidreModel(S,options);
%         filename = sprintf('%s%s', options.folder_destination, chanel_folder, '/', ...
%             sprintf('cidre_model_optical_section_%i.mat',optical_section));
%         save(filename, 'model', 'options');
%         % correct the source images, if requested
% %                 cdr_correct(model,options, optical_section);
% 
%         cdr_gui_toggle(options)
%     end
% else

    % same model for all the optical sections
    optical_section = 0;
    [S, options] = cdr_loadImages(source, options,optical_section, S_all);
    [model, options] = cdr_cidreModel(S,options);

    % save the correction model to the destination folder
    chanel = str2double(filter(end));
    filename = sprintf('%s%s', options.folder_destination, ...
        sprintf('/cidre_chanel%i_optical_section_%i.mat',chanel,optical_section));
    save( filename, 'model');
%              cdr_gui_toggle(options)
% end
    
    
    

function [R_working C_working] = determine_working_size(image_size, N_desired)
% determines a working image size based on the original image size and
% the desired number of pixels in the working image, N_desired


R_original = image_size(1);
C_original = image_size(2);

scale_working = sqrt( N_desired/(R_original*C_original));

R_working = round(R_original * scale_working);
C_working = round(C_original * scale_working);
