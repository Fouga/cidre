function model = cidre(source, varargin)
% Modified in order to load images from Tissue Vision microscopy.
% The shading image is calculated for each channel separately.
%----------------------------------------------------------

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
%
% The function was modified by Natalia Chicherova, Basel 2017

if nargin==0	
    help(mfilename)
	return
end

% add necessary paths
if (~isdeployed)
    p = mfilename('fullpath');
    p = p(1:end-5);
    addpath(fullfile(p, '3rdparty')); 
%     addpath([p '/gui/']);
    addpath(fullfile(p, 'io'));
    addpath(fullfile(p, 'main'));
end


% parse the input arguments, return a structure containing parameters
options = cdr_parseInputs(varargin);
options.cidrePath = p;

% load all the images 
% break source into a path, filter, and extension
[pth filter ext] = fileparts(source);
pth = fullfile(pth, '**');

% store the source path in the options structure
options.folder_source = pth;  
options.source = source; 
% read all the files in d
d = rdir(fullfile(options.folder_source, [filter ext]));
% create a new field for all file names of the channel in options.
options.ALLfilenames = cell(numel(d),1);

for i = 1:numel(d)
    options.ALLfilenames{i} = d(i).name;
end
% store the number of source images into the options structure
options.ALLnum_images_provided = numel(options.ALLfilenames);

% load one image to check the downsampling size
I = imread( options.ALLfilenames{1});
if numel(size(I)) == 3
    error('CIDRE:loadImages', 'Non-monochromatic image provided. CIDRE is designed for monochromatic images. Store each channel as a separate image and re-run CIDRE.'); 
end
options.image_size = size(I);
[R C] = determine_working_size(options.image_size, options.target_num_pixels);
options.working_size = [R C];
clear I

% preallocate memory for loaded images
S_all = zeros([options.working_size options.ALLnum_images_provided]);
NAMES = options.ALLfilenames;
fprintf(' Reading %d images from %s, chanel %s\n', ...
    options.ALLnum_images_provided, options.folder_source(1:end-3), filter(end-1:end));

% load all the images from the given channel into the memory
parfor z = 1:options.ALLnum_images_provided
       I = imread(NAMES{z});
       S_all(:,:,z) = double(imresize(I, [R C]));
end


% same model for all the optical sections
maxI = max(max(max(S_all)));
[S options] = cdr_preprocessData(S_all, maxI, options);
[model, options] = cdr_cidreModel(S,options);

% save the correction model to the destination folder
channel = str2double(filter(end));
% filename = sprintf('%s%s', options.folder_destination, ...
%     sprintf('/cidre_channel_%i.mat',channel));
filename = fullfile(options.folder_destination,sprintf('cidre_channel_%i.mat',channel));
save( filename, 'model');    

    
    
    

function [R_working C_working] = determine_working_size(image_size, N_desired)
% determines a working image size based on the original image size and
% the desired number of pixels in the working image, N_desired


R_original = image_size(1);
C_original = image_size(2);

scale_working = sqrt( N_desired/(R_original*C_original));

R_working = round(R_original * scale_working);
C_working = round(C_original * scale_working);
