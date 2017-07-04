function [S options] = cdr_loadImages(source, options, S_all)  
% 

% Loads the raw images into a STACK of images, or accepts a stack of images
% passed as an argument. Performs resizing, sorting, and compressing of the
% image data to keep the model generation fast.
% 
% Usage:  [S options] = cdr_loadImages(source) 
%
% Input: SOURCE  can be a path to a folder containing images, a path with a
%                filter (eg "/images/*.tif"), or an RxCxN array containing 
%                the images to be corrected (R = height, C = width, N = 
%                number of images)
%
% Output: S      an array containing the image data to be used by CIDRE
%         
%         OPTIONS returns the options structure with appropriatetly 
%                 modified parameters: the size of the original images
%                 [height width], the number of input images provided, the
%                 working image size, the source folder, the filenames of 
%                 the source images
%
% See also: cidre, cidreGui

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

% S_all ;   %#ok<NASGU> the array containing the image data returned to CIDRE
% maxI        = 0;    % the max intensity found in the supplied images




%% obtain an array S, containing the source image data at the working image
%% size. S can be loaded from a source path (+filter) or from an array 
%% passed as an argument
    % source is a char defining the path (+filter) of files to open

    
    % break source into a path, filter, and extension
    [pth filter ext] = fileparts(source);
    
    
    % check if a file filter is provided, if so use it to generate a list
    % of source filenames
    mosaicFile=getTiledAcquisitionParamFile;
    param=readMetaData2Stitchit(mosaicFile);
%     numTiles = param.numTiles.X*param.numTiles.Y;
%     numOptsec = param.mosaic.numOpticalPlanes;
%     indsec = 0:numOptsec-1;
%     start = indsec.*numTiles;
    d = options.ALLfilenames;
    

%     ind = 1:options.ALLnum_images_provided;
%         else
%             tile_start = start(optical_section):numOptsec*numTiles:options.ALLnum_images_provided-1; % number of data sets
%             tile_array = repmat(tile_start,numTiles,1);
%             array_plus = repmat(0:numTiles-1,size(tile_array,2),1);
%             tile_num = tile_array+array_plus';
%             tile_num = reshape(tile_num,1,size(tile_num,1)*size(tile_num,2));
%             opt_section_load = zeros(1,options.ALLnum_images_provided);
%             for sec_num = 1:options.ALLnum_images_provided
%                 str = d{sec_num};
%                 index_chanel = strfind(str, [filter(2:end) ext]);
%                 index_tire = strfind(str, '-');
%                 opt_section_load(sec_num) = abs(sscanf(str(index_tire(end):index_chanel), '%i', str2num(filter(end))));
%             end
%             [ind val] = ismember(opt_section_load,tile_num);
%         end
%      d = d(ind);
            % save filenames of the same optical sections or all the images
            % in case of no separation of the sections
%             options.filenames = cell(numel(d),1);
%             for i = 1:numel(d)
%                 options.filenames{i} = d{i};
%             end
        
    % check if we have found any valid files in the folder
%     if numel(options.filenames) == 0
%         error('CIDRE:loadImages', 'No files found.'); 
%     end
    
    % store the number of source images into the options structure
%     options.num_images_provided = numel(options.filenames);
   
%      S = S_all(:,:,ind);
     maxI = max(max(max(S_all)));
%     fprintf('finished in %1.2fs.\n', toc(t1));

%% apply several processing steps to the image stack S
% Now that we have loaded the stack as an RxCxN array (where R*C ~=
% options.target_num_pixels), we must do check if there is sufficient
% intensity information in the stack, sort the intensity values at each
% (x,y) image location, and compress the stack in the 3rd dimension to keep
% the computation time manageable
% 
% options = get_bit_depth(options, maxI); % store the bit depth of the images in options, needed for entropy measurement
% entropy = get_entropy(S, options);      % compute the stack's entropy
% 
% N       = size(S,3);    % number of images in the stack
% a       = 7.838e+06;    % parameters of a fitted exponential function
% b       = -1.948;       % parameters of a fitted exponential function
% c       = 20;           % parameters of a fitted exponential function
% N_required = a*exp(b*entropy) + c; 
% if N < N_required
% %     S = S_all;
% %     maxI = max(max(max(S)));
%     fprintf('There is enough contrast in the channel\n')
% end

[S options] = cdr_preprocessData(S_all, maxI, options);










% 
% 
% 
% 
% function [R_working C_working] = determine_working_size(image_size, N_desired)
% % determines a working image size based on the original image size and
% % the desired number of pixels in the working image, N_desired
% 
% 
% R_original = image_size(1);
% C_original = image_size(2);
% 
% scale_working = sqrt( N_desired/(R_original*C_original));
% 
% R_working = round(R_original * scale_working);
% C_working = round(C_original * scale_working);
% 
% 
% 
% function options = get_bit_depth(options, maxI)
% % Sets options.bit_depth describing the provided images as 8-bit, 12-bit, or 
% % 16-bit. If options.bit_depth is provided, it is used. Otherwise the bit 
% % depth is estimated from the max observed intensity, maxI.
% 
% if ~isempty(options.bit_depth)
%     if ~ismember(options.bit_depth, [2^8 2^12 2^16])
%         error('CIDRE:loadImages', 'Provide bit depth as max integer value, eg 2^12');
%     else
%         fprintf(' %d-bit depth\n', log2(options.bit_depth));
%     end
% else    
%     if maxI > 2^12
%         options.bit_depth = 2^16;
%     elseif maxI > 2^8
%         options.bit_depth = 2^12;
%     else
%         options.bit_depth = 2^8;
%     end    
%     fprintf(' %d-bit depth images (estimated from max intensity=%1.0f)\n', log2(options.bit_depth), maxI);
% end


% 
% function H = get_entropy(S, options)
% % gets the entropy of an image stack. A very low entropy indicates that
% % there may be insufficient intensity information to build a good model.
% % This can happen when only a few images are provided and the background
% % does not provide helpful information. For example in low confluency
% % fluorescence images from a glass slide, the background pixels have nearly
% % zero contribution from incident light and do not provide useful
% % information.
% 
% 
% % set one bin for every potential intensity level
% bins = 1:options.bit_depth;
% 
% % get a distrubtion representing all of S
% P = hist(S(:),bins);
% P = P/sum(P);
% 
% % compute the entropy of the distribution
% if sum(~isfinite(P(:)))
%    error('CIDRE:loadImages', 'the inputs contain non-finite values!') 
% end
% P = P(:) ./ sum(P(:));
% P(P == 0) = []; % In the case of p(xi) = 0 for some i, the value of the 
%                 % corresponding sum and 0 logb(0) is taken to be 0
% temp = P .* log2(P);
% H = -sum(temp);
% fprintf(' Entropy of the stack = %1.2f\n', H);








