function NotEnough = check_entropy(S_all, options)       

[pth filter ext] = fileparts(options.source);

optical_section = 1;
% check if a file filter is provided, if so use it to generate a list
% of source filenames
mosaicFile=getTiledAcquisitionParamFile;
param=readMetaData2Stitchit(mosaicFile);
numTiles = param.numTiles.X*param.numTiles.Y;
numOptsec = param.mosaic.numOpticalPlanes;
indsec = 0:numOptsec-1;
start = indsec.*numTiles;
d = options.ALLfilenames;

tile_start = start(optical_section):numOptsec*numTiles:options.ALLnum_images_provided-1; % number of data sets
tile_array = repmat(tile_start,numTiles,1);
array_plus = repmat(0:numTiles-1,size(tile_array,2),1);
tile_num = tile_array+array_plus';
tile_num = reshape(tile_num,1,size(tile_num,1)*size(tile_num,2));
opt_section_load = zeros(1,options.ALLnum_images_provided);
for sec_num = 1:options.ALLnum_images_provided
    str = d{sec_num};
    index_chanel = strfind(str, [filter(2:end) ext]);
    index_tire = strfind(str, '-');
    opt_section_load(sec_num) = abs(sscanf(str(index_tire(end):index_chanel), '%i', str2num(filter(end))));
end
[ind val] = ismember(opt_section_load,tile_num);
d = d(ind);

options.filenames = cell(numel(d),1);
for i = 1:numel(d)
    options.filenames{i} = d{i};
end
S = S_all(:,:,ind);

maxI = max(max(max(S)));
options = get_bit_depth(options, maxI); % store the bit depth of the images in options, needed for entropy measurement
entropy = get_entropy(S, options);      % compute the stack's entropy

N       = size(S,3);    % number of images in the stack
a       = 7.838e+06;    % parameters of a fitted exponential function
b       = -1.948;       % parameters of a fitted exponential function
c       = 20;           % parameters of a fitted exponential function
N_required = a*exp(b*entropy) + c; 
NotEnough = 0;
if N < N_required
    fprintf('NOT enough images...\n The optical sections are NOT separated.\n')
    NotEnough = 1;
end


function options = get_bit_depth(options, maxI)
% Sets options.bit_depth describing the provided images as 8-bit, 12-bit, or 
% 16-bit. If options.bit_depth is provided, it is used. Otherwise the bit 
% depth is estimated from the max observed intensity, maxI.

if ~isempty(options.bit_depth)
    if ~ismember(options.bit_depth, [2^8 2^12 2^16])
        error('CIDRE:loadImages', 'Provide bit depth as max integer value, eg 2^12');
    else
        fprintf(' %d-bit depth\n', log2(options.bit_depth));
    end
else    
    if maxI > 2^12
        options.bit_depth = 2^16;
    elseif maxI > 2^8
        options.bit_depth = 2^12;
    else
        options.bit_depth = 2^8;
    end    
    fprintf(' %d-bit depth images (estimated from max intensity=%1.0f)\n', log2(options.bit_depth), maxI);
end



function H = get_entropy(S, options)
% gets the entropy of an image stack. A very low entropy indicates that
% there may be insufficient intensity information to build a good model.
% This can happen when only a few images are provided and the background
% does not provide helpful information. For example in low confluency
% fluorescence images from a glass slide, the background pixels have nearly
% zero contribution from incident light and do not provide useful
% information.


% set one bin for every potential intensity level
bins = 1:options.bit_depth;

% get a distrubtion representing all of S
P = hist(S(:),bins);
P = P/sum(P);

% compute the entropy of the distribution
if sum(~isfinite(P(:)))
   error('CIDRE:loadImages', 'the inputs contain non-finite values!') 
end
P = P(:) ./ sum(P(:));
P(P == 0) = []; % In the case of p(xi) = 0 for some i, the value of the 
                % corresponding sum and 0 logb(0) is taken to be 0
temp = P .* log2(P);
H = -sum(temp);
fprintf(' Entropy of one optical section stack = %1.2f\n', H);


