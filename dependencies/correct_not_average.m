function imCorrected = correct_not_average(im,model,varargin)

if nargin==2
    method = 'all_corrections';
else
    method = varargin{1};
end
switch method
    case 'all_corrections'
        z  = model.z;
        v = model.v;
        imCorrected = ((double(im) - repmat(z,1,1,size(im,3)))./repmat(v,1,1,size(im,3)));
        imCorrected = uint16(imCorrected);
    case 'flatfield'
        v = model.v;
        imCorrected = (double(im)./repmat(v,1,1,size(im,3)));
        imCorrected = uint16(imCorrected);
end