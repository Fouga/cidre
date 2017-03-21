function imCorrected = correct_cidre(im,model)

z  = model.z;
v = model.v;
imCorrected = ((double(im) - repmat(z,1,1,size(im,3)))./repmat(v,1,1,size(im,3)));
imCorrected = uint16(imCorrected);

