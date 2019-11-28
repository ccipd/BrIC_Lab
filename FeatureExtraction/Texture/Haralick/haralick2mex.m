%HARALICK2MEX MEX-ed Haralick feature extraction code in 2D
%   features = haralick2mex(IMAGE,GRAYLEVELS,WS,DIST,BACKGROUND)
%   
%   IMAGE is a 2D input image
%
%   GRAYLEVELS is how many graylevels you are considering for co-occurrence
%   NOTE: You should "rescale" your IMAGE to utilize the range of
%   graylevels you are utilizing
%   e.g. image = round(rescale_range(image,0,127));
%   In this example, GRAYLEVELS = 128.
%
%   WS is the window-size within which statistics are calculated
%   e.g. WS = 5 will use a 5 x 5 window, i.e. 2 pixels in all directions
%   around the "center pixel"
%
%   DIST is how many co-occurring neighbors you want to utilize around 
%   each pixel within the window (defined by WS), for determining the
%   co-occurence matrix
%   e.g. DIST = 1 will use a 3 x 3 block around each pixel inside the
%   window to calculate co-occurrences
%   NOTE: You can utilize different values for WS and DIST.
%
%   BACKGROUND is the value of "background" intensities which will be
%   ignored during all calculations.
%
%   FEATURES will be a X x Y x 13 matrix, where X,Y are the input
%   dimensions of IMAGE
%
% Notes: 
% - haralick2mexmt is the multi-threaded version of this code
% - haralick3mex,haralick3mexmt does this in 3D, i.e. a N x N x N window is used
% - haralick3smex,haralick3smexmt does this in 3D, single precision
% - Names of features: {'entropy','energy','inertia','idm','correlation','info1','info2','sum_av','sum_var','sum_ent','diff_av','diff_var','diff_ent'}
%
% See also: graycomatrix.m, graycoprops.m
%
% Example
% I = double(imread('circuit.tif'));
% datarange=max(I(:))-min(I(:));
% N1 = 0; N2 = 127; % rescale between 0-127 (128 graylevels)
% wantedrange=N2-N1; 
% I = round(N1+(I-min(I(:)))/(datarange/wantedrange));
% Features = haralick2mex(I,128,5,1,0);
