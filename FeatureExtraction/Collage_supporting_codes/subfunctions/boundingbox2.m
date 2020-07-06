function varargout = boundingbox2(vol,mask,n,slice_opt,disp_opt)
% Create a cropped volume (& mask) using a bounding box of n-pixels outside the mask boundary

% example call:
% [cv,cm] = boundingbox2(img,mask)

%{
INPUTS
vol = original volume/image
mask = corresponding mask
n = number of pixels away from outer edges of mask (width of bounding box) in x-y plane
slice_opt (optional) use all slices or crop  or in z-direction as well ('allz' [default]) or 'cropz')
disp_opt = (optional) display output or not ('on' [default] or 'off')
%}

%{
OUTPUTS
varargout{1} = croppedV (the cropped version of vol)
varargout{2} = croppedM (the cropped version of masks)
%}

%{
 ____________________________________________
|(1,1)          ...           (1,size(vol,2)|
|                                           |
|                                           |
|                                           |
|                                           |
|(size(vol,1),1)...(size(vol,1),size(vol,2))|
_____________________________________________

%}

% @Jacob Antunes, 2015
% Last Updated: 04-23-2018

%% check inputs
if nargin < 4
    slice_opt = 'allz';
end

if ~exist('slice_opt','var') || (~strcmp(slice_opt,'allz') && ~strcmp(slice_opt,'cropz'))
        error('Invalid or missing SLICE_OPT option.');
end

if nargin < 5
    disp_opt = 'on';
end

if ~exist('disp_opt','var') || (~strcmp(disp_opt,'on') && ~strcmp(disp_opt,'off'))
        error('Invalid or missing DISP_OPT option.');
end

%% make all zero values in the original volume nonzero
volcopy = vol;
minval = max(vol(:))/1000;
vol(find(vol == 0)) = minval;

if ~isempty(mask) %use the mask to create a bounding box
    [row, col, slice] = ind2sub(size(mask),find(mask ~= 0)); 
else
    [row, col, slice] = ind2sub(size(vol),find(vol ~= minval));
end

%% find dimensions of bounding box
bottom = min(max(row)+n,size(vol,1));
top = max(min(row)-n,1);
left = max(min(col)-n,1);
right = min(max(col)+n,size(vol,2));
if strcmp(slice_opt,'allz')
    front = 1;
    back = size(vol,3);
elseif strcmp(slice_opt,'cropz')
    front = max(min(slice)-n,1);
    back = min(max(slice)+n,size(vol,3));
end

%% apply bounding box
if ~strcmp(disp_opt,'off')
    fprintf('Cropping volume and mask with padding of %i pixels\n',n);
end
croppedV = volcopy(top:bottom,left:right,front:back);
croppedM = mask(top:bottom,left:right,front:back);


varargout{1} = croppedV;
if nargout == 2
    varargout{2} = croppedM;
end

%testing purposes
% subplot(2,2,1);
% imshow(mask(:,:,20));
% subplot(2,2,3);
% imshow(vol(:,:,20));
% subplot(2,2,4);
% imshow(croppedV(:,:,20));

%view results
% v(:,:,:,1) = vol;
% v(:,:,:,2) = croppedV;
% mosaic4D(v);

end


