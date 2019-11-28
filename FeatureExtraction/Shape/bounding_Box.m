%% boundingBox
%Jacob Antunes
%20151009

%Create a cropped volume mask using a bounding box of n-pixels outside the
%annotated volume mask

function [croppedV,croppedM] = boundingBox(vol,mask,n)

%{
INPUTS
vol = original 3D volume
mask = original 3D mask
n = number of pixels away from outer edges of mask (width of bounding box)
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
%make all zero values in the original volume nonzero
minval = max(vol(:))/1000;
vol(find(vol == 0)) = minval;

%{
5 = nodes
4 = fat
3 = [unused]
2 = lumen/healthy rectum
1 = cancer
%}
[row, col, slice] = ind2sub(size(mask),find(mask ~= 0)); 

%find dimensions of bounding box
bottom = min(max(row)+n,size(vol,2));
top = max(min(row)-n,1);
left = max(min(col)-n,1);
right = min(max(col)+n,size(vol,2));
front = max(min(slice)-1,1);
back = min(max(slice)+1,size(vol,3));

%apply bounding box
fprintf('Cropping volume and mask with padding of %i pixels\n',n);
croppedV = vol(top:bottom,left:right,front:back);
croppedM = mask(top:bottom,left:right,front:back);

varargout{1} = croppedV;
if nargout == 2
    varargout{2} = croppedM;
end

end


