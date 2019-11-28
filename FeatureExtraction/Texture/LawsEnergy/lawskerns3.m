function [K3, K3names] = lawskerns3(ws)
% LAWSKERNS Generate the Law's 3D kernels from the 1D kernel bases.
%
%   K3 = LAWSKERNS    returns the kernels in a 5-by-5-by-5-by-125 array K3   (default: ws=5).
%
%   K3 = LAWSKERNS(5) returns the kernels in a 5-by-5-by-5-by-125 array K3.
%
%   K3 = LAWSKERNS(3) returns the kernels in a 3-by-3-by-3-by-27 array K3.
%
%   [K3, K3names] = LAWSKERNS(...) returns the names of kernel sets from multiplying permutations of 1D kernels
%
%   Uses mtimesx for 3D (e.g. E3E3E3, E3S3S3, L3S3E3, etc.):
%   L3L3E3=mtimesx(L3',mtimesx(permute(L3,[3 1 2]),E3))  % 3x1x1 * 1x1x3 * 1x3x1
%   OR
%   L3L3E3=mtimesx(permute(mtimesx(L3',L3),[1 3 2]),E3)  % (3x1x1 * 1x3x1)^Transpose(2<->3) * 1x3x1
%
%   ***mtimesx folder found in same folder as lawskerns3.m! (just addpath)
%
%
%   1D Kernels based on Suzuki, M.T., Yaginuma, Y., Kodama, H., 2009. A texture energy measurement technique for 3D volumetric data. In: IEEE International Conference on Systems, Man and Cybernetics, pp. 3779–3785.


if nargin~=1 || isempty(ws)
    ws = 5;
end

K3 = zeros(ws,ws,ws,ws^3); %set of 3D kernels to be created
K3names = cell(1,ws^3); %set of 3D kernels to be created

if ws==3
    %if ws = 3, 1D form of kernels:
    K(1,:)=[-1,0,1]; %E3 (Edge)
    K(2,:)=[1,2,1]; %L3(Level)
    K(3,:)=[-1,2,-1];  %S3 (Spot)
    Knames = {'E3','L3','S3'};  
    
elseif ws==5
    %if ws = 5, 1D form of kernels:
    K(1,:) = [ 1  4  6  4  1]; % L5 (Level) --> L5 = L3*L3
    K(2,:) = [-1 -2  0  2  1]; % E5 (Edge) --> E5 = L3*E3 = E3*L3
    K(3,:) = [-1  0  2  0 -1]; % S5 (Spot) --> S5 = L3*S3 = E3*E3
    K(4,:) = [ 1 -4  6 -4  1]; % R5 (Ripple) --> R5 = S3*L3 = S3*S3
    K(5,:) = [-1  2  0 -2  1]; % W5 (Wave) --> W5 = E3*S3 = S3*E3
    Knames = {'L5','E5','S5','R5','R5','W5'};    

else
    error('Invalid window size! ws must be 3 or 5.');
end


%create 3D kernels
for k=1:ws,
    for i=1:ws,
        for j=1:ws,
            K3(:,:,:,j+(i-1)*ws+(k-1)*ws^2) = mtimesx(K(i,:)',mtimesx(permute(K(j,:),[3 1 2]),K(k,:))); % 2D: K(i,:)'*K(j,:);
            K3names{j+(i-1)*ws+(k-1)*ws^2} = strcat(Knames{i},Knames{j},Knames{k});
        end
    end
end


