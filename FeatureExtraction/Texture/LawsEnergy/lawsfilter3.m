function [lawsresponses, K3names] = lawsfilter3(I,ws,shortflag)
% LAWSFILTER3 Apply 3D Law's filters to an image.
%
%   LAWSRESPONSES = LAWSFILTER3(I,OPTIONS) returns in LAWSRESPONSES a M-by-M-by-K 
%       array containing filter responses from volume I due convolution 
%       with the the laws kernels.
%       
%       OPTIONS:
%           WS: window size of kernels. Either 3 or 5 (default:5)
%           SHORTFLAG: option to reduce number of performed convolutions to
%                      10.Either 0 or 1 (default: 0)
%
%    [LAWSRESPONSES K3NAMES] = LAWSFILTER3(I,OPTIONS) additionally returns
%       the  names of kernel sets from multiplying permutations of 1D kernels

%   See also: lawskerns2, convn, lawskerns3

[nrows, ncols, nplanes, highdims]=size(I);
if highdims>1, error('3D grayscale volumes only.'); end

if nargin>=2
    [K3, K3names]=lawskerns3(ws);
else
    [K3, K3names]=lawskerns3;
end

% [ks1 ks2 ks3 nkerns]=size(K3);
if nargin==3 && shortflag, 
    nkerns=10;
else
    nkerns=size(K3,4); % nkerns 3D kernels
end

% Calculate filter responses
lawsresponses=zeros([nrows ncols nplanes nkerns]);
for i=1:nkerns,
    lawsresponses(:,:,:,i)=convn(I,K3(:,:,:,i),'same');
    if mod(i,4)==0,fprintf('.');end
end
fprintf('\n');