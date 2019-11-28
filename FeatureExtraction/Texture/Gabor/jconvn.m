function cI = jconvn(I,K)
% JCONVN N-D convolution in Fourier domain.
%   CI = JCONVN(In,G) convolves In with kernel G and returns CI. The
%   dimensions of CI are the same as I. The array G can be of equal or
%   smaller size than I.
%
%   This function is the same as convn(I,G,'same'), but computed in the
%   Fourier domain. The advantages are speed and robustness to edge
%   effects.
%
%JC

Isize = size(I);
Ksize = size(K);

% nffts=2.^(4:16);
% for i=1:length(Isize),
%     goodnfft(i)=nffts(min(find((nffts-Isize(i)-Ksize(i))>=0)));
% end
% goodnfft(i+1)=Isize(end)+Ksize(end)-1;
% fI = fftn(I,goodnfft);
% fK = fftn(K,goodnfft);

fI = fftn(I,Isize+Ksize-1);
fK = fftn(K,Isize+Ksize-1);
cI = abs(ifftn(fI.*fK));
dims=1:ndims(I);
for i=dims,
    cropinds{i}=1+floor(Ksize(i)/2):Isize(i)+floor(Ksize(i)/2); %#ok<AGROW,NASGU>
end
argstring=sprintf('cropinds{%d}, ',dims);argstring(end-1:end)=[];
eval(['cI = cI(' argstring ');']);
