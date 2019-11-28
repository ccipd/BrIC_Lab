function [gb_c, gb_s]=gabor3dkerns(xy_orient,xz_orient,wavelength,bandwidth)
% GABOR	3DKERNS Generate 3D Gabor kernels.
%   [GB_C,GB_S]=GABOR3DKERNS(XY_ORIENT,XZ_ORIENT,LAMBDA,B) calculates the cosine
%   and sine (pi/2 phase difference) Gabor kernels for the specified
%   parameter values.
%	Combine these RESPONSES to get the magnitude (ABS) and PHASE (ATAN(sin/cos)) features for a particular filter
%
%	Vary XY_ORIENT and XZ_ORIENT to acheive proper rotation in the
%corresponding planes
%   LAMBDA = Wavelength
%   B = Bandwidth (set B=1 to reduce params, meaning SIGMA=0.56*LAMBDA)
%	To reduce complication, this functions assumes you only want to construct an isotropic filter (sigma is same in X,Y,Z)
%
%Satish Viswanath, Jul 2008
%
% see also: gabor2dkerns.m

theta = xy_orient;
phi = xz_orient;
lambda = wavelength;
b = bandwidth;

%Sigma based off lambda from
%http://matlabserver.cs.rug.nl/edgedetectionweb/web/edgedetection_params.html
%If "bandwidth" (b) is 1, sigma = 0.56*lambda;
%This is an isotropic filter, so sigma is the same in X,Y,Z 
sigma = (1/pi)*sqrt(log(2)/2)*((2^b+1)/(2^b-1))*lambda;

%Define half window size of filter based off sigma
hws = 3*sigma;
if mod(hws,2)==0, hws=hws+1;end
[x y z]=meshgrid(-hws:hws,-hws:hws,-hws:hws);

% Rotation first by theta, then by phi
x_theta=x*cos(theta)-y*sin(theta);
y_theta=x*sin(theta)+y*cos(theta);
 
x_theta_phi=x_theta*cos(phi)-z*sin(phi);
z_phi=x_theta*sin(phi)+z*cos(phi);

%Only "real" part of filter calculated, imaginary results in weird filter
%2D
%gb=exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);

%3D
gb_c=exp(-.5*(x_theta_phi.^2/sigma^2+y_theta.^2/sigma^2+z_phi.^2/sigma^2)).*cos(2*pi/lambda*x_theta_phi);
gb_s=exp(-.5*(x_theta_phi.^2/sigma^2+y_theta.^2/sigma^2+z_phi.^2/sigma^2)).*sin(2*pi/lambda*x_theta_phi);

%gb = sqrt(gb_r.^2+gb_i.^2); % sample combined response

% To view the filter
% g = abs(fftshift(fftn(gb,[200 200 200])));
% g1 = uint8(255*rescale(g));
% for i=1:200
% image(g1(:,:,i));colormap gray; axis off;pause;
% end

end
