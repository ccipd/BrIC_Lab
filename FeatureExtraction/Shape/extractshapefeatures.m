%This code has been released by CCIPD, Case Western Reserve University
%(Cleveland OH)
function [Ftrs]=extractshapefeatures(Mask,PixelDim)
%filehdr1=' 1-rows, 2-cols, 3-slices, 4-Area_mu, 5-Area_std';
%filehdr2= '6-Perim_mu, 7-Perim_std, 8-EquivDiam_mu, 9-EquivDiam_std, 10-Eccen_mu, 11-Eccen_std';
%filehdr3= '12-Ext_mu, 13-Ext_std, 14-Compactness_mu, 15-Compactness_std, 16-RdDistStd_mu, 17-RdDistStd_std';
%filehdr4= '18-Roughness_mu, 19-Roughness_std, 20-Elongation_mu, 21-Elongation_std';
%filehdr5= '22-Convexity_mu, 23-Convexity_std, 24-sphe ';

%-- Pixel Dimensions
pdimx=PixelDim(1);
pdimy=PixelDim(2);
pdimz=PixelDim(3);
nFeature=24;    % Number of Features

Ftrs=zeros(1,nFeature);
x=double(Mask);

boundrycount=0;
areacount=0;
CC3d = bwconncomp(x);
S3d = regionprops(CC3d,'Centroid','BoundingBox');
for numobj=1:1 %CC3d.NumObjects
    if S3d(numobj).BoundingBox(6)>1
        feat=zeros(S3d(numobj).BoundingBox(6),10);
        for ly=1:S3d(numobj).BoundingBox(6)
            x1=x(:,:,round(S3d(numobj).BoundingBox(3))+ly-1);
            CC = bwconncomp(x1);
            %        objoldflag=0;
            
            if(CC.NumObjects==1)
                %            objflag=1;
                S = regionprops(CC,'Eccentricity','Extent','ConvexArea','EquivDiameter','Area','MajorAxisLength','MinorAxisLength','Perimeter','ConvexImage','ConvexHull','BoundingBox','Centroid');
                blobMeasurements =S.Centroid;
                
                CVh = bwconncomp(S.ConvexImage);
                Cvx = regionprops(CVh,'Area','MajorAxisLength','MinorAxisLength','Perimeter','BoundingBox','Centroid');
                
                % Get the centroid.
                centroidX = blobMeasurements(1,1);
                centroidY = blobMeasurements(1,2);
                % Get the boundary.
                boundaries = bwboundaries(x1);
                thisBoundary = boundaries{1};
                ConvexImageboundry = bwboundaries(S.ConvexImage);
                % Get the distances of the boundary pixels from the centroid.
                distances = sqrt((thisBoundary(:,1) - centroidY).^2 + (thisBoundary(:,2) - centroidX).^2);
                
                RdialDistanceStd=std(distances);
                Roughness=S.Perimeter/Cvx.Perimeter;
                Elongation=S.MajorAxisLength/S.MinorAxisLength;
                Compactness=(S.Perimeter^2/(4*pi*S.Area)) ;
                Convexity=(S.Area/Cvx.Area);
                % feature vector of each slice
                feat(ly,:)=[S.Area*pdimx*pdimy, S.Perimeter*pdimx,S.EquivDiameter *pdimx, S.Eccentricity*pdimx*pdimy, S.Extent*pdimx*pdimy,Compactness*pdimx*pdimy,RdialDistanceStd*pdimx,Roughness*pdimx*pdimy,Elongation*pdimx*pdimy,Convexity*pdimx*pdimy];
                
                areacount =areacount+ sum(sum(x1(:,:)) );
                boundrycount=boundrycount+size(thisBoundary,1);
            end
        end
        sphe= ((pi^(1/3))* (6 * areacount)^(2/3) /boundrycount)*pdimx*pdimy*pdimz;
        %final feature vector of each nodule
        Ftrs=[round(S3d(numobj).BoundingBox(4))*pdimy,round(S3d(numobj).BoundingBox(5))*pdimx,round(S3d(numobj).BoundingBox(6))*pdimz,mean(feat(:,1)),std(feat(:,1)),mean(feat(:,2)),...
            std(feat(:,2)),mean(feat(:,3)),std(feat(:,3)),mean(feat(:,4)),std(feat(:,4)),...
            mean(feat(:,5)),std(feat(:,5)),mean(feat(:,6)),std(feat(:,6)),mean(feat(:,7)),...
            std(feat(:,7)),mean(feat(:,8)),std(feat(:,8)),mean(feat(:,9)),std(feat(:,9)),...
            mean(feat(:,10)),std(feat(:,10)),sphe];
        clear feat;clear x1;clear CC;clear S;clear CVh;clear Cvx;
    end
end




