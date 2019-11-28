function [featstats,statnames]=compute_FeatStats(feats,roi,featnames)
%feats must be either a 4D matrix (multiple 3D volumes!) OR 1 3D volume
%ROI: should be 3D label mask
%featnames: optional; use [] if you dont care about knowing statnames

fprintf('Extracting 6 feature statistics\n');
featstats=[]; %mxn matrix (m = # statistics, n = # texture features)
    
if ndims(feats)==4
    for i=1:size(feats,4)
        temp1=feats(:,:,:,i);
        temp1=temp1(roi>0); %currently only works for 1 label          
                     featstats(:,i)= [median(temp1);std(temp1);...
                        kurtosis(temp1); skewness(temp1)];
        clear temp1
    end
    if ~isempty(featnames)
       for i=1:size(feats,4)
            statnames{1,i} = strcat('median-',featnames{i});
            statnames{2,i} = strcat('std-',featnames{i});
            statnames{3,i} = strcat('kurtosis-',featnames{i});
            statnames{4,i} = strcat('skewness-',featnames{i});
       end
    end 
elseif ndims(feats)==3
        temp1=feats;
        temp1=temp1(roi>0); %currently only works for 1 label
                      featstats = [median(temp1);std(temp1); ...
                         kurtosis(temp1);skewness(temp1)];
    if ~isempty(featnames) %and assuming only 1 featname for 3D input
            statnames{1} = strcat('median-',featnames);
            statnames{2} = strcat('std-',featnames);
            statnames{3} = strcat('kurtosis-',featnames);
            statnames{4} = strcat('skewness-',featnames);
            statnames = statnames';
     end
else
    error('Invalid inputted "feats" variable! Must be 1 3D volume or a 4D matrix of 3D volumes!');
end

end