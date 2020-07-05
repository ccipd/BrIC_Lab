%% Test code for feature extraction functions

close all;clear;clc;

% parameters to set
class_options = {'collage'};
ws_options = [3,5];%:2;11; %(optional) window sizes from which to extract


%% Add path of supporting codes
codepath = pwd
addpath(genpath([codepath '../Collage_supporting_codes']));

%% Read data and extract Collage features
dirname = 'DATA PATH';
listing = dir(dirname);
patient_names = listing(3:end);

for i =1:length(patient_names)
 new_path = strcat(dirname,'\', patient_names(i).name);
    image = load_nii(strcat(new_path,'\T1.nii'));%T1
    vol = double(abs(image.img));
    label = load_nii(strcat(new_path,'\annotated_file.nii'));
    label = label.img;
    mask = (label>0);
    
[featints, featnames, featstats, statnames] = extract3DFeatureInfo(vol,mask,class_options, ws_options);
 patient_stats= reshape(featstats,[1,numel(featstats)]); 
 Texture_features(:,i) = patient_stats';
end
Texture_feature_names = statnames(:);
