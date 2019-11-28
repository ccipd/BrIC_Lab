% Shape Feature Extraction
%This code has been released by CCIPD, Case Western Reserve University
%(Cleveland OH)
%This script is used for shape feature extraction
%Written by Mahdi Orooji
% Orooji et al., "Combination of computer extracted shape and texture features enables discrimination of granulomas from adenocarcinoma on chest computed tomography," J. Med. Imag. 5(2) 024501 (18 April 2018)
%Edited by Niha Beig

clc;
clear;
close all;
%add necessary files needed
addpath(genpath('~\ccipd')); %add path to github repository with dependencies
addpath(genpath('~\data')); %add path to data
cd('~\data');
atlas_mask=load_untouch_nii('~\sri24_spm8\templates\mask-label.nii.gz'); %add path to the atlas mask
mask=logical(atlas_mask.img);

%% input studies
shape_studies={}; %load patient folder names here
%%
i=size(shape_studies,1);
tic

for j=1:i
    s=shape_studies(j,1);
    new_folder = char(s);
    shape_final{j,1}=char(s);
    y=sprintf('~\\%s',new_folder);
    cd(y);
    fprintf('##########################################################\n');
    disp(new_folder)
    b=struct2cell(dir(y));
    %This script runs for specified MR sequence. 
    %For edema and enhancing tumor region
       if (any(ismember(b(1,:),'edema-label.nii.gz'))==1)
            fprintf('\n%%%%%%%%%%%%\nEdema label found');
            label2=load_untouch_nii('edema-label.nii.gz'); %download this Nifti reader MATLAB library
            tmp2=logical(label2.img);
            featstats2=extractshapefeatures(tmp2,[1 1 1]);
            shape_final{j,2}=featstats2;
       else
            fprintf('\n%%%%%%%%%%%%\nNo Edema label found');
       end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (any(ismember(b(1,:),'et-label.nii.gz'))==1)
            label3=load_untouch_nii('et-label.nii.gz');
            fprintf('\n%%%%%%%%%%%%\nEnhancing Tumor label found');
            tmp3=logical(label3.img);

            featstats=extractshapefeatures(tmp3,[1 1 1]);
            shape_final{j,3}=featstats;
       else
            fprintf('\n%%%%%%%%%%%%\nNo enhancing tumor label found');
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new = strcat(new_folder,'_shape');
    shape_train{1,1}=shape_final{j,1};
    shape_train{1,2}=shape_final{j,2};
    shape_train{1,3}=shape_final{j,3};
    save(new,'shape_train')
        cd ..
        cd ..
   % else
       % fprintf('FLAIR absent');
end
%end

%saving the results
cd('~\results');
        save('shape_train','shape_final');
        
toc
