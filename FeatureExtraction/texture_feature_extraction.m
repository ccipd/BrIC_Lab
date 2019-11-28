 %Feature Extraction
%This code has been released by CCIPD, Case Western Reserve University
%(Cleveland OH)

%This script is used for feature extraction:
% haralick features with window size of 5
% gabor filter response features
% laws filter response features (5x5x5 kernel)

%Written by Jacob Antunes
%Edited by Niha Beig

clc;
clear;
close all;
%add necessary files needed
addpath(genpath('~\ccipd')); %add path to github repository with dependencies
%haralickfun=@haralick2mexmt; %needed to run haralick
addpath(genpath('~\data'));
cd('~\data');
atlas_mask=load_untouch_nii('~\sri24_spm8\templates\mask-label.nii.gz');
mask=logical(atlas_mask.img);


%% input studies
studies={}; %load patient folder names here
%%
i=size(studies,1);
tic

for j=1:i
    T1={}; %change for Gd-T1w, T2w and T2w-FLAIR images
    s=studies(j,1);
    new_folder = char(s);
    T1_final{j,1}=char(s);
    y=sprintf('~\\%s',new_folder);
    cd(y);
    fprintf('##########################################################\n');
    disp(new_folder)
    b=struct2cell(dir(y));
    %This script runs for a specific MR sequence (i.e., Gd-T1w, T2w or
    %FLAIR)
  if (any(ismember(b(1,:),'t1Gd.nii.gz'))==1)
        main=load_untouch_nii('t1Gd.nii.gz');
        if (any(ismember(b(1,:),'nec-label.nii.gz'))==1)
            label1=load_untouch_nii('nec-label.nii.gz');
            tmp1=logical(label1.img);
            [X1,X2] = bounding_Box(main.img,tmp1,5);
            [featstats,namestats] = extract3DFeatures(double(X1),X2);
            T1_final{j,2}=featstats;
       else
            fprintf('\n%%%%%%%%%%%%\nNo Necrosis label found');
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (any(ismember(b(1,:),'edema-label.nii.gz'))==1)
            label2=load_untouch_nii('edema-label.nii.gz');
            tmp2=logical(label2.img);
            [X1,X2] = bounding_Box(main.img,tmp2,5);
            [featstats] = extract3DFeatures(double(X1),X2);
            T1_final{j,3}=featstats;
       else
            fprintf('\n%%%%%%%%%%%%\nNo Edema label found');
       end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (any(ismember(b(1,:),'et-label.nii.gz'))==1)
            label3=load_untouch_nii('et-label.nii.gz');
            tmp3=logical(label3.img);
            [X1,X2] = bounding_Box(main.img,tmp3,5);
            [featstats] = extract3DFeatures(double(X1),X2);
            T1_final{j,4}=featstats;
       else
            fprintf('\n%%%%%%%%%%%%\nNo enhancing tumor label found');
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new = strcat(new_folder,'_T1');
    T1{1,1}=T1_final{j,1};
    T1{1,2}=T1_final{j,2};
    T1{1,3}=T1_final{j,3};
    T1{1,4}=T1_final{j,4};
    save(new,'T1')
        cd ..
        cd ..
    else
        fprintf('T1 absent');
  end
end

%saving the results
cd('~\results');
        save('T1_final','T1_final');
        
toc
