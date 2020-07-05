function [featints,featnames,featstats,statnames] = extract2DFeatureInfo(img,mask,class_options,ws_options)
    % Feature extraction for a series of 2D images
    
    %INPUTS:
    %img = 3D matrix of set of 2D slices OR a single 2D image (double format)
    %mask = corresponding set of 2D label(s)same size as img (ROI = 1)
    %class_options = (OPTIONAL) cell array of strings corresponding to desired feature classes:
    %                 DEFAULT: class_options = {'raw','gray','gradient','haralick','gabor','laws','collage'};
    %ws_options = (OPTIONAL) array of integers corresponding to desired window levels:
    %                 DEFAULT: ws_options = [3, 5, 7, 9, 11]
   
    %OUTPUTS:
    %featints = 1D cell array of feature intensity values from within ROI of feature slice representations of img
    %featnames = 1D cell array of feature names
    %featstats = 2D matrix of statistics for each column of featints
    %statnames = 2D matrix of statistics names

    %see subfunctions:
    % function varargout = extract2DFeatIntensities(featVol, mask, vals)
    % function [featstats,statnames]=   compute2DFeatStats(featints,featnames)
    
    % @ Jacob Antunes, 2016
    % Last Updated: 06-18-2018
    
%% Setup overarching pathway needed

funcname = 'extract2DFeatureInfo.m';
funcpath = which(funcname);
codepath = funcpath(1:end-length(funcname));

addpath([codepath '../images']);
addpath([codepath '../mha']);
addpath([codepath '../pre-processing']);

%% Check inputs

if isempty(img) || isempty(mask)
    error('Check inputs. Either IMG or MASK is an empty matrix.');
end

if ~isequal(size(img),size(mask))
    error('Check inputs. Size of MASK does not equal size of IMG');
end

if isempty(find(mask>0,1))
    error('Mask is empty. Nothing to extract.');
end

%% Initialization
rawints = []; rawstats = [];
grayints = []; graystats = [];
gradints = []; gradstats = [];
harints = []; harstats = [];
gabints = []; gabstats = [];
lawints = []; lawstats = [];
collageints = []; collagestats = [];

if exist('class_options','var')~=1 || isempty(class_options)
    class_options = {'raw','gray','gradient','haralick','gabor','laws','collage'};
end

if ~isa(class_options,'cell')
    class_options = mat2cell(class_options,1,length(class_options));
end

if exist('ws_options','var')~=1 || isempty(ws_options)
    ws_options = 3:2:11;
end

ws1 = ws_options(1);
ws2 = ws_options(end);

%% RECOMMENDED: CROP IMAGE AND MASK (saves time and memory)!
[img,mask] = boundingbox2(img,mask,max(ws_options),'allz','off');

%% 2D Feature Intensity Extraction


for i=1:size(img,3)
    
    if size(img,3) > 1
        fprintf('slice #%i of %i:\n',i,size(img,3));pause(0.5);
    end
    
    if isempty(find(mask(:,:,i)==1,1))
        fprintf('No mask on this slice. Continuing.\n');
        continue;
    end
            
%--------------Raw Intensity -----------------%
if ~all(cellfun('isempty',strfind(class_options,'raw')))
    fprintf('\nExtracting raw intensity information:\n');
    
    rawints = [rawints; cell2mat(extract2DFeatIntensities(img(:,:,i), mask(:,:,i), 1))];
    rawnames = {'Raw intensity'};
end

%--------------Gray Level Statistics----------------%
if ~all(cellfun('isempty',strfind(class_options,'gray')))
    fprintf('\nExtracting Gray Level Statistics:\n');

    addpath([codepath,'grey']);   
    grayfeats = [];
    graynames = {};
    for ws = ws_options
        fprintf('\tUsing a window size of %i...\n',ws)
        [gf,gfn] = grayfilts2(img(:,:,i),ws);
        grayfeats=cat(3,grayfeats,gf);
        for s = 1:length(gfn)
            gfn{s} = strcat(gfn{s},' ws=',num2str(ws));
        end
        graynames(end+1:end+length(gfn))=gfn;
    end
    grayints = [grayints; cell2mat(extract2DFeatIntensities(grayfeats, mask(:,:,i), 1))];
    clear grayfeats gfn ws gf s;
    rmpath([codepath,'grey']);   
end

%--------------Gradient--------------%
if ~all(cellfun('isempty',strfind(class_options,'gradient')))
    fprintf('\nExtracting Gradient Features:\n');
    
    if any(ws_options>3)
        fprintf('\tPlease note that Gradient features cannot use a ws > 3 currently.\n');
    end
    
    addpath([codepath,'grey']); 
    
    if ws1 == 3
        [gradfeats, gradnames] = gradfilts2(img(:,:,i));
        gradints = [gradints; cell2mat(extract2DFeatIntensities(gradfeats, mask(:,:,i), 1))];
    end
    
    clear gradfeats;
    rmpath([codepath,'grey']);   
end

%--------------Haralick--------------%
if ~all(cellfun('isempty',strfind(class_options,'haralick')))
    fprintf('\nExtracting Haralick Features:\n');

    addpath([codepath,'haralick']);
    haralickfun = @haralick2mex; % haralick2mexmt or, haralick2mex

    hfeatnames3 = {'Haralick entropy ws=3' 'Haralick energy ws=3' 'Haralick inertia ws=3' 'Haralick idm ws=3' 'Haralick correlation ws=3' 'Haralick info1 ws=3' ...
    'Haralick info2 ws=3' 'Haralick sum-av ws=3' 'Haralick sum-var ws=3' 'Haralick sum-ent ws=3' 'Haralick diff-av ws=3' 'Haralick diff-var ws=3' 'Haralick diff-ent ws=3'};
    hfeatnames5 = {'Haralick entropy ws=5' 'Haralick energy ws=5' 'Haralick inertia ws=5' 'Haralick idm ws=5' 'Haralick correlation ws=5' 'Haralick info1 ws=5' ...
    'Haralick info2 ws=5' 'Haralick sum-av ws=5' 'Haralick sum-var ws=5' 'Haralick sum-ent ws=5' 'Haralick diff-av ws=5' 'Haralick diff-var ws=5' 'Haralick diff-ent ws=5'};
    hfeatnames7 = {'Haralick entropy ws=7' 'Haralick energy ws=7' 'Haralick inertia ws=7' 'Haralick idm ws=7' 'Haralick correlation ws=7' 'Haralick info1 ws=7' ...
    'Haralick info2 ws=7' 'Haralick sum-av ws=7' 'Haralick sum-var ws=7' 'Haralick sum-ent ws=7' 'Haralick diff-av ws=7' 'Haralick diff-var ws=7' 'Haralick diff-ent ws=7'};
    hfeatnames9 = {'Haralick entropy ws=9' 'Haralick energy ws=9' 'Haralick inertia ws=9' 'Haralick idm ws=9' 'Haralick correlation ws=9' 'Haralick info1 ws=9' ...
    'Haralick info2 ws=9' 'Haralick sum-av ws=9' 'Haralick sum-var ws=9' 'Haralick sum-ent ws=9' 'Haralick diff-av ws=9' 'Haralick diff-var ws=9' 'Haralick diff-ent ws=9'};
    hfeatnames11 = {'Haralick entropy ws=11' 'Haralick energy ws=11' 'Haralick inertia ws=11' 'Haralick idm ws=11' 'Haralick correlation ws=11' 'Haralick info1 ws=11' ...
    'Haralick info2 ws=11' 'Haralick sum-av ws=11' 'Haralick sum-var ws=11' 'Haralick sum-ent ws=11' 'Haralick diff-av ws=11' 'Haralick diff-var ws=11' 'Haralick diff-ent ws=11'};
    
    harnames = {};
    for ws = ws_options
        harnames=cat(2,harnames,eval(['hfeatnames' int2str(ws)]));    
    end
    
    nbins = 128; % a power of 2
    bkgrd = -1; %b/c rescaling from 0 to 1
    dist = 1; %d = 0 (co-occurence with just self, d = 1 (co-occurence with just immediate neighbors), d = 2 (co-occurence with more distant neighbors)
    
    norm_img = round(rescale_range(img(:,:,i),0,nbins-1));
%     if length(unique(norm_img))<nbins
%         error('Not enough gray levels present in IMG. Try resampling or not cropping.');
%     end
    hfeats = [];
    for ws = ws_options
        fprintf('\tUsing a window size of %i...\n',ws)
        hfeats = cat(3,hfeats,haralickfun(norm_img,nbins,ws,dist,bkgrd));
    end

    harints = [harints; cell2mat(extract2DFeatIntensities(hfeats, mask(:,:,i), 1))];
    clear hfeats hfeatnames3 hfeatnames5 hfeatnames7 hfeatnames9 hfeatnames11 norm_img ws;
    rmpath([codepath,'haralick']);   
end

%--------------Gabor--------------%
if ~all(cellfun('isempty',strfind(class_options,'gabor')))
    fprintf('\nExtracting Gabor Wavelets: ');
    
    addpath([codepath,'gabor']);
    load gaborbank2D; %gabor_data
    
    fcount = 0;
    for f=1:length(gabor_data) %to save memory, dont run all of Gabor!)
        if any(ismember(ws_options, gabor_data(f).ws))
            fcount = fcount+1;
            g_c = jconvn(img(:,:,i),gabor_data(f).cos);
            g_s = jconvn(img(:,:,i),gabor_data(f).sin);
            gabfeats(:,:,fcount) = sqrt(g_c.*g_c + g_s.*g_s);
            gabnames(fcount) = {gabor_data(f).title};
        end
        if mod(f,4) == 0, fprintf('.'); end
    end
    fprintf('\n');
    gabints = [gabints; cell2mat(extract2DFeatIntensities(gabfeats, mask(:,:,i), 1))];
    clear gabfeats g_c g_s gabor_data f;
    rmpath([codepath,'gabor']);
end

%--------------Laws--------------%
if ~all(cellfun('isempty',strfind(class_options,'laws')))

    fprintf('\nExtracting Laws Filters:\n');

    addpath([codepath,'laws']);
    addpath([codepath,fullfile('laws','mtimesx')]);
    
    lawfeats = [];
    lawnames = {};
    
    if any(ws_options>5)
        fprintf('\tPlease note that Laws features cannot use a ws > 5 currently.\n');     
    end
    
    for ws = ws1:2:ws2
        if ws <= 5
            fprintf('\tUsing a window size of %i:\n',ws)
            [lf, lfn] = lawsfilter(img(:,:,i),ws);
            lawfeats = cat(3,lawfeats,lf);
            lawnames=cat(2,lawnames,lfn);
        end
    end
    
    lawints = [lawints; cell2mat(extract2DFeatIntensities(lawfeats, mask(:,:,i), 1))];
    clear lf lfn lawfeats lawfeatnames ws;
    rmpath([codepath,'laws']); 
    rmpath([codepath,fullfile('laws','mtimesx')]);
end

%--------------Collage--------------%
if ~all(cellfun('isempty',strfind(class_options,'collage')))

  fprintf('\nExtracting COLLAGE features:\n');
  
  addpath([codepath,'collage']);
  addpath([codepath,'haralick']);  % we need the haralick code to calculate collage features  

  if any(ws_options>5)
        fprintf('\tPlease note that COLLAGE features cannot use a ws > 5 currently.\n');     
  end

collagefeats = [];
collagenames = {};

  for ws=ws1:2:ws2
      if ws <= 5
          fprintf('\tUsing a window size of %i:\n',ws)
          switch ws
              case 3
                  k = 1; 
                  C_hfeatnames = {'CoLlAGe entropy ws=3' 'CoLlAGe energy ws=3' 'CoLlAGe inertia ws=3' 'CoLlAGe idm ws=3' 'CoLlAGe correlation ws=3' 'CoLlAGe info1 ws=3' ...
                                  'CoLlAGe info2 ws=3' 'CoLlAGe sum-av ws=3' 'CoLlAGe sum-var ws=3' 'CoLlAGe sum-ent ws=3' 'CoLlAGe diff-av ws=3' 'CoLlAGe diff-var ws=3' 'CoLlAGe diff-ent ws=3'};
              case 5
                  k = 2;
                  C_hfeatnames = {'CoLlAGe entropy ws=5' 'CoLlAGe energy ws=5' 'CoLlAGe inertia ws=5' 'CoLlAGe idm ws=5' 'CoLlAGe correlation ws=5' 'CoLlAGe info1 ws=5' ...
                                  'CoLlAGe info2 ws=5' 'CoLlAGe sum-av ws=5' 'CoLlAGe sum-var ws=5' 'CoLlAGe sum-ent ws=5' 'CoLlAGe diff-av ws=5' 'CoLlAGe diff-var ws=5' 'CoLlAGe diff-ent ws=5'};
           end
          
          [o1,~]=compute_CoLlAGe2D(img(:,:,i),mask(:,:,i),k,1:13);
          
          collagefeats = cat(3,collagefeats,o1);
          collagenames = cat(2,collagenames,C_hfeatnames);
          end
    end
         collageints = [collageints; cell2mat(extract2DFeatIntensities(collagefeats,mask(:,:,i),1))]; 
         rmpath([codepath,'haralick']); 
         rmpath([codepath,'collage']);
         clear collagefeats C_hfeatnames o1 k
end

end %move on to next slice

fprintf('\n');

%% append featints and featnames
startz = 0;
endz   = 0;

if ~all(cellfun('isempty',strfind(class_options,'raw'))) && ~isempty(rawints)          
    startz = endz+1;
    endz   = endz+size(rawints,2);
    featints(startz:endz) = num2cell(rawints,1);
    featnames(startz:endz) = rawnames;
    disp(['Added intensity for original image!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gray'))) && ~isempty(grayints)
    startz = endz+1;
    endz   = endz+size(grayints,2);
    featints(startz:endz) = num2cell(grayints,1); 
    featnames(startz:endz) = graynames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gray features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gradient'))) && ~isempty(gradints)
    startz = endz+1;
    endz   = endz+size(gradints,2);
    featints(startz:endz) = num2cell(gradints,1); 
    featnames(startz:endz) = gradnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gradient features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'haralick'))) && ~isempty(harints)
    startz = endz+1;
    endz   = endz+size(harints,2);
    featints(startz:endz) = num2cell(harints,1); 
    featnames(startz:endz) = harnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' haralick features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gabor'))) && ~isempty(gabints)
    startz = endz+1;
    endz   = endz+size(gabints,2);
    featints(startz:endz) = num2cell(gabints,1); 
    featnames(startz:endz) = gabnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gabor features!']) 
end

if ~all(cellfun('isempty',strfind(class_options,'laws'))) && ~isempty(lawints)
    startz = endz+1;
    endz   = endz+size(lawints,2);
    featints(startz:endz) = num2cell(lawints,1); 
    featnames(startz:endz) = lawnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' laws features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'collage')))
    startz = endz+1;
    endz   = endz+size(collageints,2);
    featints(startz:endz) = num2cell(collageints,1); 
    featnames(startz:endz) = collagenames;
    disp(['Added intensities for ', int2str(endz-startz+1),' collage features!'])
end

disp(['In total: extracted intensities for ', int2str(length(featints)),' texture features!']); 
    
%% Compute, Append 2D Feature statistics and statnames

if nargout <3 %user does not want to compute statistics of feature distributions
    return;
end

featstats = zeros(size(rawstats,1),size(rawstats,2)+size(graystats,2)+size(gradstats,2)+size(harstats,2)+size(gabstats,2)+size(lawstats,2));

startz = 0;
endz = 0;
    
if ~all(cellfun('isempty',strfind(class_options,'raw'))) && ~isempty(rawints)
    [rawstats,rawstatnames] = compute2DFeatStats(num2cell(rawints,1),rawnames);
    startz = endz + 1;
    endz   = size(rawstats,2);
    featstats(:,startz:endz) = rawstats; 
    statnames(:,startz:endz) = rawstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' raw intensity!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gray'))) && ~isempty(grayints)
    [graystats,graystatnames] = compute2DFeatStats(num2cell(grayints,1),graynames);
    startz = endz+1;
    endz   = endz+size(graystats,2);
    featstats(:,startz:endz) = graystats; 
    statnames(:,startz:endz) = graystatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gray features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gradient'))) && ~isempty(gradints)
    [gradstats,gradstatnames] = compute2DFeatStats(num2cell(gradints,1),gradnames);
    startz = endz+1;
    endz   = endz+size(gradstats,2);
    featstats(:,startz:endz) = gradstats; 
    statnames(:,startz:endz) = gradstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' grad features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'haralick'))) && ~isempty(harints)
   
    [harstats,harstatnames] = compute2DFeatStats(num2cell(harints,1),harnames);
    startz = endz+1;
    endz   = endz+size(harstats,2);
    featstats(:,startz:endz) = harstats; 
    statnames(:,startz:endz) = harstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' haralick features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gabor'))) && ~isempty(gabints)
    [gabstats,gabstatnames] = compute2DFeatStats(num2cell(gabints,1),gabnames);
    startz = endz+1;
    endz   = endz+size(gabstats,2);
    featstats(:,startz:endz) = gabstats; 
    statnames(:,startz:endz) = gabstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gabor features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'laws'))) && ~isempty(lawints)
    [lawstats,lawstatnames] = compute2DFeatStats(num2cell(lawints,1),lawnames);
    startz = endz+1;
    endz = endz+size(lawstats,2);
    featstats(:,startz:endz) = lawstats; 
    statnames(:,startz:endz) = lawstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' laws features!'])

end

if ~all(cellfun('isempty',strfind(class_options,'collage'))) 
    [collagestats,collagestatnames] = compute2DFeatStats(num2cell(collageints,1),collagenames);
    startz = endz+1;
    endz = endz+size(collagestats,2);
    featstats(:,startz:endz) = collagestats; 
    statnames(:,startz:endz) = collagestatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' collage features!'])
end

disp(['In total: computed ' int2str(size(featstats,1)) ' statistics for ', int2str(size(featstats,2)),' texture features!']); 

clear rawints grayints gradints harints gabints lawints;
clear rawnames graynames gradnames hnames gabnames lawnames;
clear rawstats graystats gradstats harstats gabstats lawstats;
clear rawstatnames graystatnames gradstatnames harstatnames gabstatnames lawstatnames;

end
 
%% extract2DFeatIntensities subfunction
% Jacob Antunes
% April 11, 2017

function varargout = extract2DFeatIntensities(featVol, mask, vals)
%INPUTS
% featVol = 3-D volume of texture features
% mask = 2-D annotation image
% vals = vector of mask values to find intensities within

% OUTPUT
% varargout = cell array(s) of intensities. Number based on length of vals


    fprintf('\n\tExtracting intensities\n');
    
    if isempty(featVol) || isempty(mask) || isempty(vals)
        varargout = {};
        return;
    end
    
    ints = cell(size(featVol,3),length(vals)); %cell array of intensities
    
    for i = 1:size(featVol,3) % loop through all features
        f = featVol(:,:,i); %single feature volume
        for j = 1:length(vals) % loop through all mask values
            ints{i,j} = f(mask == vals(j));
        end
    end
    
    for k = 1:length(vals)
        varargout{k} = ints(:,k)';
    end
    
end

%% compute2DFeatStats subfunction
% Jacob Antunes
% April 11, 2017

function [featstats,statnames]=compute2DFeatStats(featints, featnames)
%INPUTS
% featints: a 1xM cell array, where each column cell is a vector containing feature values
% featnames: (optional) a 1xM cell array of the names corresponding to each element of featints; use [] if you dont care about knowing statnames

%OUTPUT
% featstats: a 4xM double matrix, where each element in a column a single of the feature intensity vector in the corresponding column
% statnames: a 4xM cell matrix, where each element is the statistic and feature name corresponding to featstats

% fprintf('Extracting 4 feature statistics\n');
featstats=zeros(4,length(featints)); %mxn matrix (m = # statistics, n = # texture features)
statnames=cell(4,length(featints));
for i=1:length(featints)
    temp1=featints{i};
    featstats(:,i)= [median(temp1); ...
                     std(temp1)^2; ...
                     kurtosis(temp1); ...
                     skewness(temp1)];
                     %previously used: mean(temp1);range(temp1); std(temp1); harmmean(temp1);mad(temp1); iqr(temp1)
    clear temp1
end
if ~isempty(featnames)
   for i=1:length(featnames)
        statnames{1,i} = strcat('median-',featnames{i});
        statnames{2,i} = strcat('var-',featnames{i});
        statnames{3,i} = strcat('kurtosis-',featnames{i});
        statnames{4,i} = strcat('skewness-',featnames{i});
%          previously used:
%          statnames{3,i} = strcat('range-',featnames{i});
%         statnames{7,i} = strcat('harmean-',featnames{i});
%         statnames{4,i} = strcat('std-',featnames{i});
%          statnames{5,i} = strcat('mean-',featnames{i});
%         statnames{9,i} = strcat('mad-',featnames{i});
%         statnames{10,i} = strcat('iqr-',featnames{i});
   end
end 

end

 
