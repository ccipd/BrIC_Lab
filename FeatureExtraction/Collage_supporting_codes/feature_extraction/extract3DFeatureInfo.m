function [featints,featnames,featstats,statnames] = extract3DFeatureInfo(vol,mask,class_options, ws_options)
    % Feature extraction for a 3D image volume
    
    warning('Not fully tested. Collage may not be correctly concatenating');
    
    %INPUTS:
    %vol = 3D matrix of image volume (double format)
    %mask = corresponding 3D label same size as vol vol(ROI = 1)
    %class_options = (OPTIONAL) cell array of strings corresponding to desired feature classes:
    %                 DEFAULT: class_options = {'raw','gray','gradient','haralick','gabor','laws','collage'};
    %ws_options = (OPTIONAL) array of integers corresponding to desired window levels:
    %                 DEFAULT: ws_options = [3, 5, 7, 9, 11]
   
    %OUTPUTS:
    %featints = 1D cell array of feature intensity values from within ROI of feature volume representations of vol
    %featnames = 1D cell array of feature names
    %featstats = 2D matrix of statistics for each column of featints (10 statistics each)
    %statnames = 2D matrix of statistics names

    
    %see subfunctions:
    % function varargout = extract3DFeatIntensities(featVol, mask, vals)
    % function [featstats,statnames]=compute3DFeatStats(featints,featnames)

    % ******RECOMMENDED: USE ISOTROPIC VOLUMES AND MASKS! (otherwise 3D extraction doesnt make sense) *****

    % @ Jacob Antunes, 2016
    % Last Updated: 08-16-2018

%% Setup overarching pathway needed

funcname = 'extract3DFeatureInfo.m';
funcpath = which(funcname);
codepath = funcpath(1:end-length(funcname));

addpath([codepath '../images']);
addpath([codepath '../subfunctions']);

%% check inputs

if isempty(vol) || isempty(mask)
    error('Check inputs. Either VOL or MASK is an empty matrix.');
end

if numel(size(vol))~=3 || numel(size(mask))~=3
    error('Check inputs. VOL and MASK must be 3-D matrices.');
end

if ~isequal(size(vol),size(mask))
    error('Check inputs. Size of mask does not equal size of img');
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

if exist('class_options','var')~=1 || isempty(class_options)
    class_options = {'raw','gray','gradient','haralick','gabor','laws','collage'};
end

if exist('ws_options','var')~=1 || isempty(ws_options)
    ws_options = 3:2:11;
end

ws1 = ws_options(1);
ws2 = ws_options(end);

%% RECOMMENDED: CROP VOLUME AND MASK (saves time and memory)!
[vol,mask] = boundingbox2(vol,mask,max(ws_options),'cropz','off');

%% 3D Feature Intensity Extraction

%--------------Raw Intensity -----------------%
if ~all(cellfun('isempty',strfind(class_options,'raw')))
    fprintf('\nExtracting raw intensity information:\n');
    
    rawints = {vol(mask==1)};
    rawnames = {'Raw intensity'};
end

%--------------Gray Level Statistics----------------%
if ~all(cellfun('isempty',strfind(class_options,'gray')))
    fprintf('\nExtracting Gray Level Statistics:\n');
   
    addpath([codepath,'grey']);   
    
    grayfeats = [];
    graynames = {};
    for ws = ws1:2:ws2
        fprintf('\tUsing a window size of %i...\n',ws)
        [gf,gfn] = grayfilts3_nomex(vol,ws);
        grayfeats=cat(4,grayfeats,gf);
        for s = 1:length(gfn)
            gfn{s} = strcat(gfn{s},' ws=',num2str(ws));
        end
        graynames(end+1:end+length(gfn))=gfn;
    end
    grayints = extract3DFeatIntensities(grayfeats, mask, 1);
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
        [gradfeats, gradnames] = gradfilts3(vol);
        gradints = extract3DFeatIntensities(gradfeats, mask, 1);
    end
    
    clear gradfeats;
    rmpath([codepath,'grey']);   
end

%--------------Haralick--------------%
if ~all(cellfun('isempty',strfind(class_options,'haralick')))
    fprintf('\nExtracting Haralick Features:\n');
    
    addpath([codepath,'haralick']);   
    haralickfun = @haralick3mex; % haralick3mexmt or, haralick3mex
    
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
    for ws = ws1:2:ws2
        harnames=cat(2,harnames,eval(['hfeatnames' int2str(ws)]));    
    end
      
    nbins = 128; % a power of 2
    bkgrd = -1; %b/c rescaling from 0 to 1
    dist = 1; %d = 0 (co-occurence with just self, d = 1 (co-occurence with just immediate neighbors), d = 2 (co-occurence with more distant neighbors)

    norm_vol = round(rescale_range(vol,0,nbins-1));
%     if length(unique(norm_vol))<nbins
%         error('Not enough gray levels present in VOL. Try resampling or not cropping.');
%     end
    hfeats = [];
    for ws = ws1:2:ws2
        fprintf('\tUsing a window size of %i...\n',ws)
        hfeats = cat(4,hfeats,haralickfun(norm_vol,nbins,ws,dist,bkgrd));
    end
    harints = extract3DFeatIntensities(hfeats, mask, 1);
    clear hfeats hfeatnames3 hfeatnames5 hfeatnames7 img ws;
    rmpath([codepath,'haralick']);   
end

%--------------Gabor--------------%
if ~all(cellfun('isempty',strfind(class_options,'gabor')))
    fprintf('\nExtracting Gabor Wavelets: ');
    
    addpath([codepath,'gabor']);
    load gaborbank3D; %gabor_data

    fcount = 0;
    for f=1:length(gabor_data) %to save memory, dont run all of Gabor!)
        if any(ismember(ws_options, gabor_data(f).ws))
            fcount = fcount+1;
            g_c = jconvn(vol,gabor_data(f).cos);
            g_s = jconvn(vol,gabor_data(f).sin);
            gabfeats(:,:,:,fcount) = sqrt(g_c.*g_c + g_s.*g_s);
            gabnames(fcount) = {gabor_data(f).title};
        end
        if mod(f,4) == 0, fprintf('.'); end
    end
    fprintf('\n');
    gabints = extract3DFeatIntensities(gabfeats, mask, 1);
    clear gabfeats g_c g_s gabor_data f;
    rmpath([codepath,'gabor']);
end

%--------------Laws--------------%
if ~all(cellfun('isempty',strfind(class_options,'laws')))
    
    fprintf('\nExtracting Laws Filters:\n');

    addpath([codepath,'laws']);
    addpath([codepath,'laws\mtimesx']);
    
    lawfeats = [];
    lawnames = {};
    
    if any(ws_options>5)
        fprintf('\tPlease note that Laws features cannot use a ws > 5 currently.\n');     
    end
    
    for ws = ws1:2:ws2
        if ws <=5
            fprintf('\tUsing a window size of %i:\n',ws)
            [lf, lfn] = lawsfilter3(vol,ws);
            lawfeats = cat(4,lawfeats,lf);
            lawnames=cat(2,lawnames,lfn);
        end
    end
    lawints = extract3DFeatIntensities(lawfeats, mask, 1);
    
    clear lf lfn lawfeats ws;
    
    rmpath([codepath,'laws']); 
    rmpath([codepath,'laws\mtimesx']);

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
                  C_O1_hfeatnames = {'CoLlAGe O=1 entropy ws=3' 'CoLlAGe O=1 energy ws=3' 'CoLlAGe O=1 inertia ws=3' 'CoLlAGe O=1 idm ws=3' 'CoLlAGe O=1 correlation ws=3' 'CoLlAGe O=1 info1 ws=3' ...
                                 'CoLlAGe O=1 info2 ws=3' 'CoLlAGe O=1 sum-av ws=3' 'CoLlAGe O=1 sum-var ws=3' 'CoLlAGe O=1 sum-ent ws=3' 'CoLlAGe O=1 diff-av ws=3' 'CoLlAGe O=1 diff-var ws=3' 'CoLlAGe O=1 diff-ent ws=3'};
                  C_O2_hfeatnames = {'CoLlAGe O=2 entropy ws=3' 'CoLlAGe O=2 energy ws=3' 'CoLlAGe O=2 inertia ws=3' 'CoLlAGe O=2 idm ws=3' 'CoLlAGe O=2 correlation ws=3' 'CoLlAGe O=2 info1 ws=3' ...
                                 'CoLlAGe O=2 info2 ws=3' 'CoLlAGe O=2 sum-av ws=3' 'CoLlAGe O=2 sum-var ws=3' 'CoLlAGe O=2 sum-ent ws=3' 'CoLlAGe O=2 diff-av ws=3' 'CoLlAGe O=2 diff-var ws=3' 'CoLlAGe O=2 diff-ent ws=3'};
              case 5
                  k = 2;
                  C_O1_hfeatnames = {'CoLlAGe O=1 entropy ws=5' 'CoLlAGe O=1 energy ws=5' 'CoLlAGe O=1 inertia ws=5' 'CoLlAGe O=1 idm ws=5' 'CoLlAGe O=1 correlation ws=5' 'CoLlAGe O=1 info1 ws=5' ...
                                 'CoLlAGe O=1 info2 ws=5' 'CoLlAGe O=1 sum-av ws=5' 'CoLlAGe O=1 sum-var ws=5' 'CoLlAGe O=1 sum-ent ws=5' 'CoLlAGe O=1 diff-av ws=5' 'CoLlAGe O=1 diff-var ws=5' 'CoLlAGe O=1 diff-ent ws=5'};
                  C_O2_hfeatnames = {'CoLlAGe O=2 entropy ws=5' 'CoLlAGe O=2 energy ws=5' 'CoLlAGe O=2 inertia ws=5' 'CoLlAGe O=2 idm ws=5' 'CoLlAGe O=2 correlation ws=5' 'CoLlAGe O=2 info1 ws=5' ...
                                 'CoLlAGe O=2 info2 ws=5' 'CoLlAGe O=2 sum-av ws=5' 'CoLlAGe O=2 sum-var ws=5' 'CoLlAGe O=2 sum-ent ws=5' 'CoLlAGe O=2 diff-av ws=5' 'CoLlAGe O=2 diff-var ws=5' 'CoLlAGe O=2 diff-ent ws=5'};
          end
          
          o1 = nan([size(mask),13]);
          o2 = nan([size(mask),13]);
          for j=2:size(mask,3)-1 %iterate over the patient mask
              if isempty(find(mask(:,:,j)==1,1))
                  continue
              end
              [o1(:,:,j,:),o2(:,:,j,:)]=compute_CoLlAGe3D(vol(:,:,j),vol(:,:,j +1), vol(:,:,j-1),mask(:,:,j),k);
          end
          
          collagefeats = cat(4,collagefeats,o1,o2);
          collagenames = cat(2,collagenames,C_O1_hfeatnames,C_O2_hfeatnames);
       end
   end
   collageints = extract3DFeatIntensities(collagefeats,mask,1); 
   rmpath([codepath,'haralick']); 
   rmpath([codepath,'collage']);
   clear collagefeats C_01_hfeatnames C_02_hfeatnames o1 o2 k
end 
fprintf('\n');

%% append featints and featnames   
startz = 0;
endz   = 0;
    
if ~all(cellfun('isempty',strfind(class_options,'raw')))
    startz = endz+1;
    endz   = endz+length(rawints);
    featints(startz:endz) = rawints;
    featnames(startz:endz) = rawnames;
    disp(['Added intensity for original volume!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gray')))
    startz = endz+1;
    endz   = endz+length(grayints);
    featints(startz:endz) = grayints; 
    featnames(startz:endz) = graynames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gray features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gradient')))
    startz = endz+1;
    endz   = endz+length(gradints);
    featints(startz:endz) = gradints; 
    featnames(startz:endz) = gradnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gradient features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'haralick')))
    startz = endz+1;
    endz   = endz+length(harints);
    featints(startz:endz) = harints; 
    featnames(startz:endz) = harnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' haralick features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gabor')))
    startz = endz+1;
    endz   = endz+length(gabints);
    featints(startz:endz) = gabints; 
    featnames(startz:endz) = gabnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' gabor features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'laws')))
    startz = endz+1;
    endz   = endz+length(lawints);
    featints(startz:endz) = lawints; 
    featnames(startz:endz) = lawnames;
    disp(['Added intensities for ', int2str(endz-startz+1),' laws features!'])
end


if ~all(cellfun('isempty',strfind(class_options,'collage')))
   startz = endz+1; 
   endz = endz+length(collageints);
   featints(startz:endz)= collageints;
   featnames(startz:endz) = collagenames;
   disp(['Added intensities for ', int2str(endz-startz+1),' collage features!'])
end


disp(['In total: extracted intensities for ', int2str(length(featints)),' texture features!']);       

%% Compute, Append 3D Feature Statistics featstats and statnames

if nargout <3 %user does not want to compute statistics of feature distributions
    return;
end

featstats = zeros(size(rawstats,1),size(rawstats,2)+size(graystats,2)+size(gradstats,2)+size(harstats,2)+size(gabstats,2)+size(lawstats,2));

startz = 0;
endz = 0;
    
if ~all(cellfun('isempty',strfind(class_options,'raw')))
    [rawstats,rawstatnames] = compute3DFeatStats(rawints,rawnames);
    startz = endz + 1;
    endz   = size(rawstats,2);
    featstats(:,startz:endz) = rawstats; 
    statnames(:,startz:endz) = rawstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' raw intensity!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gray')))
    [graystats,graystatnames] = compute3DFeatStats(grayints,graynames);
    startz = endz+1;
    endz   = endz+size(graystats,2);
    featstats(:,startz:endz) = graystats; 
    statnames(:,startz:endz) = graystatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gray features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gradient')))
    [gradstats,gradstatnames] = compute3DFeatStats(gradints,gradnames);
    startz = endz+1;
    endz   = endz+size(gradstats,2);
    featstats(:,startz:endz) = gradstats; 
    statnames(:,startz:endz) = gradstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gradient features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'haralick')))
    [harstats,harstatnames] = compute3DFeatStats(harints,harnames);
    startz = endz+1;
    endz   = endz+size(harstats,2);
    featstats(:,startz:endz) = harstats; 
    statnames(:,startz:endz) = harstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' haralick features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'gabor')))
    [gabstats,gabstatnames] = compute3DFeatStats(gabints,gabnames);
    startz = endz+1;
    endz   = endz+size(gabstats,2);
    featstats(:,startz:endz) = gabstats; 
    statnames(:,startz:endz) = gabstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gabor features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'laws')))
    [lawstats,lawstatnames] = compute3DFeatStats(lawints,lawnames);
    startz = endz+1;
    endz = endz+size(lawstats,2);
    featstats(:,startz:endz) = lawstats; 
    statnames(:,startz:endz) = lawstatnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' laws features!'])
end

if ~all(cellfun('isempty',strfind(class_options,'collage'))) 
    [collagestats,collagestatnames]= compute3DFeatStats(collageints,collagenames);
    startz = endz+1; 
    endz = endz+size(collagestats,2); 
    featstats(:,startz:endz) = collagestats;
    statnames(:,startz:endz)=collagestatnames; 
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' collage features!'])
end

disp(['In total: computed ' int2str(size(featstats,1)) ' statistics for ', int2str(size(featstats,2)),' texture features!']); 
    
clear rawints grayints gradints harints gabints lawints;
clear rawnames graynames gradnames hnames gabnames lawnames;
clear rawstats graystats gradstats harstats gabstats lawstats;
clear rawstatnames graystatnames gradstatnames harstatnames gabstatnames lawstatnames;        
end

%% extract3DFeatIntensities subfunction
% Jacob Antunes
% August 04, 2015

function varargout = extract3DFeatIntensities(featVol, mask, vals)
%INPUTS
%featVol = 4-D volume of texture features
%mask = 3-D annotation volume
%vals = vector of mask values to find intensities within

%OUTPUT
%varargout = cell array(s) of intensities. Number based on length of vals

    if isempty(featVol) || isempty(mask) || isempty(vals)
        varargout = {};
        return;
    end
    
    fprintf('\n\tExtracting intensities\n');
    
    ints = cell(size(featVol,4),length(vals)); %cell array of intensities
    
    for i = 1:size(featVol,4) % loop through all features
        f = featVol(:,:,:,i); %single feature volume
        for j = 1:length(vals) % loop through all mask values
            ints{i,j} = f(mask == vals(j));
        end
    end
    
    for k = 1:length(vals)
        varargout{k} = ints(:,k);
    end
    
end

%% compute3DFeatStats subfunction
% Jacob Antunes
% April 11, 2017

function [featstats,statnames]=compute3DFeatStats(featints,featnames)
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
                     var(temp1); ...
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
%          statnames{1,i} = strcat('mean-',featnames{i});
%          statnames{3,i} = strcat('range-',featnames{i});
%         statnames{7,i} = strcat('harmean-',featnames{i});
%         statnames{4,i} = strcat('std-',featnames{i});
%          statnames{1,i} = strcat('mean-',featnames{i});
%         statnames{9,i} = strcat('mad-',featnames{i});
%         statnames{10,i} = strcat('iqr-',featnames{i});
   end
end

end
       

 
 
