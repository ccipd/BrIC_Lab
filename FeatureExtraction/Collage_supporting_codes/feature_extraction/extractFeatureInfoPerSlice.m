function [slices_featinfo, valid_slices] = extractFeatureInfoPerSlice(vol, volmask, class_options, ws_options)
% Extract and store features from all slices individually

%INPUTS:
    %vol = 3D matrix of multiple 2D slices (double format)
    %volmask = corresponding 3D matrix of multiple 2D labels; same size as vol (ROI = 1)
    %class_options = (OPTIONAL) cell array of strings corresponding to desired feature classes:
    %                 DEFAULT: class_options = {'raw','gray','gradient','haralick','gabor','laws'};
    %ws_options = (OPTIONAL) array of integers corresponding to desired window levels:
    %                 DEFAULT: ws_options = [3, 5, 7, 9, 11]
   
%OUTPUTS:f
    %slices_featinfo = struct array containing the following field for each slice:
    %          .featints
    %          .featnames
    %          .featstats
    %          .statnames   [see extract2DFeatureInfo.m for descriptions]
    %valid_slices = row vector containing indices of non-empty mask slices

    %see supporting function(s):
    % function extract2DFeatureInfo.m
    
    
% @ Jacob Antunes, 2018
% Last Updated: 04-24-2018

%% Check inputs

if exist('extract2DFeatureInfo.m','file')~=2
    error('Get extract2DFeatureInfo.m into your path!');
end

if isempty(vol) || isempty(volmask)
    error('Check inputs. Either VOL or VOLMASK is an empty matrix.');
end

if ~isequal(size(vol),size(volmask))
    error('Check inputs. Size of VOLMASK does not equal size of VOL');
end

if isempty(find(volmask>0,1))
    error('VOLMASK is empty. Nothing to extract.');
end

if numel(size(vol))~=3
    error('extractFeatureInfoPerSlice.m is designed for 3D inputs (multiple 2D slices) only.');
end

%% Initialization

if exist('class_options','var')~=1 || isempty(class_options)
    class_options = {'raw','gray','gradient','haralick','gabor','laws','collage'};
end

if exist('ws_options','var')~=1 || isempty(ws_options)
    ws_options = 3:2:11;
end

ws1 = ws_options(1);
ws2 = ws_options(end);
ws_options = ws1:2:ws2; %safety check that only odd window sizes are used

valid_slices = [];

for i = 1:size(vol,3)

    fprintf('slice #%i of %i:\n',i,size(vol,3));pause(0.2);

    img = vol(:,:,i);
    mask = volmask(:,:,i);
    
    if isempty(find(mask==1,1))
        fprintf('No mask on this slice. Continuing.\n');
        slices_featinfo(i).featints = [];
        slices_featinfo(i).featnames = [];
        slices_featinfo(i).featstats = [];
        slices_featinfo(i).statnames = [];
        slices_featinfo(i).img = [];
        slices_featinfo(i).mask = [];
        continue;
    end
    

    valid_slices = [valid_slices i];

    [featints,featnames,featstats,statnames] = extract2DFeatureInfo(img,mask,class_options,ws_options);        

    slices_featinfo(i).featints = featints;
    slices_featinfo(i).featnames = featnames;
    slices_featinfo(i).featstats = featstats;
    slices_featinfo(i).statnames = statnames;
    slices_featinfo(i).img = img;
    slices_featinfo(i).mask = mask;
        
%     %checkpoint (comment/uncomment):
%     figure;vv(img,mask)
%     figure;feature_map(img,mask,featints{1});
    fprintf('\n');
end

