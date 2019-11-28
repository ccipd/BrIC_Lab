%% extract3DFeatures
% Jacob Antunes
% Sept 30, 2016
%This code has been released by CCIPD, Case Western Reserve University
%(Cleveland OH)

function [featstats,statnames] = extract3DFeatures(vol,mask)

    %INPUTS:
    %vol = 3D vol (double format)
    %mask = 3D label vol(ROI = 1)
    
    %OUTPUTS:
    %featstats =3D matrix of feature statistics from ROI within 4D feature volume representations of vol
    %statnames = 3D cell array of feature statistic names
    
    %% FIRST: Create Isotropic Volume respresentations using .sh files! RECOMMENDED: CROP VOLUME AND MASK!
    
    repopath = '~\ccipd\'; %don't forget the "/" at the end!

    %% Feature Extraction

%     %--------------Haralick--------------%
    addpath([repopath,'Haralick']);   
    fprintf('\nExtracting Haralick:\n');
    hfeatnames5 = {'entropy ws=5' 'energy ws=5' 'inertia ws=5' 'idm ws=5' 'correlation ws=5' 'info1 ws=5' ...
    'info2 ws=5' 'sum_av ws=5' 'sum_var ws=5' 'sum_ent ws=5' 'diff_av ws=5' 'diff_var ws=5' 'diff_ent ws=5'};
    hfeatnames=cat(2,hfeatnames5);

    img = round(rescale_range(vol,0,127));
    hfeats = [];
    for ws = 5
        fprintf('\tUsing a window size of %i...\n',ws)
        hfeats = cat(4,hfeats,haralick3mexmt(img,128,ws,2,0));
    end
    [hstats,hnames] = compute_FeatStats(hfeats,mask,hfeatnames);
    clear hfeats hfeatnames5 img ws;
    rmpath([repopath,'haralick']);   
    
    %--------------Gabor--------------%
    fprintf('\nExtracting Gabor: ');
    
    addpath([repopath,'Gabor']);
    load gaborbank3D; %gaborstruct

    for f=1:length(gaborstruct) %to save memory, dont run all of Gabor!)
        g_c = jconvn(vol,gaborstruct(f).cos);
        g_s = jconvn(vol,gaborstruct(f).sin);
        gabfeats(:,:,:,f) = sqrt(g_c.*g_c + g_s.*g_s);
        % collapsing via l_infinity norm, results in 9 features, comment out if NOT collapsing
        %if mod(f,norient)==0, gabfeats(:,:,:,f\n_xzorient) = max(Gabtemp(:,:,:,1:f),[],3); end
        gabfeatnames(f) = {gaborstruct(f).title};
        if mod(f,4) == 0, fprintf('.'); end
    end
    fprintf('\n');
    % if you don't want to collapse i.e. you want ALL Gabor features, uncomment the below line
    [gabstats,gabnames] = compute_FeatStats(gabfeats,mask,gabfeatnames);
    clear gabfeats gabfeatnames g_c g_s gaborstruct f;
    rmpath([repopath,'gabor']);

%     %--------------Laws--------------%
    fprintf('\nExtracting Laws:\n');

    addpath([repopath,'LawsEnergy']);
    addpath([repopath,'LawsEnergy\mtimesx']);
    lawfeats = [];
    lawfeatnames = {};
    fprintf('\tUsing a window size of 5\n')
    [lf, lfn] = lawsfilter3(vol,5);
    lawfeats = cat(4,lawfeats,lf);
    lawfeatnames=cat(2,lawfeatnames,lfn);
    [lawstats,lawnames] = compute_FeatStats(lawfeats,mask,lawfeatnames);
    clear lf lfn lawfeats lawfeatnames ws;
    rmpath([repopath,'LawsEnergy']); 
    rmpath([repopath,'LawsEnergy\mtimesx']);

    fprintf('\n');
    
%% append featstats and statnames
     featstats = zeros(size(hstats,1), size(hstats,2)+size(gabstats,2)+size(lawstats,2));
     statnames = cell(size(featstats));    
     startz = 1;
     endz=size(hstats,2);
     featstats(:,startz:endz) = hstats; 
     statnames(:,startz:endz) = hnames;
     disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' haralick features!'])

    startz = endz+1;
    endz   = endz+size(gabstats,2);
    featstats(:,startz:endz) = gabstats; 
    statnames(:,startz:endz) = gabnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' gabor features!'])

    startz = endz+1;
    endz = endz+size(lawstats,2);
    featstats(:,startz:endz) = lawstats; 
    statnames(:,startz:endz) = lawnames;
    disp(['Added ',int2str(size(featstats,1)),' stats for ', int2str(endz-startz+1),' laws features!'])
   
     disp(['In total: ', int2str(size(featstats,1)),' stats for ', int2str(size(featstats,2)),' texture features!']); 
    
    clear hstats gabstats lawstats;
    clear hnames gabnames lawnames;
        
end

 
