subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04','MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10'};%'SIC02',

cd /data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/

xdistance = 30;

columns = [4 4 4];

networkIDs = [9];% CON

noisemap = ft_read_cifti_mod('/data/nil-bluearc/GMT/Evan/Temp/CON_insula_noise_map.dtseries.nii');

groupnetworks = ft_read_cifti_mod('/data/nil-bluearc/GMT/Evan/Atlases/WashUNetworks/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');

colors = [10.8, 16.5, 9.5, 11.4, 2.5];

dilation = 0;

inMNI = [0 0 0 1 1 1 1 0 0 0 0 0 0 0 0];

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
        
        column = columns(1);
        
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Laumann/MSC/MSC02/T1/MSC02_mpr_debias_avgT_111_t88.nii.gz'];
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/T1/' subname '_mpr_debias_avgT_111_t88.nii.gz'];
            
            if strcmp(subname,'SIC02')
                dmatname =  ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC06/fsaverage_LR32k/MSC06_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            else
                dmatname = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/bold1_222/cifti_distances/' subname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];
            end
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
         
        networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        subnetworks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized.dtseries.nii']);
        
%         datafolder = [basedir subname '/bold1_222/'];
%         
%         scanstouse_inds = 1:12; %pre-cast
%         
%         scanlist = scanlist(scanstouse_inds);
%         
%         for scanindnum = 1:length(scanlist)
%             scanind = scanstouse_inds(scanindnum);
%             scanname = scanlist{scanind};
%             ciftifile = dir([datafolder scanname '*surfsmooth2.55_volsmooth2.dtseries.nii']);
%             
%             data = ft_read_cifti_mod([datafolder ciftifile(1).name]);
%             data.data = data.data(:,logical(tmasks(scanind,:)));
%             
%             if scanindnum==1
%                 alldata = data;
%                 
%             else
%                 
%                 alldata.data = [alldata.data data.data];
%                 
%             end
%             
%             clear data
%         end
%         
%         data = alldata;
%         clear alldata;
%         
%         sessions = [];
%         tmask_concat = [];
%         inds_withintmaskedconcat = [];
%         prevind = 0;
%         for s = 1:numel(scanlist)
%             tmask = tmasks(s,:)';
%             sessions = [sessions ; repmat(s,size(tmask,1),1)];
%             tmask_concat = [tmask_concat ; tmask];
%             inds_withintmasked = zeros(length(tmask),1);
%             inds_withintmasked(tmask) = [(prevind+1) : (prevind + nnz(tmask))];
%             inds_withintmaskedconcat = [inds_withintmaskedconcat ; inds_withintmasked];
%             prevind = prevind + nnz(tmask);
%         end
%         
%         tmask_concat = logical(tmask_concat);
        
        
        
    elseif strcmp(subname(1:2),'ME')
        column = columns(2);
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        
%         data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
%         
%         tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
%         tmask_concat = logical(tmask_concat);
%         sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        
        networks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        subnetworks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/rawassn_minsize10_regularized.dscalar.nii']);
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        atlasT1file = [basedir '/anat/T1w/T1w_acpc.nii.gz'];
        
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        
    elseif strcmp(subname(1:3),'MSC')
        column = columns(3);
        infomapdir = ['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/'];
        
        
        %data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' subname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
    
        dmatname = ['/data/nil-bluearc/GMT/Scott/MSC/distmat_surf_geodesic_vol_euclidean_MNInlwarp_xhem1000_uint8.mat'];
    
        networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
    
        subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized.dtseries.nii']);
        
    end
    
    
    %% subnets into common subcort space
    
    
    subnetworks.data(logical(noisemap.data(1:59412)),:) = 0;
    
    if subnum==1
        
        subnetworks_allsubs_commonspace = subnetworks;
        subnetworks_allsubs_commonspace.data = zeros(size(subnetworks_allsubs_commonspace.data,1),length(subnames));
        subnetworks_allsubs_commonspace.data(:,subnum) = subnetworks.data(:,column);
        
        
        ncortverts = nnz(subnetworks_allsubs_commonspace.brainstructure==1) + nnz(subnetworks_allsubs_commonspace.brainstructure==2);
        subcort_starting_ind = max(find(subnetworks_allsubs_commonspace.brainstructure==2)) + 1;
        
        subcort_coords = subnetworks_allsubs_commonspace.pos(subcort_starting_ind:end,:);
        
        subnetworks_sorted_bynetwork = cell(1,length(networkIDs));
        
        networks_allsubs = networks.data;
        networks_allsubs = networks_allsubs(1:ncortverts);
        
    else
        
        subnetworks_allsubs_commonspace.data(1:59412,subnum) = subnetworks.data(1:59412,column);
        
        if logical(inMNI(subnum))
            subnetworks.pos = mni2tal(subnetworks.pos);
        end
        
        for vox = 1:length(subcort_coords)
            D = pdist2(subcort_coords(vox,:),subnetworks.pos(subcort_starting_ind:end,:));
            zeroind = find(D<=1);
            if ~isempty(zeroind)
                zeroind = zeroind(D(zeroind)==min(D(zeroind)));
                subnetworks_allsubs_commonspace.data(ncortverts+vox,subnum) = subnetworks.data(ncortverts+zeroind,column);
            end
        end
        
        networks_allsubs(:,subnum) = networks.data(1:ncortverts);
         
    end
    
end
    


%%
for IDnum = 1:length(networkIDs)
    
    subnets_allsubs_thisnetwork = subnetworks_allsubs_commonspace;
    subnets_allsubs_thisnetwork.data = [];
    subnets_allsubs_thisnetwork_subtracker = [];
    
    for subnum = 1:length(subnames)
        thisnetwork_thissub = networks_allsubs(:,subnum) == networkIDs(IDnum);
        
        subnetIDs = unique(subnetworks_allsubs_commonspace.data(:,subnum));
        subnetIDs(subnetIDs==0) = [];
        
        for s = 1:length(subnetIDs)
            if nnz((subnetworks_allsubs_commonspace.data(1:ncortverts,subnum)==subnetIDs(s)) & thisnetwork_thissub) ./ nnz(subnetworks_allsubs_commonspace.data(1:ncortverts,subnum)==subnetIDs(s)) > .5
                subnets_allsubs_thisnetwork.data(:,end+1) = (subnetworks_allsubs_commonspace.data(:,subnum)==subnetIDs(s));
                subnets_allsubs_thisnetwork_subtracker(1,end+1) = subnum;
            end
        end
        
        if nnz(subnets_allsubs_thisnetwork_subtracker==subnum) < 5
            subnets_allsubs_thisnetwork.data(:,subnets_allsubs_thisnetwork_subtracker==subnum) = [];
            subnets_allsubs_thisnetwork_subtracker(subnets_allsubs_thisnetwork_subtracker==subnum) = [];
            
            thisnetwork_thissub = groupnetworks.data(1:ncortverts) == networkIDs(IDnum);
            
            subnetIDs = unique(subnetworks_allsubs_commonspace.data(:,subnum));
            subnetIDs(subnetIDs==0) = [];
            
            for s = 1:length(subnetIDs)
                if nnz((subnetworks_allsubs_commonspace.data(1:ncortverts,subnum)==subnetIDs(s)) & thisnetwork_thissub) ./ nnz(subnetworks_allsubs_commonspace.data(1:ncortverts,subnum)==subnetIDs(s)) > .5
                    subnets_allsubs_thisnetwork.data(:,end+1) = (subnetworks_allsubs_commonspace.data(:,subnum)==subnetIDs(s));
                    subnets_allsubs_thisnetwork_subtracker(1,end+1) = subnum;
                end
            end
        end
        
    end
    
    neighbors = cifti_neighbors(subnets_allsubs_thisnetwork);
    neighbors_nonan = neighbors; neighbors_nonan(isnan(neighbors_nonan)) = 1;
    for s = 1:size(subnets_allsubs_thisnetwork.data,2)
        for d = 1:dilation
            neighbors_vals = single(subnets_allsubs_thisnetwork.data(neighbors_nonan,s));
            neighbors_vals(isnan(neighbors)) = 1000;
            expansion_inds = neighbors_vals(:,1)==0 & any(neighbors_vals(:,2:end)==1,2);
            subnets_allsubs_thisnetwork.data(expansion_inds,s) = true;
        end
    end
    
    
    % set seed for
    
    % preallocate similarity matrix;
    SimMat = zeros(size(subnets_allsubs_thisnetwork.data,2));

    % sweep the rows
    for i = 1:size(SimMat,1)

        % sweep the columns
        for ii = (i+1):size(SimMat,2)

            % calculate spatial overlap 
            SimMat(i,ii) = jaccard(subnets_allsubs_thisnetwork.data(:,i),subnets_allsubs_thisnetwork.data(:,ii));
            SimMat(ii,i) = SimMat(i,ii);

        end

    end
    %%
    % reproducibility
    rng(1,'twister')
    
    Ci  = 1:size(SimMat,1); % initial community affiliations
    Q0 = -1; Q = 0; % initialize modularity values
    while Q - Q0 > 1e-5 % while modularity increases
        Q0 = Q; % perform community detection
        [Ci,Q] = community_louvain(SimMat,.75,Ci);
    end
    
    % initial unique clusters;
    uCi = unique(nonzeros(Ci));
    
    % preallocate;
    nSubjects = zeros(length(uCi),1);
    
    % sweep through the clusters;
    for i = 1:length(uCi)
        
        % calculate the number of subjects
        nSubjects(i) = length(unique(subnets_allsubs_thisnetwork_subtracker(Ci==uCi(i))));
        
    end
    
    % remove subnetworks not present
    % in at least 1/2 of individuals;
    uCi(nSubjects <= (length(subnames) * 2/3)) = []; 
    nSubjects(nSubjects <= (length(subnames) * 2/3)) = []; 
    
    % remove bad clusters;
    Ci(~ismember(Ci,uCi))=0;
    
    % sort from 
    % most frequent 
    % to less frequent;
    [nSubjects_sorted,Clusters_sortind] = sort(nSubjects,'Descend');
    
    % preallocate;
    Ci_final = zeros(length(Ci),1);
    
    % serial labels;
    for i = 1:length(uCi)
        Ci_final(Ci==uCi(Clusters_sortind(i))) = i;
    end
    
    %Ci_final(Ci_final==5) = 2;
    %Ci_final(Ci_final==6) = 0;
    
    % update unique clusters;
    uCi = unique(nonzeros(Ci_final));
    
    % if subnetworks
    % are networks;
    if ~isempty(uCi)
        
        out = subnets_allsubs_thisnetwork;
        out.data = zeros(size(out.data,1),length(uCi));
        out.dimord = 'pos_time';
        
        % sweep all the
        % unique communities
        for i = 1:length(uCi)
            out.data(:,i) = sum(subnets_allsubs_thisnetwork.data(:,Ci_final==uCi(i)),2) / length(subnames);
        end
        
        % write out the subnetworks;
        ft_write_cifti_mod(['Subnetworks_ofnetwork' num2str(networkIDs(IDnum)) '_density'],out);
        
        
        [maxval,maxi] = max(out.data,[],2);
        vertexwisecolors = colors(maxi);
        
        out.data = zeros(size(out.data,1),1);
        out.data = vertexwisecolors(:) .* (maxval(:) > (2 / length(subnames)));
        ft_write_cifti_mod(['Subnetworks_ofnetwork' num2str(networkIDs(IDnum)) '_winner'],out);
        set_cifti_powercolors(['Subnetworks_ofnetwork' num2str(networkIDs(IDnum)) '_winner.dtseries.nii']);
        
        for subnum = 1:length(subnames)
            out = subnets_allsubs_thisnetwork;
            out.data = zeros(size(out.data,1),1);
            out.dimord = 'pos_time';
            
            for i = 1:length(uCi)
                out.data(any(subnets_allsubs_thisnetwork.data(:,Ci_final==uCi(i) & subnets_allsubs_thisnetwork_subtracker'==subnum),2)) = colors(i);
            end
            
            ft_write_cifti_mod([subnames{subnum} '_con_subnetworks_autodetected.dtseries.nii'],out)
            set_cifti_powercolors([subnames{subnum} '_con_subnetworks_autodetected.dtseries.nii'])
            
        end
            
        
    else
        
        % 
        Ci_final = ones(length(Ci),1);
        
    end
    
end
        
  