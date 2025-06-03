subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04','MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10'};%'SIC02',


xdistance = 30;

columns = [4 4 4];

IDs_toinclude = [9.5 10.8 16.5 11.4 1.5 10 11 17 8 5]; 

thresholdarray = [.05 .1 .15 .2];%[.05 : .05 : .3];

minsize = 40;

cd /data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/temp/gephi/

for subnum = [6 9]
    
    subname = subnames{subnum};
    disp(subname)
    
    CONsubnets = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/temp/' subname '_con_subnetworks_autodetected.dtseries.nii']);
    
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
        
        datafolder = [basedir subname '/bold1_222/'];
        
        scanstouse_inds = 1:12; %pre-cast
        
        scanlist = scanlist(scanstouse_inds);
        
        for scanindnum = 1:length(scanlist)
            scanind = scanstouse_inds(scanindnum);
            scanname = scanlist{scanind};
            ciftifile = dir([datafolder scanname '*surfsmooth2.55_volsmooth2.dtseries.nii']);
            
            data = ft_read_cifti_mod([datafolder ciftifile(1).name]);
            data.data = data.data(:,logical(tmasks(scanind,:)));
            
            if scanindnum==1
                alldata = data;
                
            else
                
                alldata.data = [alldata.data data.data];
                
            end
            
            clear data
        end
        
        data = alldata;
        clear alldata;
        
        sessions = [];
        tmask_concat = [];
        inds_withintmaskedconcat = [];
        prevind = 0;
        for s = 1:numel(scanlist)
            tmask = tmasks(s,:)';
            sessions = [sessions ; repmat(s,size(tmask,1),1)];
            tmask_concat = [tmask_concat ; tmask];
            inds_withintmasked = zeros(length(tmask),1);
            inds_withintmasked(tmask) = [(prevind+1) : (prevind + nnz(tmask))];
            inds_withintmaskedconcat = [inds_withintmaskedconcat ; inds_withintmasked];
            prevind = prevind + nnz(tmask);
        end
        
        tmask_concat = logical(tmask_concat);
        
        
        
    elseif strcmp(subname(1:2),'ME')
        column = columns(2);
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        %motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_oneID_CS.dscalar.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        %motorspots.data(59413:end,:) = 0;
        networks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        subnetworks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/rawassn_minsize10_regularized.dscalar.nii']);
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        atlasT1file = [basedir '/anat/T1w/T1w_acpc.nii.gz'];
        
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        
    elseif strcmp(subname(1:3),'MSC')
        column = columns(3);
        infomapdir = ['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/'];
        
        data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' subname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
    
        dmatname = ['/data/nil-bluearc/GMT/Scott/MSC/distmat_surf_geodesic_vol_euclidean_MNInlwarp_xhem1000_uint8.mat'];
    
        networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
    
        subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized.dtseries.nii']);
        
    end
    
    
    %
    
    
    %
    networks.data(59413:end) = 0;
    CONsubnets.data(59413:end,:) = 0;
    networks.data(networks.data==9) = 0;
    networks.data(CONsubnets.data>0) = CONsubnets.data(CONsubnets.data>0);
    
    subnetworks.data = subnetworks.data(:,column);
    subnetworks.data(59413:end) = 0;
    
    
%     IDs = unique(CONsubnets.data); IDs(IDs==0) = [];
%     for IDnum = length(IDs) : -1 : 1
%         if ~any(abs(IDs_toinclude - IDs(IDnum))<.01)
%             IDs(IDnum) = [];
%         end
%     end
    
subnetworkIDs = unique(subnetworks.data);
subnetworkIDs(subnetworkIDs==0) = [];

    allclusters = [];
    allclusterIDs = [];
    for IDnum = 1:length(subnetworkIDs)
        ID = subnetworkIDs(IDnum);
        networkID = mode(networks.data(subnetworks.data==ID));
        if any(abs(IDs_toinclude-networkID)<.01)
            clusters_thisID = cifti_cluster(subnetworks,ID,ID,minsize);
        
            allclusters = [allclusters clusters_thisID];
            allclusterIDs = [allclusterIDs repmat(networkID,1,size(clusters_thisID,2))];
        end
    end
    
    
    cluster_sizes = sum(allclusters,1);
    
    allclusters = logical(allclusters);
    
    
    clustertcs = zeros(size(data.data,2),size(allclusters,2));
    for c = 1:size(allclusters,2)
        clustertcs(:,c) = mean(data.data(allclusters(:,c),:),1);
    end
    rmat = paircorr_mod(clustertcs);
    
    
    distances = smartload(dmatname);
    centroids = zeros(size(rmat,1),1);
    for c = 1:size(allclusters,2)
        clusterinds = find(allclusters(:,c));
        [~,mini] = min(sum(distances(clusterinds,clusterinds),2));
        centroids(c) = clusterinds(mini);
    end
    dmat = distances(centroids,centroids);
    
    rmat(dmat < xdistance) = 0;
    
    %%
    
    
    for t = 1:length(thresholdarray)
        
        sorted_corrvals = sort(rmat(:),'descend');
        thresh_val = sorted_corrvals(floor(length(sorted_corrvals).*thresholdarray(t)));
        
        graph = rmat >= thresh_val;
        
        Make_gephi_files_func_v3(graph,allclusterIDs',{'NodeID'},[subname '_CONsubnets_networks_thr' num2str(thresholdarray(t))],[cluster_sizes'],{'Size'})
        
    end
    
end
    
