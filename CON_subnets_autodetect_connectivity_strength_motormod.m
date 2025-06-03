

subnames = {'MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10','ME01','ME02','ME03','ME04','SIC01','SIC02','SIC03'};

brainstructures_touse_groups = {[1 2 9 10 13:18],[1 2 10 11 16:21]};

subnetworkIDs = [10.8, 16.5, 9.5, 11.4];
subnetworknames = {'Feedback', 'Decision','Action','ParsMarg'};
networkIDs = [1 1.5 2 3 5 7 8 10 11 12 15 16 17];
networknames = {'DMN','SCAN','Vis','FPN','DAN','Lang','Sal','SMH','SMM','Aud','PMN','CAN','SMF'};

%displayorder = [7 16 2 17 10 11 1.5 5 8 15 1 3];
displayorder = [6 12 3 13 8 9 2 5 7 11 1 4];
%%
allFC = zeros(length(subnames),length(subnetworkIDs),length(networkIDs));
%loading data
for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    CONsubnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/temp/' subname '_con_subnetworks_autodetected.dtseries.nii']);
    %['/data/nil-bluearc/GMT/Lina/subnetworks/' subname '_con_subnetworks.dtseries.nii']);
    
    if strcmp(subname(1:3),'SIC')
        
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        brainstructures_touse = brainstructures_touse_groups{1};
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
                        
        else
            
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');
            
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
        
        networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        allsubnetworks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized.dtseries.nii']); allsubnetworks.data = allsubnetworks.data(:,4);
        
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
        
        
        
        
        
    elseif strcmp(subname(1:2),'ME')
        
        brainstructures_touse = brainstructures_touse_groups{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
       
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        
        networks = ft_read_cifti_mod([infomapdir '/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        allsubnetworks = ft_read_cifti_mod([infomapdir '/rawassn_minsize10_regularized.dscalar.nii']); allsubnetworks.data = allsubnetworks.data(:,4);
        
        fslrdir = [basedir 'anat/MNINonLinear/fsaverage_LR32k/'];
        
    elseif strcmp(subname(1:3),'MSC')
        data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' subname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
        networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        allsubnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized.dtseries.nii']); allsubnetworks.data = allsubnetworks.data(:,4);

    end
    %starting looping the subnetworks
    
    
    con_subnetwork_tcs = zeros(length(subnetworkIDs),size(data.data,2));
    
    
    
    for subnetwork = 1:length(subnetworkIDs)
        if any(abs(CONsubnetworks.data-subnetworkIDs(subnetwork))<.01)
            con_subnetwork_tcs(subnetwork,:) = mean(data.data(abs(CONsubnetworks.data(1:59412)-subnetworkIDs(subnetwork))<.01,:),1);
        end
        
    end
    network_tcs = zeros(length(networkIDs),size(data.data,2));
    for network = 1:length(networkIDs)
        network_tcs(network,:) = mean(data.data(networks.data==networkIDs(network),:),1);
    end
    
    FC = FisherTransform(paircorr_mod(con_subnetwork_tcs',network_tcs'));    
    allFC(subnum,:,:) = FC;
end

save('/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/temp/network_connectivity_autodetected.mat','allFC');

%%

subnets = repmat(subnetworkIDs,size(allFC,1),1,size(allFC,3));
subs = repmat([1:size(allFC,1)]',1,size(allFC,2),size(allFC,3));
nets = zeros(size(allFC));
for n = 1:size(allFC,3)
    nets(:,:,n) = n;
end
modelterms = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
[P,T,STATS,TERMS]=anovan(allFC(:),{subnets(:),subs(:),nets(:)},'random',2,'model',modelterms);
%%

for n = 1:length(networkIDs)
    thisFC = allFC(:,:,n);
    thissubnets = repmat(subnetworkIDs,size(allFC,1),1);
    thesesubs = repmat([1:size(allFC,1)],1,size(allFC,2));
    [P,T,STATS,TERMS]=anovan(thisFC(:),{thissubnets(:),thesesubs(:)},'random',2,'display','off');
    disp(['Network ' num2str(networkIDs(n)) ': F(' num2str(T{2,3}) ',' num2str(T{4,3}) ') = ' num2str(T{2,6}) ', p = ' num2str(P(1)) ', corrected p = ' num2str(P(1) .* length(networkIDs))])
    
end

%%
ps = [];
for n = 1:length(networkIDs)
    string = [];
    for s = 1:length(subnetworkIDs)
        for s2 = (s+1):length(subnetworkIDs)
            [H,P,CI,STATS] = ttest(allFC(:,s,n),allFC(:,s2,n));
            %disp(['Network ' num2str(networkIDs(n)) '; subnetwork ' subnetworknames{s} ' vs ' subnetworknames{s2} ': T(' num2str(STATS.df) ')=' num2str(STATS.tstat) '; p=' num2str(P)])
            string = [string 'T(' num2str(STATS.df) ')=' sprintf('%02.2f',STATS.tstat) ';p=' sprintf('%02.3f',P) ', '];
            ps(end+1) = P;
        end
    end
    disp(string)
end

%%
make_network_polar_plot_general(squeeze(nanmean(allFC(:,:,displayorder),1))',subnetworkIDs, repmat({' '},1,length(displayorder)));%,network_names(displayorder))

%% DMN vs others
for s = 1:length(subnetworkIDs)
    for n = 2:length(networkIDs)
        [H,P,CI,STATS] = ttest(allFC(:,s,1),allFC(:,s,n));
        disp(['subnetwork ' subnetworknames{s} ': DMN vs Network ' num2str(networkIDs(n)) ': T(' num2str(STATS.df) ')=' num2str(STATS.tstat) '; p=' num2str(P)])
    end
end
