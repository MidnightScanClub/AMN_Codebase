

subnames = {'SIC01','SIC02','SIC03','MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10','ME01','ME02','ME03','ME04'};
trs = [1.1 1.355 2.2];
subnetworkIDs = [10.8, 16.5, 9.5, 11.4];
networkIDs = [1.5 8 10];%[1.5 5 8 10 11 17];

interval_times = [.03];


for subnum = [1:length(subnames)]
    
    subname = subnames{subnum};
    disp(subname)
    
    subnets = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/' subname '_con_subnetworks_autodetected.dtseries.nii']);
    
    
     
    
    ciftidata = cell(0,1);
    
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        tr = trs(1);
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        
        networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        %brainstructures_touse = brainstructures_touse_groups{1};
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            
            
        end
        
        
        datafolder = [basedir subname '/bold1_222/'];
        
        scanstouse_inds = 1:12; %pre-cast
        
        scanlist = scanlist(scanstouse_inds);
        
        for scanindnum = 1:length(scanlist)
            scanind = scanstouse_inds(scanindnum);
            scanname = scanlist{scanind};
            ciftifile = dir([datafolder scanname '*surfsmooth2.55_volsmooth2.dtseries.nii']);
            
            data = ft_read_cifti_mod([datafolder ciftifile(1).name]);
            
            data.data = interp_boldts_linear(data.data, logical(tmasks(scanind,:))');
            ciftidata(end+1) = {data};
            
            
            
            
        end
        
        
        
        
        
    elseif strcmp(subname(1:2),'ME')
        
        
        tr = trs(2);
        %brainstructures_touse = brainstructures_touse_groups{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        %cd([infomapdir])
        
        data_orig = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        data = data_orig;
        data.data = zeros(size(data.data,1),length(sessions));
        data.data(:,logical(tmask_concat)) = data_orig.data;
        data.data = interp_boldts_linear(data.data, tmask_concat);
        
        sessnums = unique(sessions);
        for s = 1:length(sessnums)
            thisdata = data;
            thisdata.data = data.data(:,sessions==sessnums(s));
            ciftidata(s) = {thisdata};
        end
        
        networks = ft_read_cifti_mod([infomapdir '/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        
        
        
    elseif strcmp(subname(1:3),'MSC')
        tr = trs(3);
        
        subfolder = ['/data/nil-bluearc/GMT/Laumann/MSC/' subname '/Functionals/'];
        ciftifiles = dir([subfolder 'FCPROCESS_SCRUBBED_UWRPMEAN_MSM_V2/cifti_timeseries_normalwall_native_freesurf_222_subcort/vc*']);
        for v = 1:length(ciftifiles)
            data = ft_read_cifti_mod([subfolder 'FCPROCESS_SCRUBBED_UWRPMEAN_MSM_V2/cifti_timeseries_normalwall_native_freesurf_222_subcort/' ciftifiles(v).name]);
            vcnum = ciftifiles(v).name(1:7);
            tmask = load([subfolder '/FCPROCESS_SCRUBBED_UWRPMEAN/' vcnum '/total_tmask.txt']);
            data.data = interp_boldts_linear(data.data, logical(tmask));
            ciftidata(v) = {data};
        end
        
        networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
        
        
    end
    
    ROIs = networks; ROIs.data(:) = 0;
    for n = 1:length(networkIDs)
        ROIs.data(networks.data(1:59412)==networkIDs(n)) = networkIDs(n);
    end
    ROIs.data(logical(subnets.data(1:59412))) = subnets.data(logical(subnets.data(1:59412)));
    subnetworkIDs_combined = [subnetworkIDs networkIDs];
    
    roi_inds = false(size(ROIs.data,1),length(subnetworkIDs_combined));
    roi_names = cell(1,length(subnetworkIDs_combined));
    for IDnum = 1:length(subnetworkIDs_combined)
        roi_names{IDnum} = num2str(subnetworkIDs_combined(IDnum));
        roi_inds(:,IDnum) = abs(ROIs.data-subnetworkIDs_combined(IDnum))<.001; roi_inds(59413:end,IDnum)=false;
    end
        
    outname = ['/data/egordon/data1/analysis/Evan/Controls/CON_Subnetworks/lags_altparams/' subname '_Lagmap'];
    %['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/lags/altparams/' subname '_Lagmap'];
    
    for intnum = 1:length(interval_times)
    
        lag_corrmap_maker_func_EGmod(ciftidata,outname, interval_times(intnum) ,roi_inds,3,tr,roi_names)
    
    end
        
    
    
end


