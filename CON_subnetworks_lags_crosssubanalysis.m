cd /data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/temp/

subnames = {'SIC01','SIC02','SIC03','MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10','ME01','ME02','ME03','ME04'};

subnetworkIDs = [10.8, 16.5, 9.5, 11.4];
networkIDs = [1.5 5 8 10 11 17];
networks_toinclude_ordered = [8 16.5 11.4 10.8 9.5 1.5 10];

flipsign = true;

allsubs_TDmat = zeros(length(subnetworkIDs),length(networks_toinclude_ordered),length(subnames)) .* NaN;

%abs_lag_thresh = 1;

FC_thresh = .1;

for subnum = 1:length(subnames) %
    
    subname = subnames{subnum};
    subnets = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/' subname '_con_subnetworks_autodetected.dtseries.nii']);
    
    
    if strcmp(subname(1:3),'SIC')
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
    elseif strcmp(subname(1:2),'ME')
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        networks = ft_read_cifti_mod([infomapdir '/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
    elseif strcmp(subname(1:3),'MSC')
        networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' subname '_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_recolored_wCMI.dscalar.nii']);
    end
    
    thissub_subnetworks = networks; thissub_subnetworks.data(:) = 0;
    for n = 1:length(networkIDs)
        thissub_subnetworks.data(networks.data(1:59412)==networkIDs(n)) = networkIDs(n);
    end
    thissub_subnetworks.data(logical(subnets.data(1:59412))) = subnets.data(logical(subnets.data(1:59412)));
    
    
    
    
    thissub_subnetworks.data(59413:end) = 0;
    
    
    
    
    
    for n = 1:length(subnetworkIDs)
        
        thissub_lagsmapfile = ['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/lags/' subname '_Lagmap_' num2str(subnetworkIDs(n)) '_int0.03_spline_maxcorrlag.dscalar.nii'];
        if exist(thissub_lagsmapfile,'file')
        
        lagsmap = ft_read_cifti_mod(thissub_lagsmapfile);
        if flipsign
            lagsmap.data = -lagsmap.data;
        end
        
        crosscorrmap = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/lags/' subname '_Lagmap_' num2str(subnetworkIDs(n)) '_int0.03_spline.dscalar.nii']);
        negativeinds = crosscorrmap.data(:,ceil(size(crosscorrmap.data,2)./2)) < FC_thresh;
        lagsmap.data(negativeinds) = NaN;
        
        for n3 = 1:length(networks_toinclude_ordered)
%             if networks_toinclude(n3)==10
%                 networkinds = thissub_subnetworks.data(1:59412)==10 | thissub_subnetworks.data(1:59412)==11 | thissub_subnetworks.data(1:59412)==17;
%             else
                networkinds = abs(thissub_subnetworks.data(1:59412)-networks_toinclude_ordered(n3))<.001;
%            end
            if any(networkinds)
                values = lagsmap.data(networkinds);
                allsubs_TDmat(n,n3,subnum) = nanmean(values);
            else
               allsubs_TDmat(n,n3,subnum) = NaN;
            end
            
            
        end
        end
    end
    
    
end


allsubs_meanlags = squeeze(nanmean(allsubs_TDmat,1));

make_network_bargraph(networks_toinclude_ordered,nanmean(allsubs_meanlags,2)+.02,nanstd(allsubs_meanlags,0,2)./sqrt(length(subnames)),false)
set(gca,'XAxisLocation','origin')

groupvar = repmat([1:size(allsubs_meanlags,1)],1,size(allsubs_meanlags,2));
groupvar = groupvar(:);
[P,T,STATS,TERMS]=anovan(allsubs_meanlags(:),groupvar(:));

for i = 1:length(networks_toinclude_ordered)
    for j = i+1:length(networks_toinclude_ordered)
        [~,P,~,STATS] = ttest(allsubs_meanlags(j,:),allsubs_meanlags(i,:));
        disp([num2str(i) ' vs ' num2str(j) ': t(14) = ' sprintf('%02.2f',STATS.tstat) '; p = ' sprintf('%02.3f',P)])
    end
end

% allsubs_reorglags = zeros(length(networks_toinclude), length(subnets) .* length(subnames));
% for n3 = 1:length(networks_toinclude)
%     for n = 1:length(subnets)
%         for subnum = 1:length(subnames)
%             ind = (n-1) .* length(subnets) + subnum;
%             allsubs_reorglags(n3,ind) = allsubs_TDmat(n,n3,subnum);
%         end
%     end
% end
