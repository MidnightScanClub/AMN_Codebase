subnames = {'MSC01','MSC03','MSC04','MSC05','MSC07','MSC08','MSC09','MSC10','SIC01','SIC02'};

%TASKS
%1 : 0-Tongue+1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg_zstat.dscalar.nii
%5 : 0_Tongue_zstat.dscalar.nii
%16 : 2_R_Hand_zstat.dscalar.nii
%13 : 1_L_Hand_zstat.dscalar.nii
%19 : 3_L_Leg_zstat.dscalar.nii
%20 : 4_R_Leg_zstat.dscalar.nii
%22 : GlassTask_tp3+4_zstat.dscalar.nii
%23 : WordTask_tp3+4_zstat.dscalar.nii

tasks = [5 13 16 19 20 22 23];
tasknames = {'Tongue', 'LHand','RHand','LLeg','RLeg','Spatial','Verbal'};

subnet_IDs = [10.8, 16.5, 9.5, 11.4];
subnetworknames = {'Feedback', 'Decision','Action','ParsMarg'};
%%

allsubs_anterior = []';
allsubs_posterior = []';
allsubs_lateral = []';
allsubs_paracentral = []';
activation = []';
allsubs_max_task_activation = []';

taskvals = zeros(length(tasks),length(subnames),length(subnet_IDs));

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    activation = [];
    if strcmp(subname,'SIC01')
            task = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/Tasks/MSC02_allcontrastsV2_smooth2.55.dscalar.nii']);
    end
    if strcmp(subname,'SIC02')
            task = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/Tasks/MSC06_allcontrastsV2_smooth2.55.dscalar.nii']);
    end
            
    if strcmp(subname(1:3),'MSC')  
            task = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/Tasks/' subname '_allcontrastsV2_smooth2.55.dscalar.nii']);
    end
  
    
    CONsubnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/' subname '_con_subnetworks_autodetected.dtseries.nii']);

    for eachtask = 1:length(tasks)
        
        for subnetnum = 1:length(subnet_IDs)
            taskvals(eachtask,subnum,subnetnum) = mean(task.data(abs(CONsubnetworks.data(1:59412)-subnet_IDs(subnetnum))<.01,tasks(eachtask)));
        end
    end
    
    
    
end

save('/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/CON_subnetworks/autodetected/task_activation_autodetected.mat','taskvals');

%%
for n = 1:length(tasknames)
    thistask = squeeze(taskvals(n,:,:));
    thissubnets = repmat(subnet_IDs,size(thistask,1),1);
    thesesubs = repmat([1:size(thistask,1)],1,size(thistask,2));
    [P,T,STATS,TERMS]=anovan(thistask(:),{thissubnets(:) thesesubs(:)},'display','off','random',2);
    disp(['Task ' tasknames{n} ': F(' num2str(T{2,3}) ',' num2str(T{4,3}) ') = ' num2str(T{2,6}) ', p = ' num2str(P(1)) ', corrected p = ' num2str(P(1) .* length(tasknames))])
end

%%
ps = [];
for n = 1:length(tasknames)
    string = [tasknames{n} ', '];
    for s = 1:length(subnet_IDs)
        for s2 = (s+1):length(subnet_IDs)
            [H,P,CI,STATS] = ttest(taskvals(n,:,s),taskvals(n,:,s2));
            %disp([tasknames{n} '; subnetwork ' subnetworknames{s} ' vs ' subnetworknames{s2} ': T(' num2str(STATS.df) ')=' num2str(STATS.tstat) '; p=' num2str(P)])
            string = [string 'T(' num2str(STATS.df) ')=' sprintf('%02.2f',STATS.tstat) ';p=' sprintf('%02.3f',P) ', '];
            ps(end+1) = P;
        end
    end
    disp(string)
end

%%
make_network_polar_plot_general(squeeze(nanmean(taskvals,2)),subnet_IDs,tasknames)

