function lag_corrmap_maker_func_EGmod(cifti_times,outname_prefix, int_sec,roi_inds,maxlag_sec,tr,roi_names)
%% Computes lagged correlation relative to grayordinate in cifti
% ciftilist = list of cifti timeseries to concatenate
% int = interpolation factor (0.1)
% index = grayordinate index
% maxlag = max number of timepoints to compute temporal lag (100)
% TL, 03/2024

% [ciftis] = textread(ciftilist,'%s');
% % time_concat = [];
% for r = 1:length(ciftis)
%      cifti_times{r} = ft_read_cifti_mod(ciftis{r});
% %     time_concat = [time_concat cifti_times{r}.data];
% end

int = int_sec ./ tr;
maxlag = maxlag_sec ./ int_sec;

roi_inds = logical(roi_inds);

% Interpolate and upsample
runs_use = 1:length(cifti_times);
t_concat = 1:1:(size(cifti_times{1}.data,2)*length(runs_use));
tq_concat = 1:int:(size(cifti_times{1}.data,2)*length(runs_use));

time_concat_vq = zeros(size(cifti_times{1}.data,1),length(tq_concat));
start = 1;
prevtext = [];
for r = 1:length(runs_use)
    clear vq
    t = 1:1:size(cifti_times{r}.data,2);
    tq = 1:int:size(cifti_times{r}.data,2);
    time_use = cifti_times{runs_use(r)}.data;
    
    vq = zeros(size(time_use,1),length(tq));
    for i = 1:size(time_use,1)
        disptext = ['temporally upsampling run #' num2str(r) ' index #' num2str(i)];
        fprintf([repmat('\b', 1, length(prevtext)) disptext])
        prevtext = disptext;
        vq(i,:) = interp1(t,time_use(i,:),tq,'spline');
        
    end
    stop = start+size(vq,2)-1;
    time_concat_vq(:,start:stop) = vq;
    start = start+size(vq,2);
end
disp(' ')

%figure

%hold

%plot(tq_concat,time_concat_vq(1,:),'k')

% Lag calc

if ~any(size(gcp('nocreate')))
    pool = parpool(10);
end



for roinum = 1:size(roi_inds,2)
    
    if any(roi_inds(:,roinum))
    
    cifti_temp = cifti_times{1};
    r = zeros(size(cifti_temp.data,1),2*maxlag+1);
    prevtext = [];
    roi_mean = mean(time_concat_vq(roi_inds(:,roinum),:),1);
    
    blocks = 50; %divide into blocks to avoid passing parfor a giant matrix
    blocklength = ceil(size(cifti_temp.data,1) ./ blocks);
    for b = 1:blocks
        blockinds = blocklength.*(b-1)+1 : min([blocklength*b size(cifti_temp.data,1)]);
        
        disptext = ['roi ' num2str(roinum) ': running lagged correlations for grayordinates #' num2str(blockinds(1)) ' through #' num2str(blockinds(end))];
        fprintf([repmat('\b', 1, length(prevtext)) disptext])
        prevtext = disptext;
        
        blockdata = time_concat_vq(blockinds,:);
        result = zeros(length(blockinds),2*maxlag+1);
        
        parfor i = 1:length(blockinds)
            result(i,:) = xcorr(roi_mean,blockdata(i,:),maxlag,'coeff');
        end
        
        r(blockinds,:) = result;
    end
    
    % for i = 1:size(cifti_temp.data,1)
    %     disptext = ['grayordinate #' num2str(i)];
    %     fprintf([repmat('\b', 1, length(prevtext)) disptext])
    %     prevtext = disptext;
    %
    %     r(i,:) = xcorr(roi_mean,time_concat_vq(i,:),maxlag,'coeff');
    % end
    disp(' ')
    
    if exist('roi_names','var')
        outname = [outname_prefix '_' roi_names{roinum}];
    else
        outname = [outname_prefix '_roi' num2str(roi_num)];
    end
    
    cifti_temp.data = r;
    cifti_temp.dimord = 'pos_scalar';
    sec_intervals = [-maxlag_sec : int_sec : maxlag_sec];
    for i = 1:length(sec_intervals)
        cifti_temp.mapname(i) = {[num2str(sec_intervals(i)) ' seconds']};
    end
    
    ft_write_cifti_mod([outname '_int' num2str(int_sec) '_spline'],cifti_temp)
    
    [~,maxi] = max(r,[],2);
    maxlag_times = sec_intervals(maxi);
    cifti_temp.data = maxlag_times(:);
    cifti_temp.mapname = {'Lag with maximum correlation'};
    ft_write_cifti_mod([outname '_int' num2str(int_sec) '_spline_maxcorrlag'],cifti_temp)
    
    end
    
end

