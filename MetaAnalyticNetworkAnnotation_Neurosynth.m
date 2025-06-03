function MetaAnalyticNetworkAnnotation_Neurosynth(mapname,pct_clusters_represented,cifti_intalspace,colormapping)
% MetaAnalyticNetworkAnnotation_Neurosynth(mapname,pct_clusters_represented,cifti_intalspace,[colormapping])
warning off

peak_dilation = 2; %how far to dilate single-coordinate neurosynth peaks
min_cluster_size = 100; %minimum cluster size within map that will be matched to neurosynth peaks
    

if isstruct(mapname)
    cifti = mapname;
else
    cifti = ft_read_cifti_mod(mapname);
end



if ~exist('cifti_intalspace','var')
    cifti_intalspace = false;
end

if ~exist('pct_clusters_represented','var')
    pct_clusters_represented = .5;
end


neurosynthdir = 'Neurosynth/';

load([neurosynthdir 'PreLoadedData.mat'])



[uniquetitles, uniqueinds] = unique(titles);
uniqueids = id(uniqueinds);

neurosynthcoords = single([x y z]);

talspace = strcmp(space,'TAL');
mniconverted = tal2mni(neurosynthcoords(talspace,:));
neurosynthcoords(talspace,:) = mniconverted;

unknownspace = strcmp(space,'UNKNOWN');
neurosynthcoords(unknownspace,:) = [];
titles(unknownspace) = [];

termweights(:,[false (~goodterms')]) = 0;
    


%%

surfL = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii');
surfR = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii');

coordsLR = [surfL.vertices ; surfR.vertices];

neighbors = cifti_neighbors(cifti);

coords_surf = coordsLR(cifti.brainstructure(cifti.brainstructure<3)>0,:);
coords_subcort = cifti.pos(cifti.brainstructure>2,:);
if cifti_intalspace
    coords_subcort = tal2mni(coords_subcort);
end
coords_cifti = single([coords_surf ; coords_subcort]);

ciftiinds_touse = cell(length(uniquetitles),1);
for t = 1:length(uniquetitles)
    thisstudyinds = strcmp(titles,uniquetitles{t});
    thisstudy_coords = neurosynthcoords(thisstudyinds,:);
    [~,thisstudy_cifti_inds] = pdist2(coords_cifti,thisstudy_coords,'euclidean','Smallest',1);
    
    ciftiinds_touse{t} = thisstudy_cifti_inds;
    if exist('peak_dilation','var')
        for d = 1:peak_dilation
            ciftiind_neighbors = neighbors(ciftiinds_touse{t},:); ciftiind_neighbors(isnan(ciftiind_neighbors)) = [];
            ciftiinds_touse{t} = ciftiind_neighbors;
        end
    end
end


mapvalues = unique(cifti.data); mapvalues(mapvalues<1) = [];
anovavalindnums = [];

for v = 1:length(mapvalues)
    value = mapvalues(v);
    thismap = cifti; thismap.data = single(thismap.data==value);

    
    clusters = cifti_cluster(thismap,1,1,min_cluster_size);
    numclusters_min = ceil(size(clusters,2) .* pct_clusters_represented);
    
    
    titles_withmap = cell(0,1);
    ids_withmap = [];
    
    for t = 1:length(uniquetitles)
        
        num_clusters_with_peak = sum(single(any(clusters(ciftiinds_touse{t},:),1)));
        
        if num_clusters_with_peak >= numclusters_min
            titles_withmap(end+1) = uniquetitles(t);
            ids_withmap(end+1) = uniqueids(t);
        end
    end
    titles_withmap = titles_withmap';
    
    %%
    
    study_inds = zeros(length(ids_withmap),1);
    for d = 1:length(ids_withmap)
        study_inds(d) = find(termweights(:,1)==ids_withmap(d));
    end
    
    studies_found_termweights{v} = termweights(study_inds,2:end);
    anovavalindnums = [anovavalindnums ; ones(nnz(study_inds),1) .* v];
end
    
goodterminds = find(goodterms);
anovaPs = zeros(size(goodterminds));
for t = 1:length(goodterminds)
    allvalweights = [];
    for v = 1:length(mapvalues)
        allvalweights = [allvalweights ; studies_found_termweights{v}(:,goodterminds(t))];
    end
    [P,T,STATS,TERMS]=anovan(allvalweights,{anovavalindnums},'display','off');
    anovaPs(t) = P;
    
end
   

p_fdr = FDR(anovaPs,.05);
significant_inds = goodterminds(anovaPs<=.05);
significances = anovaPs(anovaPs<=.05);
fdrsignificant_inds = anovaPs(anovaPs<=.05) <= p_fdr;

top_significant_terms = cell(1,length(mapvalues));
top_significant_term_normweights = cell(1,length(mapvalues));
top_significant_terms_passfdr = cell(1,length(mapvalues));
for t = 1:length(significant_inds)
    meannormweights = zeros(length(mapvalues),1);
    for v = 1:length(mapvalues)
        meannormweights(v) = mean(studies_found_termweights{v}(:,significant_inds(t))) ./ mean(termweights(:,significant_inds(t)+1));
    end
    [maxweight,maxi] = max(meannormweights);
    top_significant_terms{maxi}(end+1) = {[all_termnames{significant_inds(t)} ' ']};
    top_significant_term_normweights{maxi}(end+1) = atanh(1-significances(t));%maxweight;%
    top_significant_terms_passfdr{maxi}(end+1) = fdrsignificant_inds(t);
        
end


figure;
set(gcf,'Position',[318 67 1477 1107]);
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
nonsig_color = [.75 .75 .75];

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.5 .5 .5],50,1)];

if length(mapvalues) == 2
    dims = [1,length(mapvalues)];
elseif length(mapvalues) == 3
    dims = [1,length(mapvalues)];
    set(gcf,'Position',[64 514 2430 660])
elseif length(mapvalues) ==4
    dims = [2,2];
elseif length(mapvalues) <=6
    dims = [2,3];
elseif length(mapvalues) <=8
    dims = [2,4];
elseif length(mapvalues) ==9
    dims = [3,3];
elseif length(mapvalues) <=12
    dims = [3,4];
elseif length(mapvalues) <=15
    dims = [3,5];
else
    dims = [4,4];
end

for v = 1:length(mapvalues)
    h = subplot(dims(1),dims(2),v);
    
    if exist('colormapping','var')
        
        colormapping_index = colormapping(:,1)==mapvalues(v);
        thiscolor = colormapping(colormapping_index,2:end);
        
    else
        decimalval = mod(mapvalues(v),1);
        if decimalval==0
            thiscolor = power_surf_colormap(mapvalues(v),:);
        else
            thiscolor = sum([power_surf_colormap(floor(mapvalues(v)),:) .* (1-decimalval) ; power_surf_colormap(ceil(mapvalues(v)),:) .* (decimalval)],1);
        end
    end
    
    colors = repmat(nonsig_color,length(top_significant_terms{v}),1);
    colors(logical(top_significant_terms_passfdr{v}),:) = repmat(thiscolor,nnz(top_significant_terms_passfdr{v}),1);
    theseweights = top_significant_term_normweights{v};
    theseterms = top_significant_terms{v};
    [minweight,mini] = min(theseweights);
    theseweights(end+1) = minweight ./ 2;%max(theseweights) ./ 5;
    theseterms(end+1) = {' '};
    colors = [colors ; nonsig_color];
    wc = wordcloud(theseterms,theseweights,'Color',colors);
    
1;
end
