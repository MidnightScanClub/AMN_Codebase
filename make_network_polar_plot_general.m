function make_network_polar_plot_general(networkavgs,line_IDs,network_names)
%% make_network_polar_plot_general(networkavgs,line_IDs,category_names)
%networkavgs is a #networks X 1 vector with average values per network
%line_IDs is a #networks X 1 vector with network IDs associated with each network (power_colors style)
%network_names is a #networks X 1 cell array of network names

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];

line_colors = zeros(length(line_IDs),3);
for linenum = 1:length(line_IDs)
    decimalval = mod(line_IDs(linenum),1);
    if decimalval==0
        thiscolor = power_surf_colormap(line_IDs(linenum),:);
    else
        thiscolor = sum([power_surf_colormap(floor(line_IDs(linenum)),:) .* (1-decimalval) ; power_surf_colormap(ceil(line_IDs(linenum)),:) .* (decimalval)],1);
    end
    line_colors(linenum,:) = thiscolor;
end



networkavgs(networkavgs<0.001) = .001;

thetas = [0:(size(networkavgs,1)-1)]' ./ size(networkavgs,1) .* 2 .* pi;

networkavgs(end+1,:) = networkavgs(1,:);
thetas(end+1) = 2.*pi;

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);
Rmaxlim = ceil(max(networkavgs(:).*10))./10 + .05;
mmpolar_evan(thetas,networkavgs,'lineColors',line_colors,'lineWidth',5,'TTickDelta',360/length(network_names),'TTickLabel',network_names,'TTickScale','degrees','TGridVisible','off','RTickAngle',359.9,'RGridLineStyle','--','FontSize',25,'TZeroDirection','south','RTickLabelValign','bottom','RLimit',[0 Rmaxlim])%,'RTickValue',[.1 : .1 : Rmaxlim],'RTickLabel',{'0.20','0.40','0.60'}