% This script makes all of the components for figures 2 and 3 of 
% Zarin, 2019 with the exception of the muscle length vs muscle gcamp panels.  
% The data used to generate the figures are data structures
% of individual muscles that are output from the mucle_gcamp.mat script.  

% %% Load files

direction = input('Which direction? 1:fwd 0:bwd')
[single_crawl,cmug_lookup,cmugs] = Load_Muscle_Files(2,direction)
muscles = [1:1:30];  clc
%% Calculate Max DF Values

clear muscle_traces
muscle_traces = struct('muscle_number',{},'muscle_data_unaligned',{});
for i = 1:30
   
    muscle_traces(i).muscle_number = num2str(i)
    muscle_traces(i).cmug = cmug_lookup(i);
   
    for ii = 1:length(single_crawl)
        if isempty(find([single_crawl(ii).muscle_data(:).Muscle_Number] == i)) == 0
                   
            muscle_traces(i).muscle_data_unaligned(ii) = single_crawl(ii).muscle_data(find([single_crawl(ii).muscle_data(:).Muscle_Number] == i));
        else
        
        end
    end
    if isempty(muscle_traces(i).muscle_data_unaligned) == 0
        muscle_traces(i).muscle_data_unaligned = muscle_traces(i).muscle_data_unaligned(~cellfun(@isempty,{muscle_traces(i).muscle_data_unaligned.Muscle_Name}));
   else
   end

end
sizes = cellfun(@length,{muscle_traces(:).muscle_data_unaligned})
muscle_traces(sizes < 1) = []



%% Plot max DF/F values
figure('rend','painters','pos',[1 1 450 952]); hold on
for i = 1:length(muscle_traces)
    clear df and max_df
    
    for ii = 1:length(muscle_traces(i).muscle_data_unaligned)
        df = muscle_traces(i).muscle_data_unaligned(ii).df_ratio_smooth;
        max_df(ii) = max(df);
    end
        errorbar(str2double(muscle_traces(i).muscle_number),mean(max_df),std(max_df),'--','Color','k','LineWidth',2)
        scatter(repmat(str2double(muscle_traces(i).muscle_number),length(max_df),1),max_df,1000,'r','.')

    end
plot([0 30],[1 1])
xticks(1:1:30)
yticks(0:2:30)
ylim([0 30])
set(gcf,'Color','w')
set(gca,'FontSize',18)
ylabel('Max DF/F')
xlabel('Muscle')

clear muscle_traces and sizes and max_df and df


%% Find and remove inactive muscles  
% Find any muscles that do not show at least a 100% DF/F signal and remove
% them from the analysis.
for i = 1:length(single_crawl)
    for ii = 1:length(single_crawl(i).muscle_data)
        
        df_max(ii) = max(single_crawl(i).muscle_data(ii).df_ratio_smooth)
        if df_max(ii) < 1 
            inactive_index(ii) = 0
        else
            inactive_index(ii) = 1
        end
        
       
    end  
    single_crawl(i).muscle_data(inactive_index == 0) = []   
    clear inactive_index and df_max
end
single_crawl(cellfun(@isempty,{single_crawl(:).muscle_data})) = []

%% Find Onset times for each crawl cycle   
% First threshold for completeness.  Only include crawls that have >40% of
% the muscles included.

single_crawl([single_crawl(:).completeness] < .4) = [] %Threshold for all crawl cycles that have 40% of the muscles or greater
clear df_values and df_values_s and map and tp and rp and r and tlin and df_mat and muscle_nums
   
figure('rend','painters','pos',[1,1,1000,1000]); hold on;

% First run PCA on the calcium traces, then take the first two PCs and
% "unwrap" them.  When plotted in 2D space, the crawl cycle appears as a
% trajectory away from, and back to some origin.  If we plot this in polar
% coordinates, there is a peak that corrisponds to the trace being furthest
% from the origin.  This is peak we want to detect and call the midpoint of
% the crawl cycle.  
for i = 1:length(single_crawl);
    
    df_values_s = [single_crawl(i).muscle_data(:).normalized_df_ratio_smooth];
   

    [coeff,score,latent,tsquared,explained,mu] = pca(df_values_s); 
    map = jet(size(score,1));
    ii = 1;
    while ii < size(score,1)-1;
        plot([score(ii,1),score(ii+1,1)],[score(ii,2),score(ii+1,2)],'Color',map(ii,:),'LineWidth',3);
        ii = ii+1;
    end
    [~, rp] = cart2pol(score(:,1)', score(:,2)'); % Convert (unwrap) first two PCs to polar coord.
    r{i} = rp; clear p
    
    
    scatter(score(1,1),score(1,2),100,'k','.');
end
axis off;
set(gcf,'Color','w');
c = colorbar;
c.Label.String = 'Normalized Time';
colormap jet;
set(gca,'FontSize',18)

MINPROM = .5;
MINDISTANCE = 5;
MINHEIGHT = 2;
% 
for i = 1:length(r)
    figure; hold on; clear norm_time and p1 and p2 and pks and locs and w and p
    norm_time = transpose(single_crawl(i).muscle_data(1).normalized_time)*100;
    df_values = [single_crawl(i).muscle_data(:).normalized_df_ratio_smooth];
    plot(norm_time,df_values)
    findpeaks(r{i},1/norm_time(2),'minpeakdistance', MINDISTANCE,'MinPeakProminence',MINPROM,'Annotate','Extents','WidthReference','HalfHeight')
    [pks,locs,w,p,wxpk] = findpeaks_new(r{i},1/norm_time(2), 'minpeakdistance', MINDISTANCE,'MinPeakProminence',MINPROM,'Annotate','Extents','WidthReference','HalfHeight');
    pk{i} = pks;
    if length(pks) > 1
        wxpk = wxpk(pks == max(pks),:);
    else
    end
    p1 = wxpk(1)
    p2 = wxpk(2)
   [~,crawl_end] = min(abs(norm_time-p2));
   [~,crawl_start] =  min(abs(norm_time-p1));
       scatter(norm_time(crawl_start),1,100,'r')
       scatter(norm_time(crawl_end),1,100,'b')
    
     single_crawl(i).onset_point = norm_time(crawl_start);
     single_crawl(i).last_onset = norm_time(crawl_end);
     xlim([0 100])
end
close all
%clear p1 and p2 and pks and locs and w and p and wxpk and crawl_start and crawl_end and MINPROM and MINDISTANCE and MINHEIGHT and c and r and df_values_s 
%% Plot crawl speed variability

for i = 1:length(single_crawl)
    crawl_duration = single_crawl(i).last_onset - single_crawl(i).onset_point
    for ii = 1:length(single_crawl(i).muscle_data)
        if isfield(single_crawl(i).muscle_data(ii).Muscle_Data,'time') == 1
            crawl_time(i) = (crawl_duration/100)*single_crawl(i).muscle_data(ii).Muscle_Data.time(end)
        else
            ii = ii-1
        end
    end
end


figure('pos',[1,1,200,1000]); hold on
scatter(ones(length(crawl_time),1),crawl_time,1000,'r','.')
errorbar(mean(crawl_time),std(crawl_time),'k','LineWidth',2)
ylim([0 10])
set(gcf,'Color','w')
ylabel('Time (s)')
set(gca,'FontSize',18)
clear crawl_time and crawl_duration
%% Pick a Direction and align to onset and offset points

crawl_alignment
clear crawl_interval and crawl_length and crawl_trace and crawl_trace_s and df_post_crawl and df_post_crawl_s and df_pre_crawl and df_pre_crawl_s and post_crawl and post_crawl_s and pre_crawl and pre_crawl_s and new_crawl_time
  close all
%% Compile Single Muscle Traces

clear muscle_traces
muscle_traces = struct('muscle_number',{},'muscle_data_aligned',{},'muscle_data_unaligned',{});
for i = 1:30
   
    muscle_traces(i).muscle_number = num2str(i)
    muscle_traces(i).cmug = cmug_lookup(i);
   
    for ii = 1:length(crawls)
        if isempty(find([crawls(ii).muscle_data_aligned(:).muscle_number] == i)) == 0
        muscle_traces(i).muscle_data_aligned(ii) = crawls(ii).muscle_data_aligned(find([crawls(ii).muscle_data(:).Muscle_Number] == i));
        
        muscle_traces(i).muscle_data_unaligned(ii) = crawls(ii).muscle_data(find([crawls(ii).muscle_data(:).Muscle_Number] == i));
        else
        
        end
    end
    
   if isempty(muscle_traces(i).muscle_data_aligned) == 0
        muscle_traces(i).muscle_data_aligned = muscle_traces(i).muscle_data_aligned(~cellfun(@isempty,{muscle_traces(i).muscle_data_unaligned.Muscle_Name}));
   else
   end
   if isempty(muscle_traces(i).muscle_data_aligned) == 0
        muscle_traces(i).muscle_data_unaligned = muscle_traces(i).muscle_data_unaligned(~cellfun(@isempty,{muscle_traces(i).muscle_data_unaligned.Muscle_Name}));
   else
   end
   
   if length(muscle_traces(i).muscle_data_unaligned) < 1
       muscle_traces(i).muscle_data_aligned = []
       muscle_traces(i).muscle_data_unaligned = []
   else
   end
end

% Remove muslce 7 since it it's trace is unreliable
muscle_traces(7).muscle_data_aligned = []
muscle_traces(7).muscle_data_unaligned = []

%  Redo cmugs
 if direction == 1     
    cmugs = [3,1,2,2,2,1,1,3,2,1,1,2,2,1,3,3,3,2,2,3,4,4,4,4,2,2,2,3,2,1]  % clustering defined by fwd 
 else 
    cmugs =  [2,3,2,2,3,2,2,3,2,1,4,2,2,1,1,1,1,4,3,3,4,4,4,4,4,3,3,2,3,1] % clustering defined by bwd
 end
cmug_lookup = containers.Map(muscles,cmugs);  

sizes = cellfun(@length,{muscle_traces(:).muscle_data_aligned})

cmugs(sizes <3) = []
muscle_traces(sizes < 3) = []
names = {muscle_traces(:).muscle_number}
clear crawls
%% Resample everything
clrs = {'r','y','g','b','m'};
clear all_data
figure; hold on
for i = 1:length(muscle_traces)
    
    if isempty(muscle_traces(i).muscle_data_aligned) == 0
        for j = 1:length(muscle_traces(i).muscle_data_aligned);
            trace_length(i,j) = length(muscle_traces(i).muscle_data_aligned(j).norm_time);
        end
    end
end
sample_size = ceil(mean(trace_length(trace_length>0)));
for i = 1:length(muscle_traces)
   %figure; hold on
    if isempty(muscle_traces(i).muscle_data_aligned) == 0;
        for ii = 1:length(muscle_traces(i).muscle_data_aligned);
            
            clear ts and ts_new and df_s and df and time
            df_s = (muscle_traces(i).muscle_data_aligned(ii).crawl_trace_s);
            df = (muscle_traces(i).muscle_data_aligned(ii).crawl_trace);
            time = muscle_traces(i).muscle_data_aligned(ii).norm_time;
            ts{ii} = timeseries(df_s,time);
            ts_new{i}{ii} = resample(ts{ii},[0:100/(sample_size-1):100]);
            muscle_traces(i).muscle_data_aligned(ii).crawl_trace_n = ts_new{i}{ii}.Data;
            muscle_traces(i).muscle_data_aligned(ii).norm_time_n = ts_new{i}{ii}.Time;
           %plot(muscle_traces(i).muscle_data_aligned(ii).norm_time_n,muscle_traces(i).muscle_data_aligned(ii).crawl_trace_n, clrs{5})
            %plot(muscle_traces(i).muscle_data_unaligned(ii).normalized_time,muscle_traces(i).muscle_data_unaligned(ii).normalized_df_ratio_smooth, clrs{4})

            %plot(time,df,'Color',[.8,.8,.8])
        end
        
        all_data(:,i) = nanmean([muscle_traces(i).muscle_data_aligned(:).crawl_trace_n],2);
        shadedErrorBar(muscle_traces(i).muscle_data_aligned(1).norm_time_n,nanmean([muscle_traces(i).muscle_data_aligned(:).crawl_trace_n],2),nanstd([muscle_traces(i).muscle_data_aligned(:).crawl_trace_n],[],2),'lineprops',clrs{cmug_lookup(str2double(muscle_traces(i).muscle_number))})
        title(muscle_traces(i).muscle_data_unaligned(1).Muscle_Name)
    else
     end
end

clear ts and ts_new and df_s and df and time
set(gcf,'Color','w')
ylabel('Normalized DF/F')
xlabel('Normalized Time')
set(gca,'FontSize',18)
plot([50,50],[0 1],'r','LineWidth',3,'LineStyle','--')





figure; hold on
for i = 1:length(muscle_traces)
    for ii = 1:length(muscle_traces(i).muscle_data_aligned)
        df = muscle_traces(i).muscle_data_aligned(ii).crawl_trace_n;
        time = muscle_traces(i).muscle_data_aligned(ii).norm_time_n;
        plot(time,df,'Color',[.8 .8 .8])
    end
end
shadedErrorBar(muscle_traces(1).muscle_data_aligned(1).norm_time_n,nanmean(all_data,2),nanstd(all_data,[],2))
set(gcf,'Color','w')
ylabel('Normalized DF/F')
xlabel('Normalized Time')
set(gca,'FontSize',18)
 plot([25,25],[0 1],'r','LineWidth',3,'LineStyle','--')
 plot([75,75],[0 1],'r','LineWidth',3,'LineStyle','--') 
figure;
for i = 1:max(cmugs)
    subplot(4,1,i); hold on
    clear mt and df and time
    mt = muscle_traces(find(cmugs == i))
    for ii = 1:length(mt)
        df = [mt(ii).muscle_data_aligned(:).crawl_trace_n];
        time = mt(ii).muscle_data_aligned(1).norm_time_n;
        shadedErrorBar(time,nanmean(df,2),nanstd(df,[],2),'lineprops',clrs{i})
    end
    ylabel('Normalized DF/F')
    xlabel('Normalized Time')
    set(gca,'FontSize',18)
    plot([25,25],[0 1],'r','LineWidth',3,'LineStyle','--')
    plot([75,75],[0 1],'r','LineWidth',3,'LineStyle','--') 

    
end
clear mt and df and time
set(gcf,'Color','w')

figure; hold on
for i = 1:max(cmugs)
    
    clear df_vals and time
    df_vals = all_data(:,find(cmugs == i))
    time = muscle_traces(5).muscle_data_aligned(1).norm_time_n;
    plot(time,df_vals,'Color',clrs{i},'LineStyle','-.')
    shadedErrorBar(time,nanmean(df_vals,2),nanstd(df_vals,[],2),'lineprops',clrs{i})
end
clear df_vals and time
set(gcf,'Color','w')
ylabel('Normalized DF/F')
xlabel('Normalized Time')
set(gca,'FontSize',18)
plot([50,50],[0 1],'r','LineWidth',3,'LineStyle','--')


 %% Cluster average muscle traces
all_data(isnan(all_data)) = 0;
ad = all_data'
[ h t perm dist coph eva_gap eva_sil] = muscle_clustering(ad(:,1:round(end*.5)),names)

reordered_edges = perm;

subplot(3,4,[1 5 9]);
ax = get(gca);  
lab = ax.YAxis.TickLabels;
set(gca,'FontSize',24);
xticks([0:1:5]);


loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-0.01,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',clrs{cmug_lookup(str2double(lab(ind)))},'FontSize',24);
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
end
disp(coph)
 
%clear h and t and perm and dist and coph and eva_gap and eva_sil and lab and ax and ad
if direction == 1
    save('fwd_traces.mat')
elseif direction == 0
    save('bwd_traces.mat')
else
    error('wrong_direction')
end
%% Muscle clustering algorithm
function [ h t perm dist coph eva_gap eva_sil ] = muscle_clustering(x,Names)
    [s_c] = PCA_Analysis(x)
    dist = pdist(s_c,'euclidean');
    d_mat = squareform(dist);
    
tic, Z = linkage(dist,'ward'); toc

op = optimalleaforder(Z,dist,'Transformation','linear')

eva_gap = evalclusters(x,'linkage','gap','KList',[1:6])
eva_sil = evalclusters(x,'linkage','silhouette','KList',[1:6])

coph = cophenet(Z,pdist(x,'euclidean'))
[h t perm] = dendrogram(Z,4,'Orientation','Left','Labels',Names,'ColorThreshold',1,'reorder',fliplr(op))
figure; 
silhouette(x,t)
figure('rend','painters','pos',[1 1 550 952]); hold on
subplot(3,4,[1 5 9 ]);

[h t perm] = dendrogram(Z,0,'Orientation','Left','Labels',Names,'ColorThreshold',1,'reorder',fliplr(op));

subplot(3,4,[2 3 4 6 7 8 10 11 12 ]);
if nargin == 3
imagesc(x(perm,perm)); axis xy;
else
imagesc(x(perm,:)); axis xy
end
axis off
c = colorbar ; c.Label.String = 'Normalized DF/F';

colormap plasma
set(gcf,'Color','w')
set(gca,'FontSize',18)

end

%%
function [Selected_Comps] = PCA_Analysis(data)
dataCorr = corrcoef(data');
figure
imagesc(dataCorr,[-1 1]); xlabel(''); ylabel(''); axis xy; colorbar; title('correlation matrix'); set(gcf,'Color','w')
[coeff,score,latent,tsquared,explained,mu] = pca(data); %looks different if not transposed
colormap(parula)
figure; 
for i = 1:4
subplot(2,2,i)
scatter(score(:,i),score(:,i+1))
title(strcat('PC ',sprintf('%.0f',i)))
set(gca,'FontSize',18)
end
set(gcf,'Color','w')

hold on
figure
plot(latent); xlabel ' component #'; ylabel 'eigenvalue'; %most variance is explained by the number of components where curve drops off 
low_lat = find(latent<.05);
display(explained);
display(low_lat)
Selected_Comps = score(:,1:min(low_lat))
set(gca,'FontSize',18)
set(gcf,'Color','w')
end
