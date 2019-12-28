% This script will generate the length vs gcamp plots.  First it loads all
% the gcamp data, then all of the length data, and then plots the gcamp vs
% length data for muscles in which both exist.  

%% Load files
input_directory = uigetdir
disp('Load Gcamp data')
gcamp_files = dir([input_directory,'/*.mat']);
% Compile muscle data into structures
muscles = [1:1:30];
cmugs = [1,1,1,1,3,1,1,5,1,1,2,1,1,1,4,4,4,5,3,3,5,5,5,5,5,3,3,4,3,1];      % These were defined by Aref
%Dictionary for muscle to cmugs
cmug_lookup = containers.Map(muscles,cmugs);     

 %include muscle files that have ratiometric measurements or all data 2:ratiometric 1:all
ratiometric_only = 2                                                       

for i = 1:length(gcamp_files)
    m_f = load(strcat(input_directory,'/',gcamp_files(i).name));
        if size(m_f.traces,2) == ratiometric_only
            if min(m_f.traces(:,1)) > 0
                muscle_files(i).Muscle_Name = m_f.Muscle_Data.Muscle_Name(1:15);

                m_num_ind = regexp(muscle_files(i).Muscle_Name,'[0-9]');
                if m_num_ind(1) == 1
                    muscle_files(i).Muscle_Name = strcat('Muscle_',muscle_files(i).Muscle_Name);
                    muscle_files(i).Muscle_Name = muscle_files(i).Muscle_Name(1:15);
                else
                end

                muscle_files(i).Muscle_Number = str2num(muscle_files(i).Muscle_Name(8:9));
                muscle_files(i).Muscle_Data = m_f.Muscle_Data;
                muscle_files(i).experiment_name = m_f.experiment_name;
                name_breaks = strfind(gcamp_files(i).name,'_');
                
                direction_ind = strfind(muscle_files(i).Muscle_Data.Muscle_Name,'wd')
                muscle_files(i).crawl_direction = muscle_files(i).Muscle_Data.Muscle_Name(direction_ind-1:direction_ind+1)
                
              
                muscle_files(i).normalized_time = [0:1/(length(m_f.Muscle_Data.df)-1):1]; % normalized time is just the time series set to 0 to 1
                window = find(muscle_files(i).normalized_time<.05);
                window = window(end);
                if mod(window,2) == 0
                    window = window+1;
                end


                muscle_files(i).raw_traces_ratio = sgolayfilt(m_f.traces(:,2)./m_f.traces(:,1),1,window); % Filter the raw traces and make ratiometric
                
                muscle_files(i).raw_traces = sgolayfilt(m_f.traces(:,size(m_f.traces,2)),1,window); % Filter raw traces non-ratiometric
                muscle_files(i).raw_traces_ratio_smooth = sgolayfilt(m_f.traces(:,2)./m_f.traces(:,1),1,window*3);
                
                muscle_files(i).df_ratio = (muscle_files(i).raw_traces_ratio - prctile(muscle_files(i).raw_traces_ratio(1:end/3),10))./ prctile(muscle_files(i).raw_traces_ratio(1:end/3),10,1); % Ratiometric Df/f.  F0 is calculated as the 10th percentile of the raw signal.
                muscle_files(i).df_ratio_smooth = (muscle_files(i).raw_traces_ratio_smooth - prctile(muscle_files(i).raw_traces_ratio_smooth(1:end/3),10))./ prctile(muscle_files(i).raw_traces_ratio_smooth(1:end/3),10,1);
                muscle_files(i).df = (muscle_files(i).raw_traces - prctile(muscle_files(i).raw_traces(1:end/3),10))./ prctile(muscle_files(i).raw_traces(1:end/3),10,1); %regular df/f




                muscle_files(i).normalized_df = normalize(muscle_files(i).df,1,'range') ; %normalized df, range normalization
                muscle_files(i).normalized_df_ratio = normalize(muscle_files(i).df_ratio,1,'range'); %normalized ratiometric DF, range normalization    
                muscle_files(i).normalized_df_ratio_smooth = normalize(muscle_files(i).df_ratio_smooth,1,'range');
                muscle_files(i).cmug = cmug_lookup(muscle_files(i).Muscle_Number);
                muscle_files(i).ID = i;

                % TO USE RATIOMETRIC FOR ANALYSIS: RUN THIS LINE.  OTHERWISE
                % COMMENT IT OUT.
                %
                muscle_files(i).df = muscle_files(i).df_ratio;
                muscle_files(i).normalized_df = muscle_files(i).normalized_df_ratio;
                %
                %
                %

                muscle_files(i).dvdt = diff([0;muscle_files(i).df])./diff([0;muscle_files(i).normalized_time(:)]); % First derivative of the DF trace

                startIndex = regexp(muscle_files(i).Muscle_Name,'(');
                %endIndex = regexp(muscle_files(i).Muscle_Name,')');
                if contains(muscle_files(i).Muscle_Name(startIndex:end),'+') == 1;
                    muscle_files(i).segment_index = str2double(muscle_files(i).Muscle_Name(startIndex+3));
                else
                    muscle_files(i).segment_index = 0;
                end
                muscle_files(i).frame_start = str2double(gcamp_files(i).name(name_breaks(end)+1:end-4))
%                 if contains(gcamp_files(i).name(end-7:end),'_') == 1;
%                     muscle_files(i).frame_start = str2double(gcamp_files(i).name(end-6:end-4));
%                 else
%                     muscle_files(i).frame_start = str2double(gcamp_files(i).name(end-7:end-4));
%                 end
               clear startIndex and endIndex
            else
            end
        else
        clear m_f
        end
end
muscle_files = muscle_files(~cellfun(@isempty,{muscle_files.Muscle_Name})); %remove muscle files without ratiometric data.

% get indices for combinations of variables that will determine individual
% experiments.  Experiment names are unique identifiers for each movie,
% wave numbers are identifiers of each segment within a movie, direction is
% for crawl direction and frame index is for the start time of a given
% crawl within a movie.
[~,experiment_index] = ismember({muscle_files(:).experiment_name},unique({muscle_files(:).experiment_name}));
%[muscle_files(:).experiment_number] = deal(experiment_index);

[~,frame_index] = ismember([muscle_files(:).frame_start],unique([muscle_files(:).frame_start]));
%[muscle_files(:).frame_index] = deal(frame_index);

wave_index = [muscle_files(:).segment_index];
direction_index = contains({muscle_files(:).crawl_direction},'fwd');
experiment_combinations = transpose([direction_index;experiment_index;frame_index;wave_index]); %define a matrix of indices for each muscle file
single_cycles = unique(experiment_combinations,'rows'); %find unique combinations of indices.  The length of this matrix is the number of n's 

for i = 1:length(muscle_files)
    muscle_files(i).exp_comb = experiment_combinations(i,:)
end

muscle_files_gcamp = muscle_files;
clearvars -except muscle_files_gcamp and experiment_combinations
%%
disp('Load muscle length data')
input_directory = uigetdir
gcamp_files = dir([input_directory,'/*.mat']);

for i = 1:length(gcamp_files)
m_f = load(strcat(input_directory,'/',gcamp_files(i).name));
    muscle_files(i).Muscle_Name = m_f.Muscle_Data.Muscle_Name(1:15);

    m_num_ind = regexp(muscle_files(i).Muscle_Name,'[0-9]');
    if m_num_ind(1) == 1
        muscle_files(i).Muscle_Name = strcat('Muscle_',muscle_files(i).Muscle_Name);
        muscle_files(i).Muscle_Name = muscle_files(i).Muscle_Name(1:15);
    else
    end

    muscle_files(i).Muscle_Number = str2num(muscle_files(i).Muscle_Name(8:9));
    muscle_files(i).Muscle_Data = m_f.Muscle_Data;
    muscle_files(i).experiment_name = m_f.experiment_name;
    name_breaks = strfind(gcamp_files(i).name,'_');
    muscle_files(i).crawl_direction = gcamp_files(i).name(name_breaks(end)-9:name_breaks(end)-7);
    muscle_files(i).length = m_f.Muscle_Data.Lenghts    
    muscle_files(i).normalized_time = [0:1/(length(muscle_files(i).length)-1):1]; % normalized time is just the time series set to 0 to 1
    
     startIndex = regexp(muscle_files(i).Muscle_Name,'n+');
    if contains(muscle_files(i).Muscle_Name(startIndex:startIndex+2),'_') == 1;
        muscle_files(i).segment_index = 0;
    else
        muscle_files(i).segment_index = str2double(muscle_files(i).Muscle_Name(startIndex+2));
    end
    if contains(gcamp_files(i).name(end-7:end),'_') == 1;
        muscle_files(i).frame_start = str2double(gcamp_files(i).name(end-6:end-4));
    else
        muscle_files(i).frame_start = str2double(gcamp_files(i).name(end-7:end-4));
    end

end
[~,experiment_index] = ismember({muscle_files(:).experiment_name},unique({muscle_files(:).experiment_name}));
%[muscle_files(:).experiment_number] = deal(experiment_index);

[~,frame_index] = ismember([muscle_files(:).frame_start],unique([muscle_files(:).frame_start]));
%[muscle_files(:).frame_index] = deal(frame_index);

wave_index = [muscle_files(:).segment_index];
direction_index = contains({muscle_files(:).crawl_direction},'fwd');



experiment_combinations_l = transpose([direction_index;experiment_index;frame_index;wave_index]); %define a matrix of indices for each muscle file
for i = 1:length(muscle_files)
    
    muscle_files(i).exp_comb = experiment_combinations_l(i,:)
end
muscle_files_length = muscle_files; clear muscle_files
clearvars -except muscle_files_length and muscle_files_gcamp and experiment_combinations
%%
for i = 1:length(muscle_files_length)
    clear experiment_index and frame_index and direction_index and segment_index and muscle_index and mfg
    experiment_index = contains({muscle_files_gcamp(:).experiment_name},muscle_files_length(i).experiment_name)
    frame_index = [muscle_files_gcamp(:).frame_start] == muscle_files_length(i).frame_start
    direction_index = contains({muscle_files_gcamp(:).crawl_direction},muscle_files_length(i).crawl_direction)
    segment_index = [muscle_files_gcamp(:).segment_index] == muscle_files_length(i).segment_index
    muscle_index = [muscle_files_gcamp(:).Muscle_Number] == muscle_files_length(i).Muscle_Number
    mfg = muscle_files_gcamp(experiment_index & frame_index & direction_index & segment_index & muscle_index);
    muscle_files_gcamp_comparison(i) = mfg(1)
end

clearvars -except muscle_files_gcamp_comparison and muscle_files_length
%% 
crawl_ind = contains({muscle_files_length(:).crawl_direction},'fwd')
fwd_length = muscle_files_length(crawl_ind == 1)
fwd_gcamp = muscle_files_gcamp_comparison(crawl_ind == 1)

for i = 1:length(fwd_gcamp)
    clear fwd_df and fwd_l and l_g and l_l
    fwd_df = fwd_gcamp(i).normalized_df_ratio
    fwd_l = transpose(fwd_length(i).length)
    
    l_g = length(fwd_df)
    l_l = length(fwd_l)
    
    if l_g > l_l
        fwd_df = fwd_df(1:l_l)
    else
        fwd_l = fwd_l(1:l_g)
    end

    fwd_df = fwd_df(2:end)
    fwd_l = fwd_l(2:end)
       
    [peak_contract,peak_ind] = min(fwd_l)
    fs = 1/(length(fwd_l)-1)
    f_rates(i) = fs;
    time = [0-(fs*peak_ind):fs:fs*(length(fwd_l)-peak_ind-1)]
    ts{i} = timeseries([fwd_df,fwd_l],time)
    figure; hold on
    yyaxis left
    plot(time,fwd_l)
    yyaxis right
    plot(time,fwd_df)
end

fs_new = max(fs);
fs_new_ind = find(max(fs));
figure; hold on
for i = 1:length(ts)
    ts_new{i} = resample(ts{i},ts{1}.Time);
    df_new(:,i) = ts_new{i}.Data(:,1);
    length_new(:,i) = ts_new{i}.Data(:,2);
    disp(length(ts_new{i}.time))
    yyaxis left
    plot(ts_new{i}.time,ts_new{i}.Data(:,1));
    yyaxis right
    plot(ts_new{i}.time,ts_new{i}.Data(:,2));
end
    
df_mean = nanmean(df_new,2);
df_std = nanstd(df_new,[],2);

length_mean = nanmean(length_new,2);
length_std = nanstd(length_new,[],2);
figure; hold on
set(gca,'FontSize',18)
set(gcf,'Color','w')
xlabel('Normalized Time')

yyaxis right
shadedErrorBar(ts_new{1}.Time,length_mean,length_std,'lineprops','r')
ylabel('Muscle Length (µm)')
yyaxis left
ylabel('DF/F')
shadedErrorBar(ts_new{1}.Time,df_mean,df_std,'lineprops','b')



%%
bwd_length = muscle_files_length(crawl_ind == 0)
bwd_gcamp = muscle_files_gcamp_comparison(crawl_ind == 0)
bwd_df = [bwd_gcamp(:).normalized_df_ratio];
clear time_b and df_new and length_new and ts_b and ts_b_new


for i = 1:length(bwd_gcamp)
    clear bwd_df and bwd_l and l_g and l_l
    bwd_df = bwd_gcamp(i).normalized_df_ratio
    bwd_l = transpose(bwd_length(i).length)
    
    l_g = length(bwd_df)
    l_l = length(bwd_l)
    
    if l_g > l_l
        bwd_df = bwd_df(1:l_l)
    else
        bwd_l = bwd_l(1:l_g)
    end

    bwd_df = bwd_df(2:end)
    bwd_l = bwd_l(2:end)
       
    [peak_contract,peak_ind] = min(bwd_l)
    fs_b = 1/(length(bwd_l)-1)
    f_rates(i) = fs_b;
    time_b = [0-(fs_b*peak_ind):fs_b:fs_b*(length(bwd_l)-peak_ind-1)]
    ts_b{i} = timeseries([bwd_df,bwd_l],time_b)
    figure; hold on
    yyaxis left
    plot(time_b,bwd_l)
    yyaxis right
    plot(time_b,bwd_df)
end

fs_new = max(fs_b);
fs_new_ind = find(max(fs_b));
figure; hold on
for i = 1:length(ts_b)
    ts_b_new{i} = resample(ts_b{i},ts_b{1}.Time);
    df_new(:,i) = ts_b_new{i}.Data(:,1);
    length_new(:,i) = ts_b_new{i}.Data(:,2);
    disp(length(ts_b_new{i}.time))
    yyaxis left
    plot(ts_b_new{i}.time,ts_b_new{i}.Data(:,1));
    yyaxis right
    plot(ts_b_new{i}.time,ts_b_new{i}.Data(:,2));
end
    
df_mean = nanmean(df_new,2);
df_std = nanstd(df_new,[],2);

length_mean = nanmean(length_new,2);
length_std = nanstd(length_new,[],2);
figure; hold on
set(gca,'FontSize',18)
set(gcf,'Color','w')
xlabel('Normalized Time')

yyaxis right
shadedErrorBar(ts_b_new{1}.Time,length_mean,length_std,'lineprops','r')
ylabel('Muscle Length (µm)')
yyaxis left
ylabel('DF/F')
shadedErrorBar(ts_b_new{1}.Time,df_mean,df_std,'lineprops','b')
