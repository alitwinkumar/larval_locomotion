crawls = single_crawl
muscle_numbers = [1:1:30];

for i = 1:length(crawls)
    clear crawl_length and crawl_end and crawl_start and crawl_interval and new_crawl_time and df and crawl_trace and pre_crawl and post_crawl
    crawl_length = crawls(i).last_onset - crawls(i).onset_point;
    crawl_end = find(crawls(i).muscle_data(1).normalized_time*100 == crawls(i).last_onset);
    crawl_start = find(crawls(i).muscle_data(1).normalized_time*100 == crawls(i).onset_point);
    crawl_interval = (crawl_end-crawl_start);
    new_crawl_time = 0:100/(crawl_interval*2):100;
    df = [crawls(i).muscle_data(:).normalized_df];
    df_s = [crawls(i).muscle_data(:).normalized_df_ratio_smooth];
    
    
    crawl_trace = nan(length(new_crawl_time),size(df,2));
    crawl_trace(ceil(crawl_interval*.5):ceil(crawl_interval*1.5),:) = df(crawl_start:crawl_end,:);
    
    crawl_trace_s = nan(length(new_crawl_time),size(df_s,2));
    crawl_trace_s(ceil(crawl_interval*.5):ceil(crawl_interval*1.5),:) = df_s(crawl_start:crawl_end,:);
    
    
    pre_crawl = crawl_trace(1:ceil(.5*crawl_interval)+1,:);
    post_crawl = crawl_trace(ceil(1.5*crawl_interval)+1:end,:);
    
    
    pre_crawl_s = crawl_trace_s(1:ceil(.5*crawl_interval)+1,:);
    post_crawl_s = crawl_trace_s(ceil(1.5*crawl_interval)+1:end,:);
    
    df_pre_crawl = df(1:crawl_start-1,:)
    df_post_crawl = df(crawl_end+1:end,:)
    
    df_pre_crawl_s = df_s(1:crawl_start-1,:)
    df_post_crawl_s = df_s(crawl_end+1:end,:)
    
    if size(pre_crawl,1)-size(df_pre_crawl,1) > 0
        crawl_trace(1:size(pre_crawl,1)-size(df_pre_crawl,1),:) = nan
        crawl_trace(size(pre_crawl,1)-size(df_pre_crawl,1)+1:size(pre_crawl,1),:) = df_pre_crawl;
        
        crawl_trace_s(1:size(pre_crawl_s,1)-size(df_pre_crawl_s,1),:) = nan
        crawl_trace_s(size(pre_crawl_s,1)-size(df_pre_crawl_s,1)+1:size(pre_crawl_s,1),:) = df_pre_crawl_s;
    
    elseif size(pre_crawl,1)-size(df_pre_crawl,1) == 0
        crawl_trace(1:size(pre_crawl,1) ,:) = df_pre_crawl;
        crawl_trace_s(1:size(pre_crawl_s,1),:) = df_pre_crawl_s;

    else
        crawl_trace(1:size(pre_crawl,1) ,:) = df_pre_crawl(end - size(pre_crawl,1)+1 : end,:);
        crawl_trace_s(1:size(pre_crawl_s,1) ,:) = df_pre_crawl_s(end - size(pre_crawl_s,1)+1 : end,:);
    end
    
     
    if size(post_crawl,1) - size(df_post_crawl,1) > 0;
        crawl_trace(ceil(crawl_interval*1.5)+(size(post_crawl,1) - size(df_post_crawl,1)):end,:) = nan;
        crawl_trace(ceil(crawl_interval*1.5)+1:ceil(crawl_interval*1.5)+size(df_post_crawl,1),:) = df_post_crawl;
        
        crawl_trace_s(ceil(crawl_interval*1.5)+(size(post_crawl_s,1) - size(df_post_crawl_s,1)):end,:) = nan;
        crawl_trace_s(ceil(crawl_interval*1.5)+1:ceil(crawl_interval*1.5)+size(df_post_crawl_s,1),:) = df_post_crawl_s;
    else
        crawl_trace(ceil(crawl_interval*1.5)+1:end,:) = df_post_crawl(1:end+(size(post_crawl,1) - size(df_post_crawl,1)),:)
    
        crawl_trace_s(ceil(crawl_interval*1.5)+1:end,:) = df_post_crawl_s(1:end+(size(post_crawl_s,1) - size(df_post_crawl_s,1)),:)

    end
    
   for ii = 1:length(crawls(i).muscle_data);
       crawls(i).muscle_data_aligned(ii).crawl_trace = crawl_trace(:,ii);
       crawls(i).muscle_data_aligned(ii).crawl_trace_s = crawl_trace_s(:,ii);
       crawls(i).muscle_data_aligned(ii).norm_time = new_crawl_time;
       crawls(i).muscle_data_aligned(ii).muscle_number = crawls(i).muscle_data(ii).Muscle_Number;
       
       crawls(i).muscle_data_aligned(ii).alignment_data.crawl_length = crawl_length;
       crawls(i).muscle_data_aligned(ii).alignment_data.crawl_start = crawl_start;
       crawls(i).muscle_data_aligned(ii).alignment_data.crawl_end = crawl_end;
       
   end
end

map = hsv(length(crawls));
% 
 for i = 1:length(crawls)
     
figure; hold on
    for ii = 1:length(crawls(i).muscle_data_aligned)
        
     plot(crawls(i).muscle_data_aligned(ii).norm_time,crawls(i).muscle_data_aligned(ii).crawl_trace_s,'Color',clrs{cmug_lookup(crawls(i).muscle_data_aligned(ii).muscle_number)})
    end
 end