%% Muscle GCAMP Analysis
% This will do the basic image analysis for the muscle GCamP imaging for
% Zarin19. It will load a czi using bioformats.  It will analyze 1 or 2
% color imaging data, and is essentially designed for manually assigning a
% defined ROI onto every single frame. The output is a .mat file that
% contains metadata, image data for QC, and the traces used for later
% analysis.

if input('Load image data?') == 1
    clear all
    Load_Muscle_Data
    % Select start and end frames, crop image
    mov = reshape(image_series_mat_c(:,:,:,size(image_series_mat_c,4)),size(image_series_mat_c,1),size(image_series_mat_c,2),1,size(image_series_mat_c,3));
    implay(mov)
    first = input('Start Frame'); 
    last = input('End Frame');  
    single_cycle = image_series_mat_c(:,:,[first:last],:); 
    clear mov and image_series_mat and image_series_mat_c

    single_cycle = rotate_and_crop(single_cycle,0);
    close all
    Muscle_Data.Pixel_Size = double(meta_data.getPixelsPhysicalSizeX(0).value) % Find pixel size 
    Muscle_Data.time = [0:1:size(single_cycle,3)-1]*double(meta_data.getPixelsTimeIncrement(0).value) % Get framerate
else
    mat_file_outpur_dir = uigetdir
end



%% Make ROIs
mg = input('Measure Gcamp? 1:yes 0:no');
if  mg == 1
    [ROI_Name,ROI,ROI_image] = muscle_ROIs(single_cycle);
    Muscle_Data.ROI_image = ROI_image;
    Muscle_Data.Muscle_Name = ROI_Name;
    Muscle_Data.ROI = ROI;
    clear ROI_Name and ROI and ROI_image
else
end


%% Measure Muscles
if input('Measure Muscle? 1:Yes 0:No') == 1
    [ROI_Name,ROI,distances,ROI_image] = muscle_length(single_cycle,Muscle_Data.Pixel_Size);
    Muscle_Data.Length_ROIs = ROI;
    Muscle_Data.Lenghts = distances;
    Muscle_Data.ROI_image = ROI_image;
    Muscle_Data.Muscle_Name = ROI_Name
    clear ROI_Name and ROI and distances
else
end

%% Measure GCamP
if mg == 1
    for chan = 1:size(single_cycle,4);
        for f_num = 1:size(single_cycle,3);
            clear m and i_temp
            i_temp = single_cycle(:,:,f_num,chan);
            m = createMask(drawpolygon('Position',Muscle_Data.ROI{f_num}),i_temp);
            traces(f_num,chan) = mean(i_temp(m));

        end
        clear f_0 
        f_0 = prctile(traces(:,chan),10)
        df(:,chan) = (traces(:,chan) - f_0)./f_0
    end
    Muscle_Data.df = df
else
end

Muscle_Data.Image_Data = single_cycle

   
%% Save Data
save(strcat(mat_file_outpur_dir,'/',experiment_name,'_',Muscle_Data.Muscle_Name,'_',sprintf('%0.0f',first),'.mat'))