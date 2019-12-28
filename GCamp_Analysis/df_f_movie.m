% Make Df/f movie
fname = uigetfile
load(fname)
image = Muscle_Data.ROI_image;
image = permute(image,[1,2,4,3]);
image = im2double(image);
f0_image = prctile(image(:,:,:,2),1,3);

df_image = (image(:,:,:,2)-f0_image);
lowthresh = prctile(df_image(:),2);
upperthresh = prctile(df_image(:),98)*1.5;

df_mov = mat2gray(df_image,[lowthresh,upperthresh]);

outputVideo = VideoWriter(fullfile('df_f_movie_fwd.avi'));
outputVideo.FrameRate = 20;
open(outputVideo)
for i = 1:length(Muscle_Data.time)
    s{i} = strcat(sprintf('%0.2f',Muscle_Data.time(i)),'s');
end
for i = 1:size(df_image,3)
    figure;
    imagesc(df_mov(:,:,i),[0 1])
    c = colorbar
    c.Label.String = 'Normalized DF/F'
    colormap plasma
    set(gcf,'Color','w')
    axis off
    text(0,5,s{i},'Color','w', 'fontsize', 14)
    set(gca,'FontSize',18)
    img = getframe(gcf);
    close all
    writeVideo(outputVideo,img);
end
close(outputVideo);
%%
disp('Movie Size')
size(df_mov)
interval = input('Choose montage interval')
df_mov_short = df_mov(:,:,interval);
montage(df_mov_short,'Size',[size(df_mov_short,3),1]) ; colormap plasma;
c = colorbar
c.Label.String = 'Normalized DF/F'
set(gcf,'Color','w')
set(gca,'FontSize',28)
s{interval}
%%
for i = 1:size(df_mov,3)
    figure;
    img = imshow(df_mov(:,:,i))
    df_ms(:,:,:,i) = insertText(img.CData, [0 0], s{i}, 'AnchorPoint','LeftTop', ...
        'BoxColor', 'white', 'fontsize', 14);
end
