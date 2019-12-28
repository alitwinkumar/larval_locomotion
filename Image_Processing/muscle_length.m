function [ROI_Name,ROI,distances,ROI_Image] = muscle_length_v2(image,pixel_size)    
if size(image,4) == 2
        ROI_channel = input('ROI Channel? 1/2')
    else
        ROI_channel = 1;
end
image = imgaussfilt(image(:,:,:,ROI_channel),1.5);
figure('Position',[10,10,1400,1200])    ;
imshow(image(:,:,1,1),[0 150],'InitialMagnification','fit');   
colormap(gca,parula(256));

ROI_Name = input('Which Muscle?');
disp('Draw ROI');
roi_object = drawline('StripeColor','c','LineWidth',3);
input('Verify ROI');
ROI{1} = roi_object.Position;
distances(1) = pdist2([ROI{1}(1,1),ROI{1}(1,2)],[ROI{1}(2,1),ROI{1}(2,2)])*pixel_size
ROI_Image(:,:,:,1) = insertShape(image(:,:,1,1),'Line',reshape(ROI{1}',1,[]));
close all
    
i = 2
while i < size(image,3)+1
    
    figure('Position',[10,10,1400,1200])
    imshow(image(:,:,i,1),[0 150],'InitialMagnification','fit');
    colormap(gca,parula(256)) ;
    disp('Position ROI');
    ROI_temp = drawline('Position',ROI{i-1},'StripeColor','c','LineWidth',3);
    input('Verify ROI');
    ROI{i} = ROI_temp.Position;
    distances(i) = pdist2([ROI{i}(1,1),ROI{i}(1,2)],[ROI{i}(2,1),ROI{i}(2,2)])*pixel_size
    ROI_Image(:,:,:,i) = insertShape(image(:,:,i,1),'Line',reshape(ROI{i}',1,[]));
    i = i + 1
end
distances(1) = pdist2([ROI{1}(1,1),ROI{1}(2,1)],[ROI{2}(1,2),ROI{2}(2,2)])*pixel_size
close all
end
