function [ROI_Name,ROI,ROI_Image] = muscle_ROIs(image)    

if size(image,4) == 2
        ROI_channel = input('ROI Channel? 1/2')
    else
        ROI_channel = 1;
end
figure('Position',[10,10,1400,1200])    ;
imshow(image(:,:,1,ROI_channel),'InitialMagnification','fit');   
colormap(gca,parula(256));

ROI_Name = input('Which Muscle? Input as "Muscle_xx (n+i)_gcamp_dir" where i is the number of segments away from the first, and dir is either fwd or bwd');
disp('Draw ROI');
roi_object = drawpolygon;
input('Verify ROI');
ROI{1} = roi_object.Position;
ROI_Image(:,:,:,1) = insertShape(image(:,:,1,ROI_channel),'Polygon',reshape(ROI{1}',1,[]));
close all
    
i = 2
while i < size(image,3)+1
    
    figure('Position',[10,10,1400,1200])
    imshow(image(:,:,i,ROI_channel),'InitialMagnification','fit');
    colormap(gca,parula(256)) ;
    disp('Position ROI');
    ROI_temp = drawpolygon('Position',ROI{i-1},'StripeColor','r','LineWidth',2);
    input('Verify ROI');
    ROI{i} = ROI_temp.Position;
    ROI_Image(:,:,:,i) = insertShape(image(:,:,i,ROI_channel),'Polygon',reshape(ROI{i}',1,[]));
    i = i + 1
end
close all
end
