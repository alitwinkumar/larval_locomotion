
%% Load Muscle Image Data
close all; clear all
if input('Stimulus Data? 1:Yes 0:No  :') == 1
    [image_series_mat_c,meta_data,f_name,stim_meta] = load_image_data(1);
else
    [image_series_mat_c,meta_data,f_name] = load_image_data(0);
end
experiment_name = strsplit(f_name,'.');
experiment_name = experiment_name{1,1};
% if input('Save Figs? Yes:1  :') == 1
%     fprintf('Choose Output Directory for Figures');
%     psfile = strcat(uigetdir,'/',experiment_name,'.ps')
%     if exist(psfile,'file')==2;delete(psfile);end
% end
fprintf('Select Output Directory for Data');
mat_file_outpur_dir = uigetdir
clear f_name