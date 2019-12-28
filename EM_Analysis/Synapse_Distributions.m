function Synapse_Distributions(nl,clr,direction,alpha,draw_surface,synapse_index)
%This script makes the synapse distribution plots for figures 2 and 3.  
% Synapse_Distributions(nl,clr,direction,alpha,draw_surface)
% Inputs:   nl:  A neuron list 
%           clr: color
%           direction:  Pre (1) or post (2) or both (3)
%           alpha: synapse alpha
%           draw_surface:  Draw an outline of the neuropil for the a/p.
%           This is an option because the neuropil mesh object can be large
%           when exported as an svg.
%           synapse_index:  This is an optional input that will subset the
%           synapses to plot.  It needs to be a cell of logical arrays
%           indicating which synapses should be plotted for each neuron in
%           the neuron list.

% A few details about how this works:  The brain is slightly offset on the
% M/L axis, so to make the left and right hemisegments aligned for the
% cross sectional plots, we rotate all of those points 12 degrees.  The
% synapses are scaled by how many inputs they get.  Basically a pre-synapse
% with one output gets a size of 100, and a pre_synapse with n partners is
% size n*100.  The hisogram plots are 1D kernel density estimates with a 1µm
% bandwith.  


    
load Neuropil_Mesh.mat 

for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 0
        if exist('synapse_index') & direction == 2
            Input{i} = nl(i).Inputs.xyz(synapse_index{i});
        else
            Input{i} = nl(i).Inputs.xyz;
        end
    else
        Input{i} = [nan,nan,nan];
    end
    if isempty(nl(i).Outputs.treenodeID) == 0
        if exist('synapse_index') & direction == 1
            output_index = find(synapse_index{i} == 1);
            outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(output_index(x),:),nl(i).Outputs.polyadics(output_index(x)),1),[1:length(output_index)],'UniformOutput',false);
        else
            outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(x,:),nl(i).Outputs.polyadics(x),1),[1:length(nl(i).Outputs.polyadics)],'UniformOutput',false);
        end    
         Output{i} = cat(1,outputs{:});
    else
        Output{i} = [nan,nan,nan];
    end
end





%% Correct for offset brain
% n_lin = length(nl);
xrot = 0;
yrot = -1;
zrot = -12;


Neuropil_mesh.vert = rotate_pointsV2(Neuropil_mesh.vert,zrot,3)
Neuropil_mesh.vert = rotate_pointsV2(Neuropil_mesh.vert,xrot,1)
Neuropil_mesh.vert = rotate_pointsV2(Neuropil_mesh.vert,yrot,2)

%Find limits each set of axes.  The distances correspond to the A2 and T3
%gaps.  This basically just sets a limit so that we are looking at cross
%sectional points only in A1.  It does not affect the A/P distributions.
axis_lims_MLDV = [min(Neuropil_mesh.vert(Neuropil_mesh.vert(:,3)> 105550 & Neuropil_mesh.vert(:,3)< 144000,1)), max(Neuropil_mesh.vert(Neuropil_mesh.vert(:,3)> 105550 & Neuropil_mesh.vert(:,3)< 144000,1)) , min(Neuropil_mesh.vert(Neuropil_mesh.vert(:,3)> 105550 & Neuropil_mesh.vert(:,3)< 144000,2)), max(Neuropil_mesh.vert(Neuropil_mesh.vert(:,3)> 105550 & Neuropil_mesh.vert(:,3)< 144000,2))];
axis_lims_APDV = [min(Neuropil_mesh.vert(:,3)), max(Neuropil_mesh.vert(:,3)),min(Neuropil_mesh.vert(:,1)), max(Neuropil_mesh.vert(:,1))];



%% Organize synapse coords
if nargin == 2
direction = input('Presynapses:1 Postsynapses:2 Both:3');
else 
end

    if direction == 1
        % Rotate Coords
        if isempty(nl(i).Outputs.treenodeID) == 0
            Output_Coords = cat(1,Output{:});
            Output_Coords = rotate_pointsV2(Output_Coords,zrot,3);
            Output_Coords = rotate_pointsV2(Output_Coords,xrot,1);
            Output_Coords = rotate_pointsV2(Output_Coords,yrot,2);


            d = 'Presynaptic';
            % Get polyadic information
            [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
            a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
            value_counts = [C, a_counts] ;
            polyadics_ap = value_counts(:,4); % Get the size of each polyadic for plotting.

            apdv_points_counts = Output_Coords(:,[1,3]);
            apdv_points_plot = Output_Coords(ia,[1,3]);     
            % Remove points beyond a2/t3 for cross-sectional analysis
            mldv_points_counts = Output_Coords(apdv_points_counts(:,2) >= 105550 & apdv_points_counts(:,2) <= 144000,1:2);
            mldv_points_plot = Output_Coords(ia,1:2);
            mldv_points_plot = mldv_points_plot(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,1:2);
            polyadics_ml = value_counts(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,4);
        else
        end
        
    elseif direction == 2
        if isempty(nl(i).Inputs.treenodeID) == 0
            Input_Coords = cat(1,Input{:});
            Input_Coords = rotate_pointsV2(Input_Coords,zrot,3);
            Input_Coords = rotate_pointsV2(Input_Coords,xrot,1);
            Input_Coords = rotate_pointsV2(Input_Coords,yrot,2);

            apdv_points_counts = Input_Coords(:,[1,3]);
            apdv_points_plot = apdv_points_counts;

            mldv_points_counts = Input_Coords(apdv_points_counts(:,2) > 105550 & apdv_points_counts(:,2) < 144000,1:2);
            mldv_points_plot = mldv_points_counts;

            polyadics_ap = ones(length(apdv_points_plot),1)*3;
            polyadics_ml = ones(length(mldv_points_plot),1)*3;
            d = 'Postsynaptic';
        else
        end
    
    elseif direction == 3
        if isempty(nl(i).Outputs.treenodeID) == 0
            Output_Coords = cat(1,Output{:});
            Output_Coords = rotate_pointsV2(Output_Coords,zrot,3);
            Output_Coords = rotate_pointsV2(Output_Coords,xrot,1);
            Output_Coords = rotate_pointsV2(Output_Coords,yrot,2);

            [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
            a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
            value_counts = [C, a_counts] ;
            polyadics_outputs_ap = value_counts(:,4); % Get the size of each polyadic for plotting.

            apdv_outputs_counts = Output_Coords(:,[1,3]);
            apdv_outputs_plot = Output_Coords(ia,[1,3]);     

            mldv_outputs_counts = Output_Coords(apdv_outputs_counts(:,2) > 105550 & apdv_outputs_counts(:,2) < 144000,1:2);
            mldv_outputs_plot = Output_Coords(ia,1:2);
            mldv_outputs_plot = mldv_outputs_plot(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,1:2);
            polyadics_outputs_ml = value_counts(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,4);
        else
        end
        if isempty(nl(i).Inputs.treenodeID) == 0
            Input_Coords = cat(1,Input{:})
            Input_Coords = rotate_pointsV2(Input_Coords,zrot,3);
            Input_Coords = rotate_pointsV2(Input_Coords,xrot,1);
            Input_Coords = rotate_pointsV2(Input_Coords,yrot,2);

            apdv_inputs_counts = Input_Coords(:,[1,3]);
            apdv_inputs_plot = apdv_inputs_counts;

            mldv_inputs_counts = Input_Coords(apdv_inputs_counts(:,2) > 105550 & apdv_inputs_counts(:,2) < 144000,1:2);
            mldv_inputs_plot = mldv_inputs_counts;

            polyadics_inputs_ap = ones(length(apdv_inputs_plot),1)*3;
            polyadics_inputs_ml = ones(length(mldv_inputs_plot),1)*3;
        else
        end
        
        mldv_points_counts = vertcat(mldv_inputs_counts,mldv_outputs_counts);
        mldv_points_plot = vertcat(mldv_inputs_plot,mldv_outputs_plot);
        polyadics_ml = vertcat(polyadics_inputs_ml,polyadics_outputs_ml);
        
        apdv_points_counts = vertcat(apdv_inputs_counts,apdv_outputs_counts);
        apdv_points_plot = vertcat(apdv_inputs_plot,apdv_outputs_plot);
        
        polyadics_ap = vertcat(polyadics_inputs_ap,polyadics_outputs_ap);
         
        d = 'Combined Synaptic ';
    
    else error('Incorrect Direction')
    end



if isempty(mldv_points_counts) == 0
    if nargin == 4
        subplot(3,7,[12,13,14,19,20,21])
        Neuropil_mesh.vert = Neuropil_mesh.vert(:,[1,3,2]);
        surfaces_v2(Neuropil_mesh,'k',.02,3,[])
        %surfaces(NPM.v(:,[1,3,2]),'k',.02,'-');
    elseif draw_surface == 1
        subplot(3,7,[12,13,14,19,20,21])
        Neuropil_mesh.vert = Neuropil_mesh.vert(:,[1,3,2]);
        surfaces_v2(Neuropil_mesh,'k',.02,3,[])
        %surfaces(NPM.v(:,[1,3,2]),'k',.02,'-');
%         subplot(3,7,[2 3 4 9 10 11])
%         surfaces([NPM.v(NPM.v(:,3)> 105550 & NPM.v(:,3)< 144000,1),NPM.v(NPM.v(:,3)> 105550 & NPM.v(:,3)< 144000,2)],'k',.05,'-')
    else
    end
    
    % Plot the cross sectional points
    subplot(3,7,[2 3 4 9 10 11]); hold on
    scatter(mldv_points_plot(:,1),mldv_points_plot(:,2),polyadics_ml*100,'.','MarkerFaceColor',clr,'MarkerEdgeColor',clr,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
    ax_scat = gca;
    ax_scat.XDir = 'reverse';
    axis off; 
    ylim([axis_lims_MLDV(3),axis_lims_MLDV(4)]);
    xlim([axis_lims_MLDV(1),axis_lims_MLDV(2)]);
    view([0,-90])
    
    % Plot the DV histogram
    subplot(3,7,[1 8]); hold on
    pd_MLDV_2 = fitdist(mldv_points_counts(:,2)*-1,'kernel','Kernel','normal','BandWidth',1000);
    x2_MLDV = axis_lims_MLDV(4)*-1:200:axis_lims_MLDV(3)*-1;
    y2_MLDV = pdf(pd_MLDV_2,x2_MLDV);
    area(x2_MLDV,y2_MLDV,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1);
    xlim([axis_lims_MLDV(4)*-1,axis_lims_MLDV(3)*-1]);
    ax_dv = gca;
    ax_dv.XAxisLocation = 'Top'; xticks([]);
    view([270, 90]); 
    axis off
   
    % Plot the M/L histogram
    subplot(3,7,[16 17 18]); hold on
    pd_MLDV_1 = fitdist(mldv_points_counts(:,1),'kernel','kernel','normal','BandWidth',1000);
    x1_MLDV = axis_lims_MLDV(1):200:axis_lims_MLDV(2);
    y1_MLDV = pdf(pd_MLDV_1,x1_MLDV);
    area(x1_MLDV,y1_MLDV,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1)

    ylh = ylabel({d
            'Density'});
    ylh.Position = [10.0091e+04 3.6002e-05 -1];
    ax_ml = gca; ; ax_ml.XAxisLocation = 'Top';xticks([]); ax_ml.YAxisLocation = 'Left';
    view([0 -90]);
    xlim([axis_lims_MLDV(1),axis_lims_MLDV(2)]);
    set(get(gca,'ylabel'),'rotation',0);
    ax_ml.XDir = 'reverse';
    linkaxes([ax_scat,ax_ml],'x');
    axis off
    
    % Plot the AP Points
    subplot(3,7,[12,13,14,19,20,21]); hold on;
    scatter(apdv_points_plot(:,1),apdv_points_plot(:,2),polyadics_ap*100,'.','MarkerFaceColor',clr,'MarkerEdgeColor',clr,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
    ax_ap1 = gca; 
    ax_ap1.XAxisLocation = 'Top'; xticks([]);
    view([90 90]);
    axis off; axis equal;
    ylim([axis_lims_APDV(1),axis_lims_APDV(2)]);
    xlim([axis_lims_APDV(3),axis_lims_APDV(4)]) ;
    
    
    % Plot the AP histogram
    ax_ap2 = subplot(3,7,[5,6,7]); hold on
    
    pd_apdv = fitdist(apdv_points_counts(:,2),'kernel','kernel','normal','BandWidth',2000);
    x1_apdv = axis_lims_APDV(1):800:axis_lims_APDV(2);
    y1_apdv = pdf(pd_apdv,x1_apdv);
    area(x1_apdv,y1_apdv,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1)
    ylabel(strcat(d,' Density'));
    xlim([axis_lims_APDV(1),axis_lims_APDV(2)]); 
    %xticks([]);
    axis off
    
    %scatter(x1_apdv(y1_apdv==max(y1_apdv)),max(y1_apdv)+.1*max(y1_apdv),20,'MarkerFaceColor',clr','MarkerEdgeColor','k','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
    hold off
    
    
    
    set(gcf,'Color','w');
    set(findall(gcf,'-property','FontSize'),'FontSize',24);
    
else
end
end
