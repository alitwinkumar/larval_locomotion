%% This will generate the panels for figure 5.  

%% Cluster the MNs by synapse similarity
load MNs
savefigs = input('Save Figures? 1:Yes 0:no')
if savefigs == 1
    d = uigetdir
else
end

% First, we want to split the left and right populations of neurons.  We
% want to do this by nerve root, not cell body location, so the RPs are
% going to need to be switched.
% first get neurons from a single hemisegment
clear left_ind and MN_left and MN_right
left_ind = find(contains([MNs(:).Names],'a1l'));
right_ind = find(contains([MNs(:).Names],'a1r'));
MN_right = MNs(right_ind);
MN_left = MNs(left_ind);
% Switch the RPs and MN12
RP_ind = find(contains([MN_right(:).Names],'RP'));
RP2_ind = find(contains([MN_right(:).Names],'RP2'));
MN12_ind = find(contains([MN_right(:).Names],'MN12'));
contra_ind = setdiff(RP_ind,RP2_ind);
RPR = MN_right([contra_ind,MN12_ind]);
MN_right([contra_ind,MN12_ind]) = [];

RP_ind = find(contains([MN_left(:).Names],'RP'));
RP2_ind = find(contains([MN_left(:).Names],'RP2'));
MN12_ind = find(contains([MN_left(:).Names],'MN12'));
contra_ind = setdiff(RP_ind,RP2_ind);
RPL = MN_left([contra_ind,MN12_ind]);
MN_left([contra_ind,MN12_ind]) = [];

MN_right = [MN_right,RPL];
MN_left = [MN_left,RPR];

% Calculate synapse similarity for left and right, then average the
% matrices for clustering.
[sim_mat_l] = synapse_similarity_v2(MN_left,2000,3,[1:3],2);
[sim_mat_r] = synapse_similarity_v2(MN_right,2000,3,[1:3],2);
sim_mat = (sim_mat_l + sim_mat_r)*.5;
group_ind_left = {MN_left(:).SMG};
% group_ind_left = arrayfun(@(x) {num2str(MN_left(x).FMG)},1:length(MN_left))
group_left = unique(group_ind_left);

%Cluster
[h1 t1 perm1 Z] = Synapse_Distance_Clustering_v2(sim_mat,[MN_left(:).Names],0);

colors = {'r', [1,.5,0], 'g' , 'm' , 'c', [.5,0,1]};
reordered_edges = group_ind_left(perm1);
subplot(3,4,[1 5 9]);
ax = get(gca);
lab = ax.YAxis.TickLabels;
for i = 1:length(lab)
    %lab{i}(1:2) = []
    end_name_ind = strfind(lab{i},' ');
    lab{i}(end_name_ind:end) = [];
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-0.01,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',colors{find(contains(group_left,reordered_edges{k}))},'FontSize',18);
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
end
set(gca,'FontSize',24);


%% Get indices for all spatial muscle groups, and plot the distributions of their post-synaptic sites.
DL = ismember({MNs(:).SMG},'DL');
DO_LO = ismember({MNs(:).SMG},'DO/LO');
VL = ismember({MNs(:).SMG},'VL');
VA = ismember({MNs(:).SMG},'VA');
VO = ismember({MNs(:).SMG},'VO');
T = ismember({MNs(:).SMG},'T');
smgs = unique({MNs(:).SMG});

isn = boolean(DL+DO_LO+VL+VO);
sn = boolean(VA+T);

figure('Position',[1100 0 2000 900],'rend','Painters'); hold on
Synapse_Distributions(MNs(isn),'b',2,.25,1);
Synapse_Distributions(MNs(sn),'m',2,.25,0);


figure('Position',[1100 0 2000 900],'rend','Painters'); hold on
%figure;
for i = 1:length(smgs)-1
    Synapse_Distributions(MNs(ismember({MNs(:).SMG},smgs{i})),colors{i},2,.5,0)
end
    Synapse_Distributions(MNs(ismember({MNs(:).SMG},smgs{end})),colors{end},2,.5,1)
    
%% Synapse Distributions of CMuGs

map = {'r',[1,.65,.25],'c','b'}
figure('rend','painters','pos',[1 1 1940 952]); hold on
for i = 1:3
Synapse_Distributions(MNs([MNs(:).FMG] == i),map{i},2,.25,0)

end
Synapse_Distributions(MNs([MNs(:).FMG] == 4),map{4},2,.25,1)
subplot(3,7,[2 3 4 9 10 11])
legend({'CMuG F1','CMuG F2','CMuG F3','CMuG F4'},'Location','SouthWest') ; legendmarkeradjust(50)
set(gca,'FontSize',16)
if savefigs == 1
    saveas(gcf,strcat(d,'/','FMG','_Distributions'),'svg')
    close all
else
end

figure('rend','painters','pos',[1 1 1940 952]); hold on
for i = 1:3
Synapse_Distributions(MNs([MNs(:).BMG] == i),map{i},2,.25,0)

end
Synapse_Distributions(MNs([MNs(:).BMG] == 4),map{4},2,.25,1)

subplot(3,7,[2 3 4 9 10 11])
legend({'CMuG B1','CMuG B2','CMuG B3','CMuG B4'},'Location','SouthWest') ; legendmarkeradjust(50)
set(gca,'FontSize',16)
if savefigs == 1
    saveas(gcf,strcat(d,'/','BMG','_Distributions'),'svg')
    close all
else
end

%% Synapse distributions for differentially active muscles compared to their co-active groups
VO_sans28 = VO
VO_sans28(37:38) = 0
figure('rend','painters','pos',[1 1 1940 952]); hold on
Synapse_Distributions(MNs([MNs(:).FMG]==3 - VO_sans28),[.5,.5,.5],2,.5,1)
Synapse_Distributions(MNs([MNs(:).BMG]==1 - VO_sans28),[.25,.25,.25],2,.5,1)
Synapse_Distributions(MNs(VO_sans28),'r',2,.5,0)


mn18_ind = zeros(1,length(MNs))
mn18_ind(17:18) = 1
figure('rend','painters','pos',[1 1 1940 952]); hold on
Synapse_Distributions(MNs([MNs(:).FMG] == 2 - mn18_ind),[.5,.5,.5],2,.5,1)
%Synapse_Distributions(MNs([MNs(:).BMG] == 4 - mn18_ind),[.25,.25,.25],2,.5,1)
Synapse_Distributions(MNs(17:18),'r',2,.5,0)


mn11_ind = zeros(1,length(MNs))
mn11_ind(5:6) = 1
figure('rend','painters','pos',[1 1 1940 952]); hold on
Synapse_Distributions(MNs([MNs(:).FMG] == 1 - mn11_ind),[.5,.5,.5],2,.5,0)
Synapse_Distributions(MNs([MNs(:).BMG] == 4 - mn11_ind),[.25,.25,.25],2,.5,0)

Synapse_Distributions(MNs(5:6),'r',2,.5,1)

mn2_ind = zeros(1,length(MNs))
mn2_ind(21:22) = 1
figure('rend','painters','pos',[1 1 1940 952]); hold on
Synapse_Distributions(MNs([MNs(:).FMG] == 1 - mn2_ind),[.5,.5,.5],2,.5,0)
Synapse_Distributions(MNs([MNs(:).BMG] == 3 - mn2_ind),[.25,.25,.25],2,.5,0)
Synapse_Distributions(MNs(21:22),'r',2,.5,1)
%%

mn_list = MNs([MNs(:).BMG] == 4 - mn18_ind)
mn_list = mn_list(1:2:end)

F1_sansmn_cat = concatinate_synapses(mn_list,'CMG_sansMN');
MN_cat = concatinate_synapses(MNs(17),'MN2');

dmat_cmgvsmn = pdist2(MN_cat.Inputs.xyz,F1_sansmn_cat.Inputs.xyz);
dmin_mnvscmg = min(dmat_cmgvsmn,[],2)
dmin_cmgvsmn = min(dmat_cmgvsmn,[],1)

figure;
histogram(dmin_mnvscmg(:),0:1000:5e4,'Normalization','probability')
hold on
histogram(dmin_cmgvsmn(:),0:1000:5e4,'Normalization','probability')
legend({'MN vs CMuG','CMuG vs MN'})

distances = discretize(dmin_mnvscmg,[0:100:5e4]);

map = plasma(max(distances));
figure; hold on
for i = 1:max(distances)
    scatter3(MN_cat.Inputs.xyz(distances==i,1),MN_cat.Inputs.xyz(distances==i,2),MN_cat.Inputs.xyz(distances==i,3),1000,map(i,:),'.')
    scatter3(F1_sansmn_cat.Inputs.xyz(:,1),F1_sansmn_cat.Inputs.xyz(:,2),F1_sansmn_cat.Inputs.xyz(:,3),200,'k','.','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
    
end
%% Statistics for distributions

% Concatinate synapses for each SMG

% put all the smg indices together
for i = 1:length(MN_left)
    smg_indices(i) = find(ismember(smgs,{MN_left(i).SMG}))
end

for i = 1:max(smg_indices)
    SMG_cat(i) = concatinate_synapses(MN_left(smg_indices(:,i)),smgs{i})
    SMG_cat(i).Inputs.xyz = rotate_pointsV2(SMG_cat(i).Inputs.xyz,0,3);
    SMG_cat(i).Inputs.xyz = rotate_pointsV2(SMG_cat(i).Inputs.xyz,-1,1);
    SMG_cat(i).Inputs.xyz = rotate_pointsV2(SMG_cat(i).Inputs.xyz,-12,2);
end


for i = 1:length(SMG_cat)
    for j = 1:length(SMG_cat)
        [h,p] = kstest2(SMG_cat(i).Inputs.xyz(:,1),SMG_cat(j).Inputs.xyz(:,1));
        p_smg_ml(i,j) = p
        
        [h,p] = kstest2(SMG_cat(i).Inputs.xyz(:,2),SMG_cat(j).Inputs.xyz(:,2));
        p_smg_dv(i,j) = p
        
        [h,p] = kstest2(SMG_cat(i).Inputs.xyz(:,3),SMG_cat(j).Inputs.xyz(:,3));
        p_smg_ap(i,j) = p
   
    end
end

% Concatinate synapses for FMGs and BMGs

for i = 1:4
    FMG_cat(i) = concatinate_synapses(MN_left([MN_left(:).FMG] == i),'FMG')
    FMG_cat(i).Inputs.xyz = rotate_pointsV2(FMG_cat(i).Inputs.xyz,0,3);
    FMG_cat(i).Inputs.xyz = rotate_pointsV2(FMG_cat(i).Inputs.xyz,-1,1);
    FMG_cat(i).Inputs.xyz = rotate_pointsV2(FMG_cat(i).Inputs.xyz,-12,2);
    
    BMG_cat(i) = concatinate_synapses(MN_left([MN_left(:).BMG] == i),'BMG')
    BMG_cat(i).Inputs.xyz = rotate_pointsV2(BMG_cat(i).Inputs.xyz,0,3);
    BMG_cat(i).Inputs.xyz = rotate_pointsV2(BMG_cat(i).Inputs.xyz,-1,1);
    BMG_cat(i).Inputs.xyz = rotate_pointsV2(BMG_cat(i).Inputs.xyz,-12,2);
end
% Calculate pairwise comparisons
for i = 1:4
    for j = 1:4
        % FMGs
        [h,p] = kstest2(FMG_cat(i).Inputs.xyz(:,1),FMG_cat(j).Inputs.xyz(:,1));
        p_fmg_ml(i,j) = p;

        [h,p] = kstest2(FMG_cat(i).Inputs.xyz(:,2),FMG_cat(j).Inputs.xyz(:,2));
        p_fmg_dv(i,j) = p;

        [h,p] = kstest2(FMG_cat(i).Inputs.xyz(:,3),FMG_cat(j).Inputs.xyz(:,3));
        p_fmg_ap(i,j) = p;

        %BMGs
        [h,p] = kstest2(BMG_cat(i).Inputs.xyz(:,1),BMG_cat(j).Inputs.xyz(:,1));
        p_bmg_ml(i,j) = p;

        [h,p] = kstest2(BMG_cat(i).Inputs.xyz(:,2),BMG_cat(j).Inputs.xyz(:,2));
        p_bmg_dv(i,j) = p;

        [h,p] = kstest2(BMG_cat(i).Inputs.xyz(:,3),BMG_cat(j).Inputs.xyz(:,3));
        p_bmg_ap(i,j) = p;
    end
end

map = cbrewer('seq','Blues',64);
figure; colormap(map)
subplot(3,1,1)
imagesc(triu(p_smg_ap))
yticks([1:6])
yticklabels(smgs)
xticks([1:6])
xticklabels(smgs)
subplot(3,1,2)
imagesc(triu(p_smg_dv))
yticks([1:6])
yticklabels(smgs)
xticks([1:6])
xticklabels(smgs)
subplot(3,1,3)
imagesc(triu(p_smg_ml))
yticks([1:6])
yticklabels(smgs)
xticks([1:6])
xticklabels(smgs)

