function [sim_mat,pre_post] = synapse_similarity_v2(nl,sigma,n_dims,dims,pre_post)
%Computes the synapse similarity score using the method from Schlegel 2016 
%The output is a mean similarity between the synapses of neurons of list 1.  Neuron lists come from load_synapses
% [sim_mat,pre_post] = synapse_similarity(nl,sigma,n_dims,dims,pre_post)
if nargin == 5
elseif nargin == 4
    pre_post = input('Presynapses (1) / Postsynapses (2) / Both(3)')
else
    error('Not enough input args')
end

for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 0
        Input_Coords{i} = nl(i).Inputs.xyz;
    else
    end
    if isempty(nl(i).Outputs.treenodeID) == 0
        outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(x,:),nl(i).Outputs.polyadics(x),1),[1:length(nl(i).Outputs.polyadics)],'UniformOutput',false);
        Output_Coords{i} = cat(1,outputs{:});
    else
    end
end


if n_dims < 3 
    remove_axis = setdiff([1:3],dims)
    for i = 1:length(nl)
        if exist('Input_Coords')
            Input_Coords{i}(:,remove_axis) = [];
        end
        if exist('Output_Coords')
            Output_Coords{i}(:,remove_axis) = [];
        end
    end
end

for i = 1:length(nl)
    clear D_isis;
    if pre_post == 2
        D_isis = pdist2(Input_Coords{i},Input_Coords{i},'euclidean');
    elseif pre_post == 1
        D_isis = pdist2(Output_Coords{i},Output_Coords{i},'euclidean');
    else error('Incorrect Direction')
    end

    if isempty(D_isis)
        D_isis = inf;
    else
    end
    for j = 1:length(nl)
        clear D_jkjk and D_isjk and k_vals and k_idx 
        if pre_post == 2
            D_jkjk = pdist2(Input_Coords{j},Input_Coords{j},'euclidean'); % intraneuron synapse distances for neuron j
        elseif pre_post == 1
            D_jkjk = pdist2(Output_Coords{j},Output_Coords{j},'euclidean'); % intraneuron synapse distances for neuron j
        else error('Incorrect Direction')
        end
        
        if isempty(D_jkjk) == 1
            D_jkjk = inf;     
        else 
        end
        
        
        if pre_post == 2
            D_isjk = pdist2(Input_Coords{i},Input_Coords{j},'euclidean'); % interneuron synapse distances for neurons i,j
        elseif pre_post == 1
            D_isjk = pdist2(Output_Coords{i},Output_Coords{j},'euclidean'); % interneuron synapse distances for neurons i,j
        else error('Incorrect Direction')
        end
        
        
        if isempty(D_isjk) == 1
            D_isjk = inf;
        else 
        end
        
        
        [k_vals,k_idx] = min(D_isjk,[],2);  % for every synapse of neuron i, smallest distance to a synapse of neuron j and corresponding index.
        clear sim and n_is and n_jk and Dsk
        for s = 1:size(D_isjk,1);
            Dsk(s) = k_vals(s);
            n_is(s) = numel(find(D_isis(s,:)<sigma))/size(D_isis,2);
            n_jk(s) = numel(find(D_jkjk(k_idx(s),:)<sigma))/size(D_jkjk,2);
            sim(s) = exp((-1*(Dsk(s)^2))/(2*(sigma^2)))*exp(-1*abs(n_is(s)-n_jk(s))/(n_is(s)+n_jk(s)));
        end
        similarity_matrix(i,j) = mean(sim);
        
    end
end
similarity_matrix(find(isnan(similarity_matrix))) = 0;
sim_mat = .5*(similarity_matrix+similarity_matrix');

end

    