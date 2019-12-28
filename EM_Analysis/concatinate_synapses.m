function [meta_neuron] = concatinate_synapses(neuron_list,meta_name)
    for i = 1:length(neuron_list)
        names{i} = neuron_list(i).Names;
        skids(i) = neuron_list(i).SkIDs;
        if isempty(find(contains(fieldnames(neuron_list),'Group_1'))) == 0
            group1{i} = neuron_list(i).Group_1;
        else
        end
        if isempty(find(contains(fieldnames(neuron_list),'Group_2'))) == 0
            group2{i} = neuron_list(i).Group_2;
        else
        end
        
        if length(neuron_list(i).Outputs.treenodeID) >= 1
            presyn_treenodeID{i} = neuron_list(i).Outputs.treenodeID;
            presyn_coords{i} = neuron_list(i).Outputs.xyz;
            presyn_polyadics{i} = neuron_list(i).Outputs.polyadics;
            presyn_conid{i} = neuron_list(i).Outputs.conid;
            presyn_skid{i} = repmat(neuron_list(i).SkIDs,length(presyn_treenodeID{i}),1)
            if isfield(neuron_list(i).Outputs,'Partner_skids')
                presyn_partner_skids{i} = neuron_list(i).Outputs.Partner_skids;
            else
            end
                
        else
            presyn_treenodeID{i} = [];
            presyn_coords{i} = [nan,nan,nan];
            presyn_polyadics{i} = [];
            presyn_conid{i} = [];
            presyn_skid{i} = repmat(neuron_list(i).SkIDs,length(presyn_treenodeID{i}),1)

        end
        
        if length(neuron_list(i).Inputs.treenodeID) > 1
            postsyn_treenodeID{i} = neuron_list(i).Inputs.treenodeID;
            postsyn_coords{i} = neuron_list(i).Inputs.xyz;
            postsyn_conid{i} = neuron_list(i).Inputs.conid;
            postsyn_skid{i} = repmat(neuron_list(i).SkIDs,length(postsyn_treenodeID{i}),1)
        else
            postsyn_treenodeID{i} = [];
            postsyn_coords{i} = [nan,nan,nan];
            postsyn_conid{i} = [];
        end
    end
    meta_neuron.Name = meta_name;
    meta_neuron.Included_Neurons = names;
    meta_neuron.SkIDs = skids;
   
    
    meta_neuron.Inputs.xyz = cat(1,postsyn_coords{:});
    meta_neuron.Inputs.xyz(isnan(meta_neuron.Inputs.xyz(:,1)),:) = [];
    meta_neuron.Inputs.conid = cat(1,postsyn_conid{:});
    meta_neuron.Inputs.treenodeID = cat(2,postsyn_treenodeID{:});
    meta_neuron.Inputs.associated_skid = cat(1,postsyn_skid{:});
    
    meta_neuron.Outputs.xyz = cat(1,presyn_coords{:});
    meta_neuron.Outputs.conid = cat(1,presyn_conid{:});
    meta_neuron.Outputs.treenodeID = cat(1,presyn_treenodeID{:});
    meta_neuron.Outputs.polyadics = cat(1,presyn_polyadics{:});
    meta_neuron.Outputs.associated_skid = cat(1,presyn_skid{:});
    
    if isfield(neuron_list(i).Outputs,'Partner_skids') & length(neuron_list(i).Outputs.treenodeID) > 1
        meta_neuron.Outputs.Partner_skids = cat(1,presyn_partner_skids{:});
    else
    end
end
