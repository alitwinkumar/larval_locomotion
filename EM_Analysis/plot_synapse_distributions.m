function plot_synapse_distributions(nl_cell,clr,direction,alpha)
        figure('pos',[1975 153 1375 715],'rend','painters'); hold on
    
    if iscell(nl_cell)
        n = length(nl_cell)

        if n > 1
            if iscell(clr)
                Synapse_Distributions(nl_cell{n},clr{n},direction,alpha,1)
                for i = n-1:-1:1
                    Synapse_Distributions(nl_cell{i},clr{i},direction,alpha,0)
                end
            else
                Synapse_Distributions(nl_cell{n},clr(n,:),direction,alpha,1)
                for i = n-1:-1:1
                    Synapse_Distributions(nl_cell{i},clr(i,:),direction,alpha,0)
                end
            end
        else
            Synapse_Distributions(nl_cell{1},clr,direction,alpha,1)
        end
    else
        Synapse_Distributions(nl_cell,clr,direction,alpha,1)
    end
    
    ax_1 = subplot(3,7,[16 17 18]);
    pos = ax_1.Position;
    posnew = pos; posnew(4) = posnew(4) +.08; set(ax_1, 'Position', posnew); clear pos and posnew;
    ax_1.YLabel.Position(2) = ax_1.YLabel.Position(2) + .00005
    
    
    
    
    ax_2 = subplot(3,7,[12,13,14,19,20,21]); hold on;
    pos = ax_2.Position;
    posnew = pos; posnew(1) = posnew(1) + 0.08; posnew(4) = posnew(4) +.08; set(ax_2, 'Position', posnew); clear pos and posnew;
    
    ax_3 = subplot(3,7,[5,6,7]); hold on
    pos = ax_3.Position
    posnew = pos; posnew(1) = posnew(1) + 0.08; set(ax_3, 'Position', posnew); clear pos and posnew
end
        
        