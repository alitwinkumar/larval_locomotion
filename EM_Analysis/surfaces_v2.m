function surfaces_v2(mesh_obj,color,alpha, dims,axes)

if dims == 3
    DT = triangulation(mesh_obj.face,mesh_obj.vert)
    trisurf(DT,'Facecolor',color,'FaceAlpha',alpha,'EdgeAlpha',alpha,'EdgeColor',color,'LineStyle','-')

elseif dims == 2
    if ismember(axes,[1,3],'rows')
        mesh_obj.vert(:,2) = []
        DT = triangulation(mesh_obj.face,mesh_obj.vert)
        triplot(DT,color)
    elseif ismember(axes,[1,2],'rows')
        mesh_obj.vert(:,3) = []
        DT = triangulation(mesh_obj.face,mesh_obj.vert)
        triplot(DT,color)
    elseif ismember(axes,[2,3],'rows')
        mesh_obj.vert(:,1) = []
        DT = triangulation(mesh_obj.face,mesh_obj.vert)
        triplot(DT,color)
        
end
    
end



