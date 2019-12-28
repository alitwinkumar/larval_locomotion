function surfaces(P,color,alpha,lines)

k = boundary(P,.6);
if size(P,2) == 3
trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor',color,'FaceAlpha',alpha,'EdgeAlpha',alpha,'EdgeColor',color,'LineStyle','-')

else
fill(P(k,1),P(k,2),color,'FaceAlpha',alpha);
end
    
end
