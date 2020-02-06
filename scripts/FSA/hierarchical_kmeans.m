function [IDx, center] = hierarchical_kmeans(Xs, K, MaxIter)

Ksmall = 50;
if(K > Ksmall)
    disp('Doing big Kmeans with 50 clusters')
    [IDx]=kmeans(Xs, Ksmall,'MaxIter',MaxIter);
    center = [];
    newIDx = IDx*0;
    for (i = 1:length(unique(IDx)))
        disp(strcat('Doing smaller kmeans on small cluster_',num2str(i)));
        sub_Xs = find(IDx == i);
        [sub_IDx, sub_center]=kmeans(Xs(sub_Xs,:), round(K/Ksmall),'MaxIter',MaxIter);
        center = [center;sub_center];
        newIDx(sub_Xs) = sub_IDx + (i-1)*round(K/Ksmall);
    end
end
IDx = newIDx;

 
 
   