function [allpermno] = XYR_PermuCluster(giiFile,surfFile,permno,outPre)
% prepare slm
tex=gifti(giiFile);surf=gifti(surfFile);
slm.t=tex.cdata';
slm.tri=surf.faces;
mask=ones(size(slm.t));

% get Cluster Number and Vertex Number in each Cluster
[ ~, clus, ~ ] = SurfStatPeakClus(slm, mask, 1);
allpermno=zeros(length(slm.t),permno);
for k=1:permno
    % single permutation
permuRes=zeros(size(slm.t));
for i=1:length(clus.clusid)
    % get one Cluster
    [mask,vID] = XYR_RandROISurface(mask,slm,clus.nverts(i));
    permuRes(1,vID)=1;
end
allpermno(:,k)=permuRes;
end
save(outPre,'allpermno');
end
