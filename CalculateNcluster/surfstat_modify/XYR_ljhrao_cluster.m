function [ peak, clus, clusid ] = XYR_ljhrao_cluster(texdata,mask,surfFile)
slm.t=texdata';
gii=gifti(surfFile);
slm.tri=gii.faces;
data=unique(texdata);
[ peak, clus, clusid ] = SurfStatPeakClus(slm, mask, data(1));
end