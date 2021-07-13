function [Remainmask,vID] = XYR_RandROISurface(mask,slm,size)
% slm.t       = l x v matrix of data, v=#vertices; the first row
%               slm.t(1,:) is used for the clusters, and the other 
%               rows are used to calculate cluster resels if k>1. See
%               SurfStatF for the precise definition of the extra rows.
% slm.tri     = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% Remainmask  = Remained Mask for next calculation
% vID         = vertex in the generated qwqcluster

edge = SurfStatEdg(slm);
maskID = find(mask==1);
% Mask Edge to ensure all edges are inside
L1 = ismember(edge(:,1),maskID);
L2 = ismember(edge(:,2),maskID);
edge = edge((L1+L2)==2,:);

% Randomize one vertex from mask
vID = maskID(randi(length(maskID),1));
% Found neighbor Vertice around the Vertex
while 1
    surV = vID; % Define surV as this round of surround vertices
    for i=1:length(surV)
        % for single vertex in surround Vertices, found surround vertices
        NeighVID = findNeighbor(edge,surV(i)); 
        % add surround vertices in cluster, see if size exceeds
        if (length(unique(cat(1,NeighVID,vID)))>=size)
            % if size exceeds, end all loop 
            Endvalue=1;
            break;
        else
            Endvalue=0;
        end
        % if size doesnt exceed, add surround vertices in cluster
        vID = unique(cat(1,vID,NeighVID));
    end
    if (Endvalue==1)
        break;
    end
end
% Including part of NeighVIDs in vID
for j=1:length(NeighVID)
   if (length(unique(cat(1,NeighVID(j),vID)))>size)
       break
   end
   vID = unique(cat(1,vID,NeighVID(j)));
end

% Get RemainMask
otter=vID;
% Get Surrounding Vertices 
for k=1:length(vID)
   surV=findNeighbor(edge,vID(k));
   otter=unique(cat(1,otter,surV));
end
maskID = setdiff(maskID,otter);
Remainmask = zeros(1,length(mask));
Remainmask(1,maskID) = 1;
end


function vertID = findNeighbor(edge,vID)
vertID = unique(cat(1,edge(edge(:,1)==vID,2),edge(edge(:,2)==vID,1)));
end