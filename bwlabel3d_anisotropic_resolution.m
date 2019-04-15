function masks = bwlabel3d_anisotropic_resolution(bw)

G = graph;
nodes = zeros(1e4,1);
labels = zeros(size(bw));
idx = 1;

label1 = bwlabel(bw(:,:,1),4);
labels(:,:,1) = label1;
for i_label = 1:max(label1(:))
    G = addnode(G, 1);
    nodes(idx) = 1e4 + i_label;
    idx = idx + 1;
end
for i_slice = 2:size(bw,3) 
    label2 = bwlabel(bw(:,:,i_slice),4);
    labels(:,:,i_slice) = label2;
    for i_label = 1:max(label2(:))
        G = addnode(G, 1);
        nodes(idx) = i_slice*1e4 + i_label;
        idx = idx + 1;
        
        for j_label = 1:max(label1(:))
            if sum(label1(:) == j_label & label2(:) == i_label) > ...
                    0.1*0.5*(sum(label1(:) == j_label)+sum(label2(:) == i_label))
                G = addedge(G, find(nodes == i_slice*1e4 + i_label),...
                    find(nodes == (i_slice-1)*1e4 + j_label));
            end
        end
    end
    label1 = label2;
end

bins = conncomp(G);

masks = false(size(bw,1),size(bw,2),size(bw,3),max(bins));
for i_bin = 1:max(bins)
    idx = find(bins == i_bin);
    for i = 1:length(idx)
        masks(:,:,floor(nodes(idx(i))/1e4),i_bin) = ...
            masks(:,:,floor(nodes(idx(i))/1e4),i_bin) | ...
            labels(:,:,floor(nodes(idx(i))/1e4)) == rem(nodes(idx(i)), 1e4);
    end
end

end