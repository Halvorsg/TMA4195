function [rowMap,colMap] = get_node_mapping_matrix(tri)
I = ones(1,3);
colMap = [I*tri(1),I*tri(2),I*tri(3)];
rowMap = [tri,tri,tri];
end