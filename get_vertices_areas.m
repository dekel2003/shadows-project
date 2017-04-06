function [ Av ] = get_vertices_areas(rows,cols)

W = ones(rows,cols);
filter = ones(2) / 4;

A = xcorr2(W,filter);
numelA = numel(A(:));
Av = sparse(1:numelA, 1:numelA, A(:));

end