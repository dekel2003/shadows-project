function [ result ] = vector_fields_multiplication( v, u )
% first vector is interpret as a m x 2m matrix
% u is 2m x k matrix
% result - m x k matrix

m = length(v) / 2;
result = sparse([v(1:m) zeros(m, m-1) v(m+1:2*m) zeros(m, m-1)]) * u;

end

