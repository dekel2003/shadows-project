function [ result ] = vector_fields_multiplication( v, u )
% first vector is interpret as a m x 2m matrix
% u is 2m x k matrix
% result - m x k matrix

m = length(v) / 2;
% result = sparse([v(1:m) zeros(m, m-1) v(m+1:2*m) zeros(m, m-1)]) * u;

% result = v(1:m) .* u(1:m) + v(m+1:2*m) .* u(m+1:2*m);

result = sum(repmat(reshape(v, m, []),1,size(u,2)) .* reshape(u, m, []) ,2);

end

%% test
% function test()
% a = [1 1 1 1 1 1]';
% b = [1 1 1 1 1 1]';
% vector_fields_multiplication( a, b )
% end