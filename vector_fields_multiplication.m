function [ result ] = vector_fields_multiplication( v )
% first vector is interpret as a m x 2m matrix
% u is 2m x k matrix
% result - m x k matrix

m = length(v) / 2;
% result = sparse([v(1:m) zeros(m, m-1) v(m+1:2*m) zeros(m, m-1)]) * u;

% result = v(1:m) .* u(1:m) + v(m+1:2*m) .* u(m+1:2*m);
% result = zeros(m,size(u,2));
% v_ = reshape(v, m, 2);
% for i = 1:size(u,2)
%     result(:,i) = sum(v_ .* reshape(u(:,i),m,2) ,2);
% end

% result = [spdiags( v(1:m)) spdiags( v(m+1:2*m))];
result = [sparse(1:m,1:m, v(1:m), m, m) sparse(1:m,1:m, v(m+1:2*m), m, m)];


end

%% test
% function test()
% a = [1 1 1 1 1 1]';
% b = [1 1 1 1 1 1]';
% vector_fields_multiplication( a, b )
% end