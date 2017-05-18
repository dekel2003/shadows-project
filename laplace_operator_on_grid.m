function result = laplace_operator_on_grid(I)
% get 2D grayscale image I (or just 2D grid..)
% and return a laplace matrix
% which its kernel is the following:
%   0   1   0
%   1  -4   1
%   0   1   0

[R,C] = size(I);
m = R * C;

result = eye(m) * -4;

% fill rows laplase differences
i = 1:m-1;
result(sub2ind(size(result), i, i+1)) = 1;
result(sub2ind(size(result), i+1, i)) = 1;

% fill rows laplase differences
i = 1:m-R;
result(sub2ind(size(result), i, i+R)) = 1;
result(sub2ind(size(result), i+R, i)) = 1;

result = sparse(result/8);

%% commented: old sparse matrix direct laplace calculation
%  u = sparse(size(v,1)+2,size(v,1)+2);
%  u(2:end-1,2:end-1) = v;
%  ii = 2:size(u,1)-1; jj = 2:size(u,2)-1;
%  result = - 4 * u(ii,jj) + u(ii+1,jj) + u(ii-1,jj) + u(ii,jj+1) + u(ii,jj-1);
end
