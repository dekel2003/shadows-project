function [ result ] = gradient_operator_on_grid( I )
% get 2D grayscale image I (or just 2D grid..)
% and return a gradient matrix
% which its kernel is the following
%   0  -1   0
%  -1   0   1
%   0   1   0


[R,C] = size(I);
m = R * C;

result_x = zeros(m);
result_y = zeros(m);

% fill rows laplase differences
i = 1:m-1;
result_x(sub2ind(size(result_x), i, i+1)) = 1;
result_x(sub2ind(size(result_x), i+1, i)) = -1;

% fill rows laplase differences
i = 1:m-R;
result_y(sub2ind(size(result_y), i, i+R)) = 1;
result_y(sub2ind(size(result_y), i+R, i)) = -1;

result = sparse([result_x;result_y]/4);
end


