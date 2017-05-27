function [ result ] = gradient_operator_on_grid( I )
% get 2D grayscale image I (or just 2D grid..)
% and return a gradient matrix
% x:             y:
%  1   0  -1      1   2   1
%  2   0  -2      0   0   0
%  1   0  -1     -1  -2  -1


[R,C] = size(I);
m = R * C;

result_x = zeros(m);
result_y = zeros(m);

i = 1:m-R;
% fill rows laplase differences
result_x(sub2ind(size(result_x), i, i+R)) = -2;
result_x(sub2ind(size(result_x), i+R, i)) = 2;

result_x(sub2ind(size(result_x), i+R, i)-1) = 1;
result_x(sub2ind(size(result_x), i+R, i)+1) = 1;
result_x(sub2ind(size(result_x), i, i+R)-1) = -1;
result_x(sub2ind(size(result_x), i, i+R)+1) = -1;

% fill rows laplase differences
i = 1:m-1;
j = 1:m-R;

result_y(sub2ind(size(result_y), i, i)+1) = -2;
result_y(sub2ind(size(result_y), i+1, i+1)-1) = 2;

result_y(sub2ind(size(result_y), j+R, j)-1) = 1;
result_y(sub2ind(size(result_y), j+R, j)+1) = -1;
result_y(sub2ind(size(result_y), j, j+R)-1) = 1;
result_y(sub2ind(size(result_y), j, j+R)+1) = -1;

result = sparse([result_x;result_y]/8);
end


