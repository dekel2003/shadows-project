function [ result ] = gradient_operator_on_grid( I )
% get 2D grayscale image I (or just 2D grid..)
% and return a gradient matrix
% x:             y:
%  0   0  0      0   0   0
%  0  -1  1      0  -1   0
%  0   0  0      0   1   0


[R,C] = size(I);

[Px,Py] = meshgrid(1:C,1:R);


m = R * C;


result_x = zeros(m);
result_y = zeros(m);


Yt = sub2ind(size(I), Py(2:R,:)  , Px(2:R,:)  );
Ys = sub2ind(size(I), Py(1:R-1,:), Px(1:R-1,:));
Y = [Ys(:) Yt(:)];


Xt = sub2ind(size(I), Py(:,2:C)  , Px(:,2:C)  );
Xs = sub2ind(size(I), Py(:,1:C-1), Px(:,1:C-1));
X = [Xs(:) Xt(:)];


% fill cols gradient differences x
i = m-R+1:m;
result_x(sub2ind(size(result_x), i, i)) =    -1;
i = 1:m-R;
result_x(sub2ind(size(result_x), i, X(i,1)')) =   -1;
result_x(sub2ind(size(result_x), i, X(i,2)')) =    1;

% result_x(sub2ind(size(result_x), j, X(j-R,1)')) =  1;

% fill rows gradient differences y

result_y(sub2ind(size(result_y), Y(:,1),   Y(:,1))) =  -1;
result_y(sub2ind(size(result_y), Y(:,1),   Y(:,2))) =   1;
i = R:R:m;
result_y(sub2ind(size(result_y), i,   i)) =   -1;

result = sparse([result_x;result_y]/2);
end


