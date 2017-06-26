function result = laplace_operator_on_grid(I)
% get 2D grayscale image I (or just 2D grid..)
% and return a laplace matrix
% which its kernel is the following:
%   1  -2   1
%  -2   4  -2
%   1  -2   1

[R,C] = size(I);
[Px,Py] = meshgrid(1:C,1:R);
m = R * C;

Yt = sub2ind(size(I), Py(2:R,:)  , Px(2:R,:)  );
Ys = sub2ind(size(I), Py(1:R-1,:), Px(1:R-1,:));
Y = [Ys(:) Yt(:)];


Xt = sub2ind(size(I), Py(:,2:C)  , Px(:,2:C)  );
Xs = sub2ind(size(I), Py(:,1:C-1), Px(:,1:C-1));
X = [Xs(:) Xt(:)];


result = eye(m) * 4;

% fill rows laplase differences
% i = 1:m-1;
result(sub2ind(size(result), Y(:,1), Y(:,2))) = -2;
result(sub2ind(size(result), Y(:,2), Y(:,1))) = -2;

% result(sub2ind(size(result), i+1, i)) = -1;

% fill rows laplase differences
% i = 1:m-R;
result(sub2ind(size(result), X(:,1), X(:,2))) = -2;
result(sub2ind(size(result), X(:,2), X(:,1))) = -2;
% result(sub2ind(size(result), i, i+R)) = -1;
% result(sub2ind(size(result), i+R, i)) = -1;


Yt = sub2ind(size(I), Py(2:R , 2:C)  , Px(2:R , 2:C)  );
Ys = sub2ind(size(I), Py(1:R-1,1:C-1), Px(1:R-1,1:C-1));
Y = [Ys(:) Yt(:)];


Xt = sub2ind(size(I), Py(1:R-1,2:C)  , Px(1:R-1,2:C));
Xs = sub2ind(size(I), Py(2:R,1:C-1), Px(2:R,1:C-1));
X = [Xs(:) Xt(:)];

result(sub2ind(size(result), X(:,1), X(:,2))) = 1;
result(sub2ind(size(result), X(:,2), X(:,1))) = 1;
result(sub2ind(size(result), Y(:,1), Y(:,2))) = 1;
result(sub2ind(size(result), Y(:,2), Y(:,1))) = 1;

result = sparse(result/16);

%% commented: old sparse matrix direct laplace calculation
%  u = sparse(size(v,1)+2,size(v,1)+2);
%  u(2:end-1,2:end-1) = v;
%  ii = 2:size(u,1)-1; jj = 2:size(u,2)-1;
%  result = - 4 * u(ii,jj) + u(ii+1,jj) + u(ii-1,jj) + u(ii,jj+1) + u(ii,jj-1);
end
