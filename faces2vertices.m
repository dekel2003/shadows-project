function [ f2v ] = faces2vertices(I)

[R,C] = size(I);
m = R*C;
n = (R + 1) * (C + 1);
f2v = sparse(n,m);

i = 1:m-R;
j = 1:m-1;

f2v(sub2ind(size(f2v), j, j) + floor((j-1)/R)) = 1;
f2v(sub2ind(size(f2v), j, j+1) + floor((j-1)/R)) = 1;
f2v(sub2ind(size(f2v), i+R, i) + floor((i-1)/R)) = 1;
f2v(sub2ind(size(f2v), i, i+R) + floor((i-1)/R)) = 1;

% i = 1:m-R;
% f2v(sub2ind(size(f2v), i, i+R)) = 1;
% f2v(sub2ind(size(f2v), i+R, i)) = 1;

f2v = f2v/4;
end

