function [ f2v ] = faces2vertices(I)

[R,C] = size(I);
m = R*C;
n = (R + 1) * (C + 1);
f2v = speye(n,m);

i = 1:m-1;
f2v(sub2ind(size(f2v), i, i+1)) = 1;
f2v(sub2ind(size(f2v), i+1, i)) = 1;
i = 1:m-R;
f2v(sub2ind(size(f2v), i, i+R)) = 1;
f2v(sub2ind(size(f2v), i+R, i)) = 1;
end

