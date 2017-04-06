function result = laplace(v)
 u = sparse(size(v,1)+2,size(v,1)+2);
 u(2:end-1,2:end-1) = v;
 ii = 2:size(u,1)-1; jj = 2:size(u,2)-1;
 result = - 4 * u(ii,jj) + u(ii+1,jj) + u(ii-1,jj) + u(ii,jj+1) + u(ii,jj-1);
end
