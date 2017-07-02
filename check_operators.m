
%%
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.5);
I1 = imresize(I1, 0.5);

% Af = sparse(1:m,1:m,ones(m,1)); % areas of faces

[imageSizeY,imageSizeX] = size(I0);
m = numel(I0);
grad_op    = gradient_operator_on_grid(I1); % 2m x m sparse gradient op
laplace_op = laplace_operator_on_grid(I1);

laplace_op_on_vector_field = [laplace_op zeros(size(laplace_op)); zeros(size(laplace_op)) laplace_op];
vec_img = @(img) reshape(full(img),imageSizeY,imageSizeX);

Af1 = speye(m);
laplace_op2 = grad_op'*grad_op;
%% check gradient operator in x direction
% [U,L] = eigs(grad_op(1:m,:),Af1,9,'sm');
% [~,ii] = sort(diag(L));
% L = (L(ii,ii));
% U = (U(:,ii));
% figure, plot(diag(L),'.');
% figure;
% for i=1:9
%     subplot(3,3,i);
%     imagesc(vec_img(U(:,i)));
%     t = sprintf('e.v. = %.4f',L(i,i));
%     title(t);
% end


%% check Laplace operator
[U,L] = eigs(laplace_op,Af1,9,'sm');

[~,ii] = sort(diag(L));
L = (L(ii,ii));
U = (U(:,ii));
figure, plot(diag(L),'.');
figure;
for i=1:9
    subplot(3,3,i);
    imagesc(vec_img(U(:,i)));
    t = sprintf('e.v. = %.4f',L(i,i));
    title(t);
end

% check null space
figure;
I = laplace_op * ones(m,1);
imagesc(vec_img(I));

% check symmetry
if norm(laplace_op-laplace_op','fro') ~= 0
    disp('laplace operator is not symmetric');
end





