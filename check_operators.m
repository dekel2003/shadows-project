
%%
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.05);
I1 = imresize(I1, 0.05);

% Af = sparse(1:m,1:m,ones(m,1)); % areas of faces

[imageSizeY,imageSizeX] = size(I0);
m = numel(I0);
grad_op    = gradient_operator_on_grid(I1); % 2m x m sparse gradient op
laplace_op = laplace_operator_on_grid(I1);
laplace_op_on_vector_field = [laplace_op zeros(size(laplace_op)); zeros(size(laplace_op)) laplace_op];
show_vec_img = @(img) imshow(reshape(full(img),imageSizeY,imageSizeX),[]);

%% check gradient operator in x direction
[U,L] = eigs(grad_op(1:m,:),16,'sm');

figure;
for i=1:16
    subplot(4,4,i);
    show_vec_img(U(:,i));
    t = sprintf('e.v. = %.4f',L(i,i));
    title(t);
end


%% check Laplace operator
[U,L] = eigs(laplace_op,16,'sm');

figure;
for i=1:16
    subplot(4,4,i);
    show_vec_img(U(:,i));
    t = sprintf('e.v. = %.4f',L(i,i));
    title(t);
end




