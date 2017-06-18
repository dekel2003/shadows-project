
%% Load circles
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.1);
I1 = imresize(I1, 0.1);

% Af = sparse(1:m,1:m,ones(m,1)); % areas of faces

[imageSizeY,imageSizeX] = size(I0);
grad_op    = gradient_operator_on_grid(I1); % 2m x m sparse gradient op
laplace_op = laplace_operator_on_grid(I1);
% laplace_op = grad_op'*Af*grad_op;
laplace_op_on_vector_field = [laplace_op zeros(size(laplace_op)); zeros(size(laplace_op)) laplace_op];

show_vec_img = @(img) imshow(reshape(full(img),imageSizeY,imageSizeX),[]);

%% get pictures params
m = numel(I0);
f2v = faces2vertices(I0);
I0 = double(I0);
I1 = double(I1);

%% show images
subplot(1,2,1);
% imshow(full(reshape(I0,imageSizeY,imageSizeX)));
imshow(I0);
subplot(1,2,2);
% imshow(full(reshape(I1,imageSizeY,imageSizeX)));
imshow(I1);

%% calculate images as vectors of vertices

I0 = sparse(I0(:));
I1 = sparse(I1(:));

n = numel(f2v*I0);

%% show gradientx, gradienty and laplacian:
figure;

grad_res = grad_op * I0;
subplot(1,3,1);
show_vec_img(grad_res(1:m));
title('grad x');

subplot(1,3,2);
show_vec_img(grad_res(m+1:2*m));
title('grad y');

subplot(1,3,3);
show_vec_img(laplace_op * I0);
title('laplace op');

%% gradient computation
Av = get_vertices_areas(imageSizeY,imageSizeX); % areas of vertices
Af1 = sparse(1:m,1:m,ones(m,1)); % areas of faces
Af = [Af1 zeros(size(Af1)) ; zeros(size(Af1)) Af1];
alpha = 0.5; % smoothness constant
lambda = 1;

% I0 = f2v * I0;
% I1 = f2v * I1;

v0 = zeros(2*m,1)/15;

CI = @(v) I0 + (vector_fields_multiplication(v)* grad_op * I0);

curr_diff = @(v) (I1 - CI(v));

E = @(v) curr_diff(v)' * Af1 * curr_diff(v) + ... % fidelity term
    lambda * v' * Af * (speye(2*m) + alpha * laplace_op_on_vector_field) * v; % regularization term

grad_E = @(v) -2 * (vector_fields_multiplication( grad_op * I0))' * Af1 * curr_diff(v) + ...
    lambda * 2 * Af * (speye(2*m) + alpha * laplace_op_on_vector_field) * v;

params.E = E;
params.grad_E = grad_E;
target = @(v) build_target_function(v, params);


%% check correctness of gradient of E

options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
% options = optimoptions('fminunc','CheckGradients',true,'SpecifyObjectiveGradient',true);
[min_v,min_E] = fminunc(target,v0,options);
% [min_v,min_E] = fsolve(target,v0,options);

% [x, fval] = fsolve(@Energy, v0, options);
% par.EPSILON = 1e-7;
% gnum = numdiff(E,v0,par);
% plot(abs(gnum-grad_E(v0)),'*r');

%% gradient descent
% 
% alpha0 = 1.1;
% beta = 0.8;
% sigma = 0.00001;
% epsilon = 5e-1;
% [ min_v, min_E ] = gradient_descent( E, grad_E, v0, alpha0, beta, sigma, epsilon );


%% results
% id = speye(size(v));
% 
% phi = @(t) (id + t*v);
% I = @(t) (C(t,min_v,I0));
% min_v = min_v * 100;
% min_v(1:m) = min_v(1:m) * imageSizeX;
% min_v(m+1:2*m) = min_v(m+1:2*m) * imageSizeY;
I = @(t) I0 + t * (vector_fields_multiplication( min_v, grad_op * I0 ));
% I = @(t) t * (vector_fields_multiplication( min_v, grad_op * I0 ));
h = figure;
for t = 0:0.2:1
    subplot(2,3,1+5*t);
    curr_I = reshape(I(t),imageSizeY,imageSizeX);
    imshow(full(curr_I),[0,1]);
    header = sprintf('t=%0.1f',t);
    title(header);
end

%%
[x,y] = meshgrid(linspace(0,1,imageSizeY),linspace(0,1,imageSizeX));
quiver(x(:),y(:),min_v(1:m),min_v(m+1:2*m));
%% save results

fileName = sprintf('results/1_simple_interpulation.jpg');
saveas(h,fileName);

