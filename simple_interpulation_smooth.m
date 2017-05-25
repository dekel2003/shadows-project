
%% Load circles
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.5);
I1 = imresize(I1, 0.5);

[imageSizeY,imageSizeX] = size(I0);
grad_op    = gradient_operator_on_grid(I1); % 2m x m sparse gradient op
laplace_op = laplace_operator_on_grid(I1);
laplace_op_on_vector_field = [laplace_op zeros(size(laplace_op)); zeros(size(laplace_op)) laplace_op];

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

%% calculate the vector field by the simple interpolation:

I0 = sparse(I0(:));
I1 = sparse(I1(:));

n = numel(f2v*I0);

%% gradient descent
Av = get_vertices_areas(imageSizeY,imageSizeX); % areas of vertices
Af = sparse(1:m,1:m,ones(m,1)); % areas of faces
Af = [Af zeros(size(Af)) ; zeros(size(Af)) Af];
alpha = 0.5; % smoothness constant

% I0 = f2v * I0;
% I1 = f2v * I1;

v0 = zeros(2*m,1);

% this term is functional map C multiplied by given I
C = @(t,v,I) f2v * I + t * f2v * (vector_fields_multiplication( v, grad_op * I ));

E = @(v) (f2v * I1 - C(1,v,I0))' * Av * (f2v * I1 - C(1,v,I0)) + ... % fidelity term
    v' * Af * (eye(2*m) + alpha * laplace_op_on_vector_field) * v; % regularization term
grad_E = @(v) -2 * (vector_fields_multiplication( grad_op * I0, speye(2*m) ))' * f2v' * Av * (f2v * I1 - C(1,v,I0)) + ...
    2 * Af * (eye(2*m) + alpha * laplace_op_on_vector_field) * v;



alpha0 = 1;
beta = 0.9;
sigma = 0.005;
epsilon = 10e-5;
[ min_v, min_E ] = gradient_descent( E, grad_E, v0, alpha0, beta, sigma, epsilon );


%% results
% id = speye(size(v));
% 
% phi = @(t) (id + t*v);
I = @(t) (C(t,min_v,I0));

h = figure;
for t = 0:0.2:1
    subplot(2,3,1+5*t);
    curr_I = reshape(I(t),imageSizeY,imageSizeX);
    imshow(full(curr_I),[]);
    header = sprintf('t=%0.1f',t);
    title(header);
end

fileName = sprintf('results/1_simple_interpulation.jpg');
saveas(h,fileName);

