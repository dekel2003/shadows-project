
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

%% get pictures params
number_of_faces = numel(I0);
I0 = double(I0);
I1 = double(I1);

% I0 = faces2vertices(double(I0));
% I1 = faces2vertices(double(I1));

%% show images
subplot(1,2,1);
% imshow(full(reshape(I0,imageSizeY,imageSizeX)));
imshow(I0);
subplot(1,2,2);
% imshow(full(reshape(I1,imageSizeY,imageSizeX)));
imshow(I1);

%% calculate the vector field by the simple interpolation:

I0 = sparse(faces2vertices(I0(:)));
I1 = sparse(faces2vertices(I1(:)));
number_of_vertices = numel(I0);

v = (I0' \ (I1-I0)')';

%% gradient descent
Av = get_vertices_areas(imageSizeY,imageSizeX); % areas of vertices
Af = sparse(1:number_of_faces,1:number_of_faces,ones(number_of_faces,1)); % areas of faces
alpha = 0.5; % smoothness constant
ID_nxn = eye(number_of_vertices);


C = @(t,v,I) ID_nxn + t * faces2vertices([v(1:m) zeros(m, m-1) v(m+1:2*m) zeros(m, m-1)] * (grad_op * I));
E = @(v) (I1 - C(1,v,I0))' * Av * (I1 - C(1,v,I0)) + ... % fidelity term
    v' * Af * (eye(number_of_faces) + alpha * grad_op) * v; % regularization term


v = (I0' \ (I1-I0)')';
sigma = 4;


alpha = 0.0005;
last_v = v;
% prevE = v' * Av * laplace2(v) - (1/2*sigma^2) * Av * (I1 - I0 + v*I0) * (I1 - I0 + v*I0)';
prevE = v' .* Av .* laplace2(v) - (1/2*sigma^2) * Av * (I1 - I0 + v*I0) * (I1 - I0 + v*I0)';
prevE = prevE.^2;
prevE = full(sum(prevE(:)));


%% results
id = speye(size(v));

phi = @(t) (id + t*v);
I = @(t) (phi(t) * I0);

h = figure;
for t = 0:0.2:1
    subplot(2,3,1+5*t);
    curr_I = reshape(I(t),imageSizeY+1,imageSizeX+1);
    imshow(full(curr_I),[]);
    header = sprintf('t=%0.1f',t);
    title(header);
end

fileName = sprintf('results/1_simple_interpulation.jpg');
saveas(h,fileName);

