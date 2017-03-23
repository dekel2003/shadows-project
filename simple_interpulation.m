
%% Load circles
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

[imageSizeY,imageSizeX] = size(I0);

%% get pictures params
I0 = faces2vertices(double(I0));
I1 = faces2vertices(double(I1));
number_of_faces = numel(I0);

%% show images
subplot(1,2,1);
% imshow(full(reshape(I0,imageSizeY,imageSizeX)));
imshow(I0);
subplot(1,2,2);
% imshow(full(reshape(I1,imageSizeY,imageSizeX)));
imshow(I1);

%% calculate data fidelity term without smothness
I0 = sparse(I0(:));
I1 = sparse(I1(:));

% grad_I = gradient(I0);
number_of_vertices = numel(I0);

%% calculate the vector field by the simple interpolation:

v = (I0' \ (I1-I0)')';
id = eye(size(v));

phi = @(t) (id + t*v);
I = @(t) (phi(t) * I0);

h = figure;
for t = 0:0.2:1
    subplot(2,3,1+5*t);
    curr_I = reshape(I(t),imageSizeY,imageSizeX);
    imshow(curr_I);
    header = sprintf('t=%0.1f',t);
    title(header);
end

fileName = sprintf('results/1_simple_interpulation.jpg');
saveas(h,fileName);
