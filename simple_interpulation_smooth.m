
%% Load circles
clear;
close all;

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.2);
I1 = imresize(I1, 0.2);

[imageSizeY,imageSizeX] = size(I0);

%% get pictures params
number_of_faces = numel(I0);
I0 = faces2vertices(double(I0));
I1 = faces2vertices(double(I1));

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
number_of_vertices = numel(I0);

v = (I0' \ (I1-I0)')';

%% gradient descent
Af = sparse(1:number_of_faces,1:number_of_faces,ones(number_of_faces,1));
v = (I0' \ (I1-I0)')';
sigma = 2;
Av = get_vertices_areas(imageSizeY,imageSizeX);

alpha = 0.0005;
last_v = v;
prevE = v' * Av * laplace2(v) - (1/2*sigma^2) * Av * (I1 - I0 + v*I0) * (I1 - I0 + v*I0)';
prevE = prevE.^2;
prevE = full(sum(prevE(:)));

while true
    dEdv = (Av * laplace2(v) - (1/sigma^2) * Av * I0 * (I1 - I0 + v*I0)');
    v = v - alpha * dEdv;

    E = v' * Av * laplace2(v) - (1/2*sigma^2) * Av * (I1 - I0 + v*I0) * (I1 - I0 + v*I0)';
    E = E.^2;
    E = full(sum(E(:)));
    
    if E < prevE
%         alpha = alpha * 1.08;
        last_v = v;
    else
        alpha = alpha / 1.1;
        v = last_v;
    end
    disp ((prevE - E) / prevE);
%     prevE = E;
    
%     Err = full(sum(sum(dEdv.^2)));
    if abs(prevE - E) / prevE < sqrt(eps)
        break;
    end
    prevE = E;
end

v = last_v;

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

