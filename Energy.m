function [E,grad_E] = Energy(v)

addpath(genpath(pwd));
load('circles.mat');

I0 = imresize(I0, 0.1);
I1 = imresize(I1, 0.1);

[imageSizeY,imageSizeX] = size(I0);
grad_op    = gradient_operator_on_grid(I1); % 2m x m sparse gradient op
laplace_op = laplace_operator_on_grid(I1);
laplace_op_on_vector_field = [laplace_op zeros(size(laplace_op)); zeros(size(laplace_op)) laplace_op];

%% get pictures params
m = numel(I0);
f2v = faces2vertices(I0);
I0 = double(I0);
I1 = double(I1);

I0 = sparse(I0(:));
I1 = sparse(I1(:));
n = numel(f2v*I0);

%% gradient computation
Av = get_vertices_areas(imageSizeY,imageSizeX); % areas of vertices
Af = sparse(1:m,1:m,ones(m,1)); % areas of faces
Af = [Af zeros(size(Af)) ; zeros(size(Af)) Af];
alpha = 1; % smoothness constant

% this term is functional map C multiplied by given I
CI = f2v * I0 + f2v * (vector_fields_multiplication( v, grad_op * I0 ));

curr_diff = (f2v * I1 - CI);

E = curr_diff' * Av * curr_diff + ... % fidelity term
    v' * Af * (speye(2*m) + alpha * laplace_op_on_vector_field) * v; % regularization term

if nargout > 1
    grad_E = -2 * (vector_fields_multiplication( grad_op * I0, speye(2*m) ))' * f2v' * Av * curr_diff + ...
        2 * Af * (speye(2*m) + alpha * laplace_op_on_vector_field) * v;
end

end

function test
%% check correctness of gradient of E
v0 = zeros(2*m,1)/15;
rng(0,'twister'); 
options = optimoptions(@fmincon, 'Algorithm', 'interior-point',... 
    'CheckGradients', 'on', 'MaxIter', 5, 'SpecifyObjectiveGradient', 'on');

[x fval exitflag output] = fmincon(@Energy, v0,[],[],[],[],[],[],@Energy,options);
% par.EPSILON = 1e-7;
% gnum = numdiff(E,v0,par);
% plot(abs(gnum-grad_E(v0)),'*r');
end
