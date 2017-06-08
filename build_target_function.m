function [val,grad] = build_target_function(v, params)
val = params.E(v);
if nargout > 1
    grad = params.grad_E(v);
end
end