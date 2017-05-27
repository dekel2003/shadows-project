function [ min_x, min_val ] = gradient_descent( f, grad_f, x0, alpha0, beta, sigma, epsilon )

x = x0;
alpha = alpha0;
while (norm(grad_f(x)) >= epsilon)
    grad_f_x = grad_f(x);
    disp(norm(grad_f_x));
disp(['current error is: ' num2str(norm(f(x)))]);
    alpha = armijo(f, -grad_f_x , x, alpha, beta, sigma);
    x = x - alpha * grad_f_x;
%     disp([grad_f(x) x]);
end

min_x = x;

if (nargout>1)
    min_val = f(x);
end

end

function alpha_new = armijo(f, grad_f_x, x, alpha, beta, sigma)
    m = sum(grad_f_x .* normc(grad_f_x));
    t = -sigma * m + sqrt(sqrt(eps))*1000;
    alpha_new = alpha;
    if (f(x) - f(x + alpha_new*grad_f_x)) < alpha_new*t
        alpha_new = alpha_new * beta;
        disp(['step size has been changed to: ' num2str(alpha_new)]);
%     else
%         alpha_new = alpha * 1.01;
    end
end