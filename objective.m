function [J,dJ] = objective(x,Q,f)

J = 0.5*x'*Q*x + f'*x;

if nargout > 1
dJ = Q*x + f;
end