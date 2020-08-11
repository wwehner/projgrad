function plot_QuadContour(Q,f,hfig,xvals,yvals)
% plot contours of J = 0.5 x'Qx + f'x

if nargin<4
    xvals = [-1.5 1.5];
end
if nargin<5
    yvals = [-1.5 1.5];
end
xvals = linspace(xvals(1),xvals(2),20);
yvals = linspace(yvals(1),yvals(2),20);

[X1, X2] = meshgrid(xvals, yvals);

F = 0.5*(Q(1,1)*(X1.*X1) + (Q(1,2) + Q(2,1))*(X1.*X2) + Q(2,2)*(X2.*X2)) + ...
  f(1)*X1 + f(2)*X2;

% V is contour values
lower = (min(min(F)));
upper = (max(max(F)));
V = linspace(lower,upper,15);

% plot contours
if nargin < 3
    figure;
else
    figure(hfig);
end
[c,h]=contour(X1, X2, F, V, '--', 'linewidth',2);
colormap(gca,gray)

end