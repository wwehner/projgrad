clc; clear; 

% test Projected Gradient Descent
% min   J   = x'Qx + f'x
% s.t.  Ax <= b

n = 2;

Q = [0.3   -0.7
    -0.7    2.0];

f  = [-1; .1];

x0 = [-1; 1];

% Feasible Set
A = [ 1   0
     -1   0
      0   1
      0  -1];
b = [1 1 1 1]';

% Alternate Feasible Domain: Define a polygon
if 1
    
    % number of polygon sides
    n_poly = 5;
    % Ax = b, defines the polygon interior
    A = zeros(n_poly,2);
    b = zeros(n_poly,1);

    % define polygon (n_poly points around a circle)
    theta_poly = (0:(2*pi)/n_poly:2*pi);
    x_poly = cos(theta_poly);
    y_poly = 1.5*sin(theta_poly-pi/10);% skew the y coord

    % segments of circle define the polygon
    for i = 1:n_poly

        % compute direction p = x2-x1
        p = [x_poly(i+1)-x_poly(i)
             y_poly(i+1)-y_poly(i)];

        % get the normal direction
        % rotate 90 degrees counter-clockwise
        a = [p(2),-p(1)];
        a = a/norm(a);

        A(i,:) = a;

        % offset of the polygon side
        b(i) = a*[x_poly(i); y_poly(i)];

    end

end

% Plot the QP
figure(1),clf, hold on, box on
    plotregion(-A,-b,[],[],[0.8,0.8,0.8]); % plot the feasible region
    plot_QuadContour(Q,f,gcf); % plot the QP contours
    set(gca,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
    axis equal tight
    drawnow

x_qp = quadprog(Q,f,A,b);

% Projected grad descent
[x_pgd,data] = projgrad(@(x)( objective(x,Q,f) ),A,b,x0);

fprintf('Solution\n  QuadProg     PGD\n')
disp([x_qp,x_pgd])


% Projected grad descent V2 
% - this algorith fails to converge if minimizer is on the edge of feasible set
%   and the feasible edge is perpendicular to the objective gradient at the minimizer.
[x_pgd,data] = projgrad_algo2(@(x)( objective(x,Q,f) ),A,b,x0);

fprintf('Solution\n  QuadProg     PGD     (approach V2)\n')
disp([x_qp,x_pgd])

% Plot the iterates and search directions
for i = 1:15
    % plot iterates
    plot(data.x(1,i),data.x(2,i),'ko','MarkerFaceColor','k')
    %plot(data.y(1,i),data.y(2,i),'ro','MarkerFaceColor','r')
    % plot the step directions (negative gradient)
    quiver(data.x(1,i),data.x(2,i),-data.g(1,i),-data.g(2,i),0)
    pause(0.2)
end
title({'Black dot = iterates','Arrow = Search Direction','Red Dot = minimizer'})

% plot minimizer
plot(x_qp(1,:),x_qp(2,:),'ro','MarkerSize',12,'Markerfacecolor','r')

