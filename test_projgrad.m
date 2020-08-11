clc; clear; 

% test Projected Gradient Descent
% min   J   = x'Qx + f'x
% s.t.  Ax <= b

n = 2;
Q = randn(2); 
Q = .5*(Q'*Q); % make it P.D.
f  = [-1; .1];

x0 = 2*randn(2,1);

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
    y_poly = 1.5*sin(theta_poly-pi/10);% rotate the y coord

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

% Plot the iterates and search directions
plot(data.x(1,:),data.x(2,:),'ko','MarkerFaceColor','k')
plot(data.x(1,:)-data.g(1,:),data.x(2,:)-data.g(2,:),'ro','MarkerFaceColor','r')
quiver(data.x(1,:),data.x(2,:),-data.g(1,:),-data.g(2,:),0)
