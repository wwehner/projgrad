# projgrad
Matlab implementation of projected gradient descent

Two versions of projected gradient descent.

the first works well (prograd.m), and the second (projgrad_algo2.m) is shown to fail in certain cases (see the doc)

projgrad.m - main algorithm
test_projgrad.m - demonstrates the algorithm

projgrad_algo2.m - alternate proj grad algo that fails
test_projgrad_algo2.m - demonstrates the failure

plotregion.m and plot_QuadContour.m are used for plotting the objective function 
and feasible region

Requires matlab optimization toolbox for quadprog.
