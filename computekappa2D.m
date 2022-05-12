%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.

function curvature = computekappa2D(xkm1, ykm1, xk, yk, xkp1, ykp1)
%
% This function returns the scalar curvature (\kappa_k) between the nodes
% [x_{k-1}, y_{k-1}], [x_k, y_k], and [x_{k+1}, and y_{k+1}].
%

curvature = 0.2e1 * tan(atan((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1);
end
