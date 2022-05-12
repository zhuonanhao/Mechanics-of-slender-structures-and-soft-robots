function [q, d1, d2, refTwist] = DER(q0, u, d1_old, d2_old, refTwist_old)

% Global variables
global m dt unconsInd tol ScaleSolver maximum_iter
global Fg

iter = 0;
err = tol * ScaleSolver * 10;

q = q0; % Guess

while err > tol * ScaleSolver    
    iter = iter + 1;
    
    [d1, d2] = computeTimeParallel(d1_old, q0, q);
    tangent = computeTangent(q);
    refTwist = getRefTwist(d1, tangent, refTwist_old);
    theta = q(4:4:end);
    [m1,m2] = computeMaterialDirectors(d1, d2, theta);
    
    % Get forces
    [Fb, Jb] = getFb(q, m1, m2); % Bending
    [Ft, Jt] = getFt(q, refTwist); % Twisting
    [Fs, Js] = getFs(q); % Stretching
    
    F = m .* (q-q0)/dt^2 - m .* u/dt - (Fb+Ft+Fs+Fg);    
    J = diag(m)/dt^2 - (Jb+Jt+Js);
    
    F_free = F(unconsInd);
    J_free = J(unconsInd, unconsInd);
    dq = J_free \ F_free;
    
    q(unconsInd) = q(unconsInd) - dq;
    
    err = sum(abs(F_free));
    
    if iter > maximum_iter
        fprintf('Error\n');
        return
    end
end

end
