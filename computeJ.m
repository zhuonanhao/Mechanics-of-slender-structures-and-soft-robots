% Compute J
for i=1:N
    for j=1:N
        ki = K(i);
        if i==N
            kip1 = 0; % k_{i+1}
        else
            kip1 = K(i+1);
        end

        bi = b(i);
        if i==N
            bip1 = 0; % b_{i+1}
        else
            bip1 = b(i+1);
        end

        if i==j
            J(i,j) = m(i)/dt^2 + (ki+kip1) + (bi+bip1)/dt;
        elseif j==i-1
            J(i,j) = -ki - bi/dt;
        elseif j==i+1
            J(i,j) = kip1 - bip1/dt;
        else
            J(i,j) = 0;
        end
    end
end
