function refTwist = getRefTwist(d1,tangent,refTwist)
[ne,~] = size(d1);
refTwist = zeros(ne+1,1);
for kk = 2:ne
    u0 = d1(kk-1,:);
    u1 = d1(kk,:);
    t0 = tangent(kk-1,:);
    t1 = tangent(kk,:);
    
    % Method 1
%     refTwist(kk) = computeReferenceTwist(u0,u1,t0,t1);
    
    % Method 2
    ut = parallel_transport(u0,t0,t1);
    ut = rotateAxisAngle(ut,t1,refTwist(kk));
    sgnAngle = signedAngle(ut,u1,t1);
    refTwist(kk) = refTwist(kk)+sgnAngle;
    
end
end

function vNew = rotateAxisAngle(v,z,theta)
if theta == 0
    vNew = v;
else
    vNew = cos(theta)*v+sin(theta)*cross(z,v)+dot(z,v)*(1-c)*z;
end
end