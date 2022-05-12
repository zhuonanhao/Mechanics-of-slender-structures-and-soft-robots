function [m1,m2] = computeMaterialDirectors(d1,d2,theta)

ne = numel(theta);
m1 = zeros(ne,3);
m2 = zeros(ne,3);

for kk = 1:ne
    cs = cos(theta(kk));
    ss = sin(theta(kk));
    d1_l = d1(kk,:);
    d2_l = d2(kk,:);
    m1_l = cs*d1_l+ss*d2_l;
    m1(kk,:) = m1_l/norm(m1_l);
    m2_l = -ss*d1_l+cs*d2_l;
    m2(kk,:) = m2_l/norm(m2_l);    
end
end