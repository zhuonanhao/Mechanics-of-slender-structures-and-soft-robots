function kappa = getkappa(q,m1,m2)
nv = (numel(q)+1)/4;
ne = nv-1;
kappa = zeros(nv,2);
for kk = 2:ne
    node0 = q(4*(kk-2)+1:4*(kk-2)+3);
    node1 = q(4*(kk-1)+1:4*(kk-1)+3);
    node2 = q(4*(kk-0)+1:4*(kk-0)+3);
    m1e = m1(kk-1,:);
    m2e = m2(kk-1,:);
    m1f = m1(kk,:);
    m2f = m2(kk,:);
    kappaL = computekappa(node0,node1,node2,m1e,m2e,m1f,m2f);
    kappa(kk,1) = kappaL(1);
    kappa(kk,2) = kappaL(2);
end
end