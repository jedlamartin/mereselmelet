function y = D(u,r)
% D: [ (1 - r)z^(-1)(1 + z^(-2) ] / [ 2(1 + rz^(-4)) ]
% u: gerjesztés, r: paraméter, y: kimenet
% y(n)=ry(n-4)+1/2*(1-r)u(n-1)+1/2*(1-r)u(n-3)
    y=zeros(1, length(u));
    yBuf=zeros(1,3);
    uBuf=zeros(1,2);
    for n=2:length(u)
        y(n)=r*yBuf(3)+1/2*(1-r)*u(n-1)+1/2*(1-r)*uBuf(2);
        yBuf=[y(n) yBuf(1:end-1)];
        uBuf=[u(n-1) uBuf(1)];
    end
end
