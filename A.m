function y = A(u,r)
% A: [ (1 - r)z^(-1) ] / [ 1 + rz^(-4) ]
% u: gerjesztés, r: paraméter, y: kimenet
% y(n)=-ry(n-4)+(1-r)u(n-1)
    y=zeros(1, length(u));
    yBuf=zeros(1,3);
    for n=2:length(u)
        y(n)=-r*yBuf(3)+(1-r)*u(n-1);
        yBuf=[y(n) yBuf(1:end-1)];
    end
end

