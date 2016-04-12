function [img] = cosinusImage(M,N,u,v,A)
%COSINUSIMAGE

x = linspace(0,M-1,M);
y = linspace(0,N-1,N);
img = zeros(N,M);

for i = 1:length(x)
    for j = 1:length(y)
        img(i,j) = A*cos(2*pi*(u.*x(i)/M+v.*y(j)/N));
    end
end

end

