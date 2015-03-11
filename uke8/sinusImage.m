function [img] = sinusImage(M,N,u,v,A)
%SINUSIMAGE Summary of this function goes here
%   Detailed explanation goes here
x = linspace(0,M-1,M);
y = linspace(0,N-1,N);
img = zeros(N,M);

for i = 1:length(x)
    for j = 1:length(y)
        img(i,j) = A*sin(2*pi*(u.*x(i)/M+v.*y(j)/N));
    end
end

end

