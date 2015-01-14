function [h,p,c,c_n] = myHist(img)
%MYHIST Calculates the histogram of an image

[n,m] = size(img);
h = zeros(1,2^8);
p = zeros(1,2^8);
c = zeros(1,2^8);

for i = 1:2^8
    h(i) = sum(sum(uint8(img == i-1)));
end

p = h/(n*m);

for i = 2:2^8
    c(i) = h(i) + c(i-1);
end

c_n = c/(n*m);
end

