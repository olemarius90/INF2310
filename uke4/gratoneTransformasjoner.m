clear all
close all

% Vi importerer et standardbilde fra MATLAB
f = imread('pout.tif');
[n,m] = size(f);

% Bruker vår funksjon myHist til å finne histogrammet
[h,p,c,c_n] = myHist(f);

% Vi plotter resultatene
figure(1)
subplot(141)
imshow(f,[]);
title('Bilde');
xlabel('Gråtone');
ylabel('Antall piksler');

subplot(142)
bar(h,'r');
axis tight
title('Histogram (myHist)');
xlabel('Gråtone');
ylabel('Antall piksler');

subplot(143)
bar(p);
axis tight
title('Normalisert histogram (myHist)')
xlabel('Gråtone');
ylabel('Antall piksler');

subplot(144)
imhist(f);
axis tight
title('MATLAB histogram');
xlabel('Gråtone');
ylabel('Antall piksler');

% Og kan fra resultatene trekke slutningen at vår implementasjon ser ut til
% å være riktig

%%
figure(2)
subplot(131)
bar(c)
axis tight
title('Kumulativt histogram');
xlabel('Gråtone');
ylabel('Antall piksler');

subplot(132)
bar(c_n)
axis tight
title('Kumulativt histogram normalisert');
xlabel('Gråtone');
ylabel('Antall piksler');

subplot(133)
bar(p)
title('Normalisert histogram');
xlabel('Gråtone');
ylabel('Antall piksler');

%%
%Endre "lysheten" brightness
a = 1;
b = 256/2.5;

g = f*a+b;

figure(3)
subplot(221)
imshow(f,[]);
subplot(222)
imshow(g,[]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));

%%Endre kontrasten
a = 1.5;
b = 0;

g = f*a+b;

figure(4)
subplot(221)
imshow(f,[]);
subplot(222)
imshow(g,[]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));


%%
%Middelverdi av gråtonene

u1 = mean(f(:));

u2 = 0;
for i = 1:2^8
   u2 = u2 + (i-1)*p(i);
end

u1 
u2 

%% Varians av gråtonene
v1 = var(double(f(:)));

v2 = 0;
v3 = 0;
for i = 1:2^8
    v3 = v3 + p(i)*((i-1)-u1)^2;
    v2 = v2 + (i-1)^2*p(i);
end

temp = 0;
for j = 1:2^8
    temp = temp+(j-1)*p(j);
end
v2 = v2 - temp^2;

v1
v2
v3

%% Justering av middelverdi og varians
u_ny = 240
v_ny = 150
a = sqrt(v_ny)/sqrt(v1);
b = u_ny - a*u1;

g = f*a+b;

figure(5)
subplot(221)
imshow(f,[]);
subplot(222)
imshow(g,[]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));


mean(g(:))
var(double(g(:)))


%% Logaritmiske transformasjoner

g = uint8(40*log(double(1+f)));
figure(6)
subplot(221)
imshow(f,[0 2^8]);
subplot(222)
imshow(g,[0 2^8]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));

%Show a ultrasound image as well

%% "Power-law" (gamma)-transformasjoner
gamma = 1.09;
c = 1;
g = uint8(c*double(f).^gamma);

figure(6)
subplot(221)
imshow(f,[0 2^8]);
subplot(222)
imshow(g,[0 2^8]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));

%% Bit-plan oppdeling

layer = 8;
for i = 1:n
    for j = 1:m
        out(i,j) = bitand(f(i,j),2^layer);
    end
end

out = out > 0;

figure(6)
subplot(221)
imshow(f,[0 2^8]);
subplot(222)
imshow(out,[]);
