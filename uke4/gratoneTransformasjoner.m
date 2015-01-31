% Gråtonetrasformasjoner
%
% Eksempelkode skrevet for å illustrere temaene gått gjennom i forelesning
% 4 i faget INF2310 - Digital bildebehandling ved Ifi UiO.
%
% Skrevet av Ole Marius Hoel Rindal
% Vennligst gi beskjed til olemarius@olemarius.net om du finner feil eller
% mangler i koden.

%%  Histogram
%   Et histogram og da menes her et gråtonehistgram til et bilde er en
%   diskret funksjon som viser antall piksler med hver gråtone.
%
%   Altså vil et histogram h være slik at h(i) indikerer antall piksler i
%   bildet med pikselverdi i.
%
%   I den følgende koden vil vi anta at vi har bilder med 256 forskjellige
%   gråtoner.
clear all
close all

% Vi importerer et standardbilde fra MATLAB
f = imread('pout.tif');
% Her er det også fint å eksprimentere med egne bilder
% f = rgb2gray(imread('dittBilde'));
[n,m] = size(f);

% Bruker vår funksjon myHist til å finne histogrammet
% Se myHist.m filen for detaljer på hvordan man kan finne histogrammet.
[p,h,c,c_n] = myHist(f);

% Vi plotter resultatene
% Først plotter vi kun bildet
figure(1)
subplot(141)
imshow(f,[0 255]);
title('Bilde');

% Deretter plotter vi histogrammet vi får fra myHist()
subplot(142)
bar(h,'r');
axis tight
title('Histogram (myHist)');
xlabel('Gratone');
ylabel('Antall piksler');

% samt det normaliserte histogrammet også returnert fra myHist()
subplot(143)
bar(p);
axis tight
title('Normalisert histogram (myHist)')
xlabel('Gratone');
ylabel('Antall piksler');

% Til slutt finner vi histogrammet med MATLAB's egen imhist() funksjon
subplot(144)
imhist(f);
axis tight
title('MATLAB histogram');
xlabel('Gratone');
ylabel('Antall piksler');

% Og kan fra resultatene trekke slutningen at vår implementasjon ser ut til
% å være riktig

%%  Kumulativt histogram
%
%   Det kumulative histogrammet, c, indikerer hvor mange piksler som har
%   gråtone mindre enn eller lik gråtone i.

% Vi plotter det kumulative histogrammet
figure(2)
subplot(131)
bar(c)
axis tight
title('Kumulativt histogram');
xlabel('Gratone');
ylabel('Antall piksler');

% samt det normaliserte kumulative histogrammet
subplot(132)
bar(c_n)
axis tight
title('Kumulativt histogram normalisert');
xlabel('Gratone');
ylabel('Antall piksler');

% plotter også det vanlige normaliserte histogrammet som referanse
subplot(133)
bar(p)
title('Normalisert histogram');
xlabel('Gratone');
ylabel('Antall piksler');

%%  Gråtonetransformasjoner 
%   Lineære transformasjoner.
%   Vi skal nå begynne å se på gråtonetransformasjoner, og ser først på
%   lineære gråtonetrasformasjoner. Altså hvor g = f*a+b.
%   Den første gråtonetransformasjonen vi skal se på er hvordan vi kan endre
%   "lysheten" (brightness) i bildet.
%   For å endre lysheten i bildet legger vi til en kostant, b, til alle
%   pikselverdiene. Dette endrer middelverdien i bildet. Forsøk med flere
%   verdier for b og se forskjellen i bildene.
a = 1;
b = 256/2.5;

g = f*a+b;

figure(3)
subplot(221)
imshow(f,[0 255]);
title('Orginal');
subplot(222)
imshow(g,[0 255]);
title('Endret "lysheten"');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(g));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histgram etter transformasjon')

%%  Endre kontrasten
%   For å endre kontrasten multipliserer vi alle pikselverdiene i
%   orginalbildet med en faktor a. Forsøk med forskjellige verdier for a.
a = 1.5;
b = 0;

g = uint8(double(f)*a+b);

figure(4)
subplot(221)
imshow(f,[0 255]);
title('Orginal');
subplot(222)
imshow(g,[0 255]);
title('Endret kontrast');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(g));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histgram etter transformasjon')

% NB! Prøv feks å multipliser med a = 2. Hva skjer?
% Vi ser da effekten av klipping, siden vi multipliserte med et tall som
% gjorde at vi fikk verdier utenfor det støttede intervallet.

%%  Fra gråtonenivå [f1,f2] til [g1,g2]
%   Endrer intervallet [f1,f2] til [g1,g2]. I dette eksempelet endret vi
%   slik at vi utnytter hele intervallet med gråtoner, fra 0 til 255.

%   Finner max og min i orignal bildet
f1 = min(f(:));
f2 = max(f(:));

%   Bestemmer max og min i det nye intervallet
g1 = 0;
g2 = 255;

%   Dette er i grunnen kun en normalisering.
g = g1 + ((g2-g1)/(f2-f1))*(f-f1);

%   Viser frem resultatet
figure(4)
subplot(221)
imshow(f,[0 255]);
title('Orginal');
subplot(222)
imshow(g,[0 255]);
title('Bruker hele gratoneskalaen');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(g));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histgram etter transformasjon')



%%  Standarisering av bildet
%   Her skal vi sørge for at alle bildene i en serie er statistisk like
%   (1. orden). Dette gjør vi ved å justere middelverdien og variansen til
%   bildet ved hjelp av en lineær gråtonetransform.


%   Først finner vi middelverdi av gråtonene

%   Bruker enkel MATLAB mean til å finne middelverdien
u1 = mean(f(:));

%   Kan også finne middelverdien ved å hjelp av det normaliserte
%   histogrammet.
u2 = 0;
for i = 1:2^8
   u2 = u2 + (i-1)*p(i);
end

%   Verifiserer at det er likt
u1 
u2 

%% Varians av gråtonene
%   Enkel MATLAB varians
v1 = var(double(f(:)));

%   Variansen kan også finnes ved hjelp av histogrammet
%   Viser tre forskjellige implementasjoner
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

% Verifiserer at de er like
v1
v2
v3

%%  Justering av middelverdi og varians
%   Videre kan vi ved hjelp av en lineær gråtonetransformasjon justere
%   middelverdien og variansen

%   Bestemmer ny middelverdi og varians
u_ny = 200
v_ny = 120

%   Se foil 30 for utledning
a = sqrt(v_ny)/sqrt(v1);
b = u_ny - a*u1;

%   Bruker en lineær grtonetrasformasjon med a og b 
g = f*a+b;

%   Plotter resultatet
figure(5)
subplot(221)
imshow(f,[0 255]);
title('Orginal');
subplot(222)
imshow(g,[0 255]);
title('Ny middelverdi og varians');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(g));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histgram etter transformasjon')


% Middelverdien og variansen ble som ønsket
mean(g(:))
var(double(g(:)))


%%  Ikke-lineære trasformer 
%   Vi skal nå se på ikke-lineære trasformasjoner, og begynnver ved å se på
%   logaritmiske transformasjoner.

x = 0:255;

r = 0;
c = 46;

figure(66)
clf
plot(x/max(x))
hold on
plot(c*log(x+r)/max(c*log(x+r)),'r')
plot(x.^(1/2)/max(x.^(1/2)),'g');
plot(x.^(2)/max(x.^(2)),'y');
g = uint8(40*log(double(1+f)));

figure(6)
subplot(221)
imshow(f,[0 255]);
subplot(222)
imshow(g,[0 255]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));

%%
load UltrasoundImageSim
figure(7)
subplot(221)
imshow(us_image,[]);
title('Simulated ultrasound image');
colorbar
subplot(222)
imshow(20*log10(us_image),[])
title('Log transform of ultrasound image');
colorbar
subplot(223)
bar(myHist(uint8(us_image)));
subplot(224)
bar(myHist(uint8(20*log10(us_image))));

%% "Power-law" (gamma)-transformasjoner
gamma = 1.09;
c = 1;
g = uint8(c*double(f).^gamma);

figure(8)
subplot(221)
imshow(f,[0 2^8]);
subplot(222)
imshow(g,[0 2^8]);
subplot(223)
bar(myHist(f));
subplot(224)
bar(myHist(g));

%%  Stykkvis lineær mapping
%   Verifiser at dette er en korrekt implementasjon
img = imread('tm4.png');

a_1 = 1.5;
b_1 = 0;
a_2 = 2;
b_2 = 0;

limit = 30;

idx = img > 0 & img < 20;
idx_2 = img > limit+1 & img < 255;

img_g = img;

img_g(idx) = img(idx)*a_1+b_1;
img_g(idx_2) = img(idx_2)*a_2+b_2;
%img_g = reshape(img_g,n,m);
%img_g = img;
figure(9)
subplot(221)
imshow(img,[0 255]);
title('Original');
colorbar
subplot(222)
imshow(img_g,[0 255])
title('Etter stykkvis linear mapping');
colorbar
subplot(223)
bar(myHist(uint8(img)));
subplot(224)
bar(myHist(uint8(img_g)));


%% Bit-plan oppdeling

layer = 7;
for i = 1:n
    for j = 1:m
        out(i,j) = bitand(f(i,j),uint8(2^layer));
    end
end

out = out > 0;

figure(10)
subplot(221)
imshow(f,[0 2^8]);
subplot(222)
imshow(out,[]);
