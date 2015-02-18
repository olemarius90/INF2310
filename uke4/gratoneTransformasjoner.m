% Gråtonetrasformasjoner
%
% Eksempelkode skrevet for å illustrere temaene gått gjennom i forelesning
% 4 i faget INF2310 - Digital bildebehandling ved Ifi UiO.
%
% Skrevet av Ole Marius Hoel Rindal
% Vennligst gi beskjed til olemarius@olemarius.net om du finner feil eller
% mangler i koden.
%
% NB! Veldig mye av koden er kun kode for å vise frem bilder og histogram.
% Denne koden er kun kopiert videre og bør muligens flyttes i en funksjon.

%%  Histogram
%   Et histogram, og da menes her et gråtonehistgram, til et bilde er en
%   diskret funksjon som viser antall piksler med hver gråtone.
%
%   Altså vil et histogram h være slik at h(i) indikerer antall piksler i
%   bildet med pikselverdi i.
%
%   I den følgende koden vil vi anta at vi har bilder med 256 forskjellige
%   gråtoner, fra gråtone 0 til gråtone 255.
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
clear g
input = 1:255;

r = 0;
c = 46;

%   Finner max og min i orignal bildet
f1 = min(f(:));
f2 = max(f(:));

%   Bestemmer max og min i det nye intervallet
g1 = 0;
g2 = 255;

%Logaritmen
log_t = log(input+r);
%Normaliserer, legg merke til at dette er det samme som når vi endret
%intervallet
log_t = g1 + ((g2-g1)/(max(log_t(:))-min(log_t(:))))*(log_t-min(log_t(:)));

%n'te rot
root_t = sqrt(input);
root_t = g1 + ((g2-g1)/(max(root_t(:))-min(root_t(:))))*(root_t-min(root_t(:)));

%n'te potens
potens_t = input.^(2);
potens_t = g1 + ((g2-g1)/(max(potens_t(:))-min(potens_t(:))))*(potens_t-min(potens_t(:)));

%inverse log
invlog_t = exp(input/50);
invlog_t = g1 + ((g2-g1)/(max(invlog_t(:))-min(invlog_t(:))))*(invlog_t-min(invlog_t(:)));

%Viser frem de forskjellige transformasjonene
figure(6)
clf
plot(input)
hold all
plot(log_t)
plot(root_t)
plot(potens_t)
plot(invlog_t);
xlabel('Input gratone')
ylabel('Output gratone')
legend('Identitet','Logaritme','n-te rot','n-te potens','invers logaritme','Location','NW');

%   Bruker så log transformasjonen på bildet
%   Her bruker vi en såkalt LUT, (look-up-table) implementasjon hvor vi
%   slår opp i transformasjonen vår og ser (look up) hvilke intensitet vi
%   skal ha ut for den aktuelle inngangs intensitete f(x,y). Dette gjøres
%   for alle x og y.
for x = 1:n
    for y = 1:m
        g(x,y) = log_t(f(x,y));
    end
end

%Viser frem resultatet
figure(7)
subplot(221)
imshow(f,[0 255]);
title('Original');
subplot(222)
imshow(g,[0 255]);
title('Etter log transform');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(uint8(g)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histogram etter transformasjon')

%%  Ultralydbilde
%   I ultralydavbildning er det vanlig å se på bildet i DB (Desibel). 
%   Dette er i prinsippet kun en log transformasjon, og tar derfor med et
%   bilde jeg har simulert i Field II med som et eksempel.
load UltrasoundImageSim

%   Bestemmer max og min i det nye intervallet
g1 = 0;
g2 = 255;

%   Her finner vi DB av bildet (legger også til minimum for å kun ha
%   positive verdier)
img = uint8(20*log10(us_image)+abs(min(20*log10(us_image(:)))));

%   Vi skal sepå bildet i to versjoner, ett hvor bildet er i DB (slik det
%   vanligvis er i ultralydavbildning) og ett bildet hvor vi har endret
%   gråtonenivået.
%   img = bildet i DB
%   g = DB bildet transformert til å dekke hele gråtonenivået

%   Finner max og min i original bildet
f1 = min(img(:));
f2 = max(img(:));

%   Normaliserer bildet fra DB til 0 - 255 i gråtoner.
g = g1 + ((g2-g1)/(f2-f1))*(img-f1);

figure(8)
subplot(231)
imshow(us_image,[0 255]);
title('Simulert ultralyd bilde');
colorbar
subplot(232)
imshow(g,[0 255])
title('Logtransform (20*log10) / Decibel transformasjon');
colorbar
subplot(233)
imshow(img,[]);
title('Logtransform (20*log10) / Decibel transformasjon');
colorbar
subplot(234)
bar(myHist(uint8(us_image)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(235)
bar(myHist(uint8(g)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histogram etter DB transformasjon')
clear g



%%  "Power-law" (gamma)-transformasjoner
%   En mer generell måte og lage ikke-lineære tranformasjoner på er såkalte
%   gamma-transformasjoner. Disse kan skrives: s = c*i^gamma hvor s er
%   ut-intensiteten, i er input, c er en skalar og gamma typisk variabelen
%   vi varierer. Dette gir en transformasjon det er lett å forholde seg
%   til, siden vi kun trenger å endre en variabel for å få forskjellige
%   resultat.

%   Vi begynner ved å velge noen verdier for gamma
gamma = [0.01 0.5 1 2 3 4];
c = 1;  %Velger å sette skalaren c til verdien 1.

%   Regner så ut og plotter de forskjellige transformasjonene
figure(9)
clf
hold all
for i = 1:length(gamma)
    t{i} = (c*double(input).^gamma(i));
    t{i} = g1 + ((g2-g1)/(max(t{i}(:))-min(t{i}(:))))*(t{i}-min(t{i}(:)));
    plot(t{i})
    legendTxt{i} = ['gamma=',num2str(gamma(i))];
end
ylabel('Input gratone');
xlabel('Output gratone');
title('"Power-law" gamma-transformasjoner');
legend(legendTxt,'Location','NW');

%   Bruker her også LUT. Se tidligere beskrivelse
%   Her kan vi velge hvilke av transformasjonene ved å velge forskjellig
%   t{i} 
for x = 1:n
    for y = 1:m
        g(x,y) = t{5}(f(x,y));
    end
end

%   Plotter resultatet
figure(10)
clf
subplot(221)
imshow(f,[0 255]);
title('Orginal');
subplot(222)
imshow(g,[0 255]);
title('Etter gamma-transformasjon');
subplot(223)
bar(myHist(f));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(uint8(g)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')

%%  Stykkvis lineær mapping
%   Stykkvis lineær mapping er en transformasjon som er "stykkvis" lineær. 
%   Her er det implementert en to-delt transformasjon hvor vi velger 
%   stigningstallet for første del opptil valgt gråtone, for a. Resten av
%   transformasjonen stiger så jevnt til øverste gråtone 255. Dette kan
%   selvfølgelig gjøres på mange måter.

%   Laster inn et bilde. Dette er ett av båndene i et mangekanals radar
%   bilde
img = imread('tm4.png');

gratone_skille = 15;    %Gråtone hvor vi "bytter" lineær mapping
a = 5;                  %Gain på første del

%   Lager transformasjonen
stykkvisT = [input(1:gratone_skille)*a linspace(input(gratone_skille)*a,255,255-gratone_skille)];

%   Plotter transformasjonen
figure(11)
clf
plot(input)
hold on
plot(stykkvisT,'r')
ylabel('Input gratone');
xlabel('Output gratone');
title('Stykkvis linear transformasjon');
legend('Indentiy','Stykkvis linear','Location','NW');

%   Implementert med LUT
for x = 1:size(img,1)
    for y = 1:size(img,2)
        img_g(x,y)=stykkvisT(img(x,y));
    end
end

%   Plot resultatet
figure(12)
subplot(221)
imshow(img,[0 255]);
title('Original');
subplot(222)
imshow(img_g,[0 255])
title('Etter stykkvis linear mapping');
subplot(223)
bar(myHist(uint8(img)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Orginal histogram')
subplot(224)
bar(myHist(uint8(img_g)));
xlabel('Gratone');
ylabel('Antall piksler');
title('Histogram etter stykkvis linear mapping')

%%  Bit-plan oppdeling
%   Her tar vi for oss hvert bit i pikslenes representasjon. Vi bruker her
%   8 bit. Da har vi feks at 00101010 som er verdien 42. Altså vil verdien 2
%   kun vises frem om vi ser på plan 1. Vi nummererer lagene som følgende:
%
%   0   =   plan 0
%   1   =   plan 1
%   0   =   plan 2
%   1   =   plan 3
%   0   =   plan 4
%   1   =   plan 5
%   0   =   plan 6
%   0   =   plan 7

layer = 1; %Velg hvilke bit (0-7)
for i = 1:n
    for j = 1:m
        %   Ved en and operasjon kan vi teste om bittet er satt.
        %   Eks om vi har tallet 42 : 00101010
        %   vil vi om vi "and'er" dette med 2: 00000010 (2^1 altså plan 1)
        %   få : 00101010 (and) 00000010 = 1 siden dette bittet er satt i
        %   plan1.
        %   For plan 2 får vi: 
        %   00101010 (and) 00000100 = 0 siden bittet her ikke er satt.
        out(i,j) = bitand(f(i,j),uint8(2^layer));
    end
end

%   Plot resultatet
figure(13)
subplot(121)
imshow(f,[0 2^8]);
title('Original');
subplot(122)
imshow(out,[]);
title(['Bit-plan: ',num2str(layer)]);
