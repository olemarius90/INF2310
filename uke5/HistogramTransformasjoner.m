% HistogramTransformasjoner
%
% Eksempelkode skrevet for å illustrere temaene gått gjennom i forelesning
% 5 i faget INF2310 - Digital bildebehandling ved Ifi UiO.
%
% Skrevet av Ole Marius Hoel Rindal
% Vennligst gi beskjed til olemarius@olemarius.net om du finner feil eller
% mangler i koden.
clear all
close all

%Legger til path'en til histogramfunksjonen vi skrev forrige gang
if isunix
    addpath ../uke4
else
    addpath ..\uke4
end

% I dag importerer vi et bilde av en bil
img = imread('car.png');
[n,m] = size(img);  %Finner størrelsen

%   Av "pedagogiske" årsaker lager vi oss en noe dårligere versjon av
%   inngangsbildet...
%   Finner max og min i orignal bildet
f1 = min(img(:));
f2 = max(img(:));

%   Bestemmer max og min i det nye intervallet
g1 = 50;
g2 = 110;
%   Dette er i grunnen kun en normalisering.
f = double(g1) + ((double(g2)-double(g1))/(double(f2)-double(f1)))*(double(img)-double(f1));
f = uint8(f);

%   Resultatet bli at f nå er et bilde med færre gråtoner og dårligere
%   kontrast enn img, så ser vi lettere forskjell.
figure(1)
subplot(221)
imshow(img,[0 255]);
title('Orginal');
subplot(222)
imshow(f,[0 255])
title('Gråtone begrenset');
subplot(223)
bar(myHist(img))
subplot(224)
bar(myHist(f))


%%  Histogramutjevning
%   Vi skal finne en transformasjon som gir oss en histogramutjevning på
%   inngangsbildet. Altså er vi ute etter en transformasjon som gir oss
%   histogram som er mest mulig likt et uniformt histogram.
%   
%   Vi trenger da en transformasjon som har bratt stigningstall når vi har
%   mange piksler per gråtone, og lavere stigningstall når vi har mindre
%   piksler per gråtone. Det normaliserte kummulative histogrammet har
%   nettopp disse egenskapene, og ved å bruke dette som transformasjon får
%   vi et histogram som er tilnærmet et uniform. Det kumulative
%   histogrammet til et uniformt histogram vil være en rett linje.
%
%   Et bilde med "maksimal kontrast" vil ha et uniformt histogram, vi kan
%   altså forvente at en bilde som er histogramutjevnet vil ha mer kontrast
%   enn orginalen.
%
% Bruker vår funksjon myHist til å finne histogrammet
% Se myHist.m filen for detaljer på hvordan man kan finne histogrammet.
[p,h,c,c_n] = myHist(f);


% Finner transformasjonen som gir histogramutjevning ved hjelp av det
% normaliserte kumulative histogrammet slik vi har sett i forelesningen
G = 2^8;
T = round((G-1)*c_n);

% Viser frem bildet, histogrammet, det kumulative histogrammet og
% transformasjonen
h1 = figure(2)
subplot(221)
imshow(f,[0 255])
set(gca(h1),'fontSize',14)
title('Bildet');
subplot(222)
bar(p)
title('Histogrammet');
set(gca(h1),'fontSize',14)
subplot(223)
bar(c_n)
title('Kumulativt histogram')
set(gca(h1),'fontSize',14)
subplot(224)
plot(T,'LineWidth',2);
axis tight
title('Transformasjonen');
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

%%  Vi gjør transformasjonen på inngangsbildet.
for i = 1:n
    for j = 1:m
        g(i,j) = T(f(i,j));
    end
end

%   For sammenligningens skyld lager vi også et bilde hvor vi normaliserer
%   bilde til å bruke hele gråtoneskalaen
%   Finner max og min i orignal bildet
f1 = min(f(:));
f2 = max(f(:));

%   Bestemmer max og min i det nye intervallet
g1 = 0;
g2 = 255;

%   Dette er i grunnen kun en normalisering.
g_2 = double(g1) + ((double(g2)-double(g1))/(double(f2)-double(f1)))*(double(f)-double(f1));
g_2 = uint8(g_2);

%   Viser frem resultatet
h1 = figure(3)
subplot(331)
imshow(f,[0 255]);
title('Original');
set(gca(h1),'fontSize',14)  %Disse linjene er bare for å lage større skrift i titlene.
                            %Som jeg brukte til forelesningen

subplot(332)
imshow(g,[0 255])
title('Histogramutjevning');
set(gca(h1),'fontSize',14)

subplot(333)
imshow(g_2,[0 255])
title('Strekking ved hjelp av normalisering')
set(gca(h1),'fontSize',14)

subplot(334)
bar(myHist(f))
title('Histogrammet')
set(gca(h1),'fontSize',14)

subplot(335)
bar(myHist(g))
title('Histogrammet');
set(gca(h1),'fontSize',14)

subplot(336)
bar(myHist(g_2))
title('Histogrammet');
set(gca(h1),'fontSize',14)

subplot(337)
[p,h,c,c_n] =myHist(f);
bar(c_n)
title('Kumulativt histogram');
set(gca(h1),'fontSize',14)


subplot(338)
[p,h,c,c_n] =myHist(g);
bar(c_n)
title('Kumulativt histogram')
set(gca(h1),'fontSize',14)


subplot(339)
[p,h,c,c_n] =myHist(g_2);
bar(c_n)
title('Kumulativt histogram')
set(gca(h1),'fontSize',14)


set(findall(h1,'type','text'),'fontSize',18)    %Også kun for å forstørre skriften
set(gca(h1),'fontSize',14)


%% Verifsier at det er likt som MATLAB's implementasjon...
g = histeq(f,2^8);
figure(4)
subplot(321)
imshow(f,[0 255]);
subplot(322)
imshow(g,[0 255])
subplot(323)
bar(myHist(f))
subplot(324)
bar(myHist(g))
subplot(325)
[p,h,c,c_n] =myHist(f);
bar(c_n)
subplot(326)
[p,h,c,c_n] =myHist(g);
bar(c_n)

%%  Histogramtilpasning
%      
%   Videre skal vi gjøre histogramtilpasning. Her spesifiserer vi formen på
%   det ønskede histogrammet og lager en transformasjon som endrer det
%   orginale bildet slik at vi får tilnærmet ønsket histogram. Dette gjøres
%   i flere steg. Vi kan også bruke histogrammet til et annet bilde som
%   ønsket histogram.
%
%   Steg 1 : Lag transformasjonen, T(i),som gjør histogramutjevning på innbildet
%   Steg 2 : Gitt det ønskede histogrammet pz(i) og finn det kumulative
%   histogrammet for å lage transformasjonen, Tg(i), som histogramutjevner
%   det ønskede histogrammet.
%   Steg 3 : Finn så inverstransformen til Tg(i) og inverstransformerer
%   innbildet til det ønskede histogrammet.
%
%   Altså, må det orginale bildet først histogramutjevnes, før det
%   histogramutjevnede bildet transformeres med inverstransformen til det 
%   ønskede histogrammet. Altså må bildet gjennom to transformer.
%
%   Alternativt kan disse to transformene kombineres for å gi en
%   transformasjon som direkte transformerer innbildet til det ønskede
%   histogrammet.
%
%   Se foilene for en noe mer pedagogisk fremstilling...
%

f = img;% Jukser oss til å bruke et bedre innbilde ;)

% Først lager vi oss et ønsket histogram
% Lager her en Gaus fordeling med følengede verdier
x = 0:255;  %Input
u = 255/2;  %Middelverdi
s = 40;     %Varians

% Lager Gaus fordelingen
gauss = (1/(s*sqrt(2*pi)))*exp(-(x-u).^2/(2*s^2));
gauss = gauss/max(gauss);

% Lager det "kumulative histogrammet" fra Gaus'en
gauss_cumsum = cumsum(gauss);
gauss_cumsum = gauss_cumsum/max(gauss_cumsum);

% Viser frem Gauss'en, den kumulative Gaus'en 
h1 = figure(5)
clf
plot(gauss,'LineWidth',2)
hold all
plot(gauss_cumsum,'LineWidth',2)
legend('Normalisert Gauss','Normalisert kummulativ Gauss','Location','best');
axis tight
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

% Steg 1 : Finner transformen som histogramutjevner innbildet f
%          Dette er som før...
[p,h,c,c_n] =myHist(f);
T_histeq = round((G-1)*c_n);

% Steg 2: Finner transformer som histogramutjevner Gauss'en vår
T = round((G-1)*gauss_cumsum);

% Steg 2 forts: Finn inverstransformen av transformasjonen som
% histogramutjevner Gauss'en
T_identity = 0:255;
T_invers = zeros(1,256);

% NB! Litt hacky implementasjon for å finne inverstransformen. Burde funnet
% en bedre implementasjon. But time is money...
last = 0;
for i = 1:length(T)    
    if T(i) == 0
        T_invers(1) = i;
    else
        T_invers(T(i)) = i;
        
    end
end

% Transformasjonen T hadde ikke alle steg mellom 0 og 255, så vi tar med
% resten
for i = 1:length(T_invers)
    if (T_invers(i) == 0)
        T_invers(i) = T_invers(i-1);
    end
end

%  Plotter Tg og invers av Tg samt identiteten for referanse
h1 = figure(6)
plot(T,'r','LineWidth',2)
hold on
plot(T_invers,'b','LineWidth',2)
plot(1:255,'k','LineWidth',2)
legend('T_g','T_g^{-1}','Identitet','Location','best');
title('Transformasjoner');
xlabel('Input gratoner');
ylabel('Output gratoner');
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

%%
% Alternativt kan vi finne den direkte transformasjonen. Dette er jo
% greiest, men mer intuitivt om vi tar med begge
for i = 1:length(T_histeq)
    if T_histeq(i) == 0
        T_final(i) = T_invers(1);
    else
        T_final(i) = T_invers(T_histeq(i));
    end
end

% Plotter alle transformasjonene
h1 = figure(12)
plot(T_histeq,'LineWidth',2);
hold all
plot(T_invers,'LineWidth',2);
plot(T_final,'LineWIdth',2);
title('Transformasjoner');
legend('T : forste','T_g^{-1}','T_{direkte}','Location','best');
xlabel('Input gratoner');
ylabel('Output gratoner');
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)


% Gjør første transformasjon, altså histogramutjevner innbildet
% Gjør samtidig den direkte transformasjonen for å verifisere at det blir
% likt
for i = 1:n
    for j = 1:m
        if f(i,j) == 0
            g(i,j) = T_histeq(1);
            g_direct(i,j) = T_final(1)
        else
            g(i,j) = T_histeq(f(i,j));
            g_direct(i,j) = T_final(f(i,j));
        end
    end
end

% Gjør andre transformasjon
for i = 1:n
    for j = 1:m
        if g(i,j) == 0
            g(i,j) = T_invers(1);
        else
            g(i,j) = T_invers(g(i,j));
        end
    end
end

% Viser frem resultatet. Tar med begge utgangsbildene for å verifisere at
% det ble like
h1 = figure(13)
subplot(331)
imshow(f,[0 255]);
title('Innbildet');

subplot(332)
imshow(g,[0 255]);
title('Utbilde to trans');

subplot(333)
imshow(g_direct,[0 255]);
title('Utbilde direkte trans')
axis tight

subplot(334)
[p_f,h,c,c_n_f] = myHist(f);
[p_g,h,c,c_n_g] = myHist(g);
[p_final,h,c,c_n_final] = myHist(g_direct);
bar(p_f)
axis tight
title('Histogram');

subplot(335)
bar(p_g)
title('Histogram');
axis tight
set(gca(h1),'fontSize',14)

subplot(336)
bar(p_final)
title('Histogram');
axis tight
set(gca(h1),'fontSize',14)

subplot(337)
bar(c_n_f)
title('Kumulativt histogram');
axis tight
set(gca(h1),'fontSize',14)

subplot(338)
bar(c_n_g)
title('Kumulativt histogram');
axis tight
set(gca(h1),'fontSize',14)

subplot(339)
bar(c_n_final)
title('Kumulativt histogram');
axis tight
set(gca(h1),'fontSize',14)
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

%%  For å verifisere resultatene så sammenligner vi med matlab's implementasjon
%   Legg merke til at vi bruker histeq funksjonen, ikke imhismatch som
%   krever to innbilder som input.
g = histeq(f,gauss);
figure;bar(gauss*n*m)

%   Vi ser at resultatene er like, men at matlab har fått et noe annerledes
%   resultat. De har nok en noe mer sofistikert implementasjon, og finner
%   kanskje den ene nødvendige transformasjonen direkte. Det ble nevt i
%   forelesningen av en observant student at dette muligens gir mindre
%   avrundingsfeil, noe dette resultatet kanskje viser er riktig?? :)
figure(14)
subplot(321)
imshow(f)
subplot(322)
imshow(g)
subplot(323)
[p_f,h,c,c_n_f] = myHist(f)
bar(p_f)
subplot(324)
[p_g,h,c,c_n_g] = myHist(g)
bar(p_g)
subplot(325)
bar(c_n_f)
subplot(326)
bar(c_n_g)


%%  Fargebilder %Fikset frem til hit! 
%   Laster inn ett RGB fargebilde
img = imread('hane.jpg');

r = uint8(img(:,:,1));
g = uint8(img(:,:,2));
b = uint8(img(:,:,3));

HSI=rgb2hsv(img);

[p,h,c,c_n_r] = myHist(r);
[p,h,c,c_n_g] = myHist(g);
[p,h,c,c_n_b] = myHist(b);
[p,h,c,c_n_i] = myHist(HSI(:,:,3)*255);
%%

[J T] = histeq(HSI(:,:,3));

h1 = figure(100)
subplot(221)
imshow(HSI(:,:,3));
title('HSI:I value');
subplot(222)
imshow(J)
title('HSI:I value histeq');
subplot(223)
imhist(HSI(:,:,3));
title('Histogram');
subplot(224)
imhist(J)
title('Histogram');%
%set(gca,'DefaultTextFontSize',18)


set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)


HSI_new = HSI;
HSI_new(:,:,3) = J;

RGB = hsv2rgb(HSI_new);


%%
input = 1:255;
T_r = round((G-1)*(c_n_r));
T_g = round((G-1)*(c_n_g));
T_b = round((G-1)*(c_n_b));
T_i = round((G-1)*(c_n_i));

figure(888)
plot(T_r);
hold all
plot(T_g);
plot(T_b);
plot(T_i);


[n,m,dummy] = size(img);

HSI_new = HSI;
for i = 1:n
    for j = 1:m
        if r(i,j) == 0
            r_2(i,j) = T_r(1);
        else
            r_2(i,j) = T_r(r(i,j));
        end
        
        if g(i,j) == 0
            g_2(i,j) = T_g(1);
        else
            g_2(i,j) = T_g(g(i,j));
        end
        
        if b(i,j) == 0
            b_2(i,j) = T_b(1);
        else
            b_2(i,j) = T_b(b(i,j));
        end
        
        if round(HSI(i,j,3)*255) == 0
            HSI_new(i,j,3) = T_i(1);
        else
            HSI_new(i,j,3) = T_i(round(HSI(i,j,3)*255));
        end
    end
end

h1 = figure(1000)
clf
subplot(211)
plot(myHist(r),'r','LineWidth',2);
hold on
plot(myHist(g),'g','LineWidth',2);
plot(myHist(b),'b','LineWidth',2);
title('RGB Histogram');

subplot(212)
bar(myHist(r_2),'r');
hold on
bar(myHist(g_2),'g');
bar(myHist(b_2),'b');
title('RGB Histogram');

img_after_rgb_hist(:,:,1) = r_2;
img_after_rgb_hist(:,:,2) = g_2;
img_after_rgb_hist(:,:,3) = b_2;


h1 = figure(1001)
subplot(131)
imshow(img)
title('Original');
subplot(132)
imshow(uint8(img_after_rgb_hist),[0 255])
% band(:,:,1) = histeq(img(:,:,1));
% band(:,:,2) = histeq(img(:,:,2));
% band(:,:,3) = histeq(img(:,:,3));
% imshow(uint8(band))

%imshow(uint8(img_after_rgb_hist),[]);
title('RGB histeq')
subplot(133)
imshow(RGB)
title('HSI histeq');
%subplot(144)

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)


%%
figure(995)
clf
figure(996)
subplot(211)
imshow(img,[])
subplot(212)
img_copy = img;
img_copy(:,:,3) = img(:,:,3)*2 + 100;
imshow(img_copy,[]);

%%
figure(997)
imshow(uint8(img_after_rgb_hist),[])

figure(998)
subplot(211)
imshow(RGB)
subplot(212)
bar(myHist(uint8(HSI_new(:,:,3)*255)));

figure(999)
bar(myHist(uint8(img_after_rgb_hist(:,:,1))),'r');
hold on
bar(myHist(uint8(img_after_rgb_hist(:,:,2))),'g');
bar(myHist(uint8(img_after_rgb_hist(:,:,3))),'b');

%%
g = histeq(f);
g_2 = adapthisteq(f,'NumTiles',[15 15],'ClipLimit',0.01,'Distribution','uniform');

h1 = figure(30)
subplot(331)
imshow(f,[0 255]);
title('Original');

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)  

subplot(332)
imshow(g,[0 255])
title('Histogramutjevning');


set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(333)
imshow(g_2,[0 255])
title('CLAHE')


set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(334)
bar(myHist(f))
title('Histogram');

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(335)
bar(myHist(g))
title('Histogram');

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(336)
bar(myHist(g_2))
title('Histogram');


set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(337)
[p,h,c,c_n] =myHist(f);
bar(c_n)

title('Kumulativt Histogram');


set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(338)
[p,h,c,c_n] =myHist(g);
bar(c_n)
title('Kumulativt Histogram');
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

subplot(339)
[p,h,c,c_n] =myHist(g_2);
bar(c_n)
title('Kumulativt Histogram');

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

f_zoom = f(110:end-40,115:end-35);
g_zoom = g(110:end-40,115:end-35);
g_2_zoom = g_2(110:end-40,115:end-35);

h1 = figure(31)
subplot(131)
imshow(f_zoom,[0 255]);
title('Original');

subplot(132)
imshow(g_zoom,[0 255])
title('Histogramutjevning');

subplot(133)
imshow(g_2_zoom,[0 255])
title('CLAHE')

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)


%%

% Først lager vi oss et ønsket histogram
% Lager her en Gaus fordeling med følengede verdier
x = 0:255;  %Input
u = 255/2;  %Middelverdi
s = 40;     %Varians

% Lager Gaus fordelingen
gauss = (1/(s*sqrt(2*pi)))*exp(-(x-u).^2/(2*s^2));
gauss = gauss/max(gauss);

% Lager det "kumuliative histogrammet" fra Gaus'en
gauss_cumsum = cumsum(gauss);
gauss_cumsum = gauss_cumsum/max(gauss_cumsum);

h1 = figure(5)
clf
plot(gauss,'LineWidth',2)
hold all
plot(gauss_cumsum,'LineWidth',2)
legend('Normalisert Gauss','Normalisert kummulativ Gauss','Location','best');
axis tight

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)
over = 0;s
gauss_2 = gauss;
for i = 1:length(gauss_2)
   if gauss_2(i) >  0.5
       over = over + gauss_2(i) - 0.5;
       gauss_2(i) = 0.5;
   end
end
gauss_2 = gauss_2 + over/255;
% Lager det "kumuliative histogrammet" fra Gaus'en
gauss_cumsum_2 = cumsum(gauss_2);
gauss_cumsum_2 = gauss_cumsum_2/max(gauss_cumsum_2);

h1 = figure(5)
clf
plot(gauss,'LineWidth',2)
hold all
plot(gauss_2,'LineWidth',2)
plot(gauss_cumsum,'LineWidth',2)
plot(gauss_cumsum_2,'LineWidth',2)
legend('Gauss','Klippet Gauss','Kumulativ Gauss','Kumulativ klippet Gauss','Location','best');
axis tight
ylim([0 1.5])

set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)
