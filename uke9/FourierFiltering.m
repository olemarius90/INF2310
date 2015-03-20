% Denne uken har vi sett mer på Fourier transformasjonen og særlig sett på
% design av Filtere i Fourierdomenet, samt hvordan vi kan analysere
% filteret vi allerede har i bildedomenet ved å ta Fouriertransformasjonen
% av filterne og se på frekvens responsen.

% Sentralt i denne forelesningen var var konvolusjonsteoremet som sier at
% sirkelkonvolusjon i bildedomenet er punktvis multiplikasjon i Fourierd-omenet.
% Og omvendt at punktvis multiplikasjon i bildedomenet er
% sirkelkonvolusjon i Fourier-domenet.

% Ideelt lavpass filter
% Vi begynner ved å se på det såkaldte ideelle lavpassfilteret. Her er det
% viktig å huske at selv om filteret heter ideelt, er det ikke ideelt ;)
% Det har fått navnet siden overgangen mellom passbåndet (hvor verdiene er
% 1) til stopbåndet (hvor verdiene er 0) er maksimal bratt. Vi skal se at
% dette introduserer ringing i bildedomenet.
clear all
close all

% Vi ser først på bildet av mister cameraman.
f = double(imread('cameraman.tif'));
% Finner størrelsen
[M,N] = size(f);

% Definerer tomme variable for H og D
H = zeros(M,N);
D = zeros(M,N);
% Velger cutoff D0 lik 0.2
D0 = 0.2;

%Regner ut det ideele lavpassfilteret
for i = 1:M
    for j = 1:N
        %Verdiene for D er som definert i slide 17 i forelesningen.
        %Se forklaring på Foil 20 om hvorfor vi har floor(M/2+1) for
        %frekvensene. Frekvensene blir altså u = i-[m/2+1] osv.
        D(i,j) = sqrt(((i-floor(M/2 + 1))/(M/2))^2+(((j-floor(N/2+1))/(N/2))^2));
        %Vi setter verdien til filteret til 1 når vi er innenfor D0.
        if D(i,j) <= D0;
            H(i,j) = 1;
        end
    end
end


%Vi plotter D og det resulterende filteret
figure(1)
subplot(211)
imshow(D)
colorbar
title('Verdiene for D');
subplot(212)
imshow(H,[]);
title('Det ideelle lavpassfilteret i Frekvensdomenet');

%   Vi ønsker å gjøre filtreringen i frekvensdomenet, finner derfor
%   Fouriertransformen til bildet.
F = fftshift(fft2(f));
%   Så punktmultipliserer vi med filteret vårt og tar invers Fourier for å
%   få det filtrerte bildet
g = real(ifft2(ifftshift(F.*H)));

figure(4)
subplot(221)
imshow(f,[]);
title('Orginal bildet');
subplot(222)
imagesc(abs(log(F+1)));
title('Fourier-spekteret til bildet');
subplot(223)
imshow(g,[]);
title('Lavpass filtrert bilde');
subplot(224)
imagesc(log(abs(fftshift(fft2(g)))+1));
title('Fourier-spekteret til det filtrerte bildet');

%Legg merke til ringingen vi har fått i bildet!!
%%  Butterworth
%   Nå vil vi se på det andre typen filteret vi lærte, nemlig Butterworth
%   filteret. Se definisjonen på foil 21. Dette filteret har en mykere
%   overgang mellom passbåndet og stoppbåndet, brattheten på overgangen
%   styres av filterordenen n. Hvis n er liten blir det lite ringing, om n
%   er 1 blir det faktisk ingen ringing. Prøv med forskjellige verdier for
%   n og se resultatet
%   D0 angir nå avstanden fra DC til punktet der H har sunket med 0.5
n = 3;
D0 = 0.2;
%   Lager filteret
H_butter = 1./(1+(D./D0).^(2*n));
%   Gjør filtreringen
g_butter = real(ifft2(ifftshift(F.*H_butter)));

figure(5)
subplot(131)
imagesc(H)
title('Det ideelle lavpassfilteret');
subplot(132)
imagesc(H_butter)
title('Butterworth lavpassfilteret');
colormap gray

%Plotter resultatet. Legg merke til at vi har fått mindre ringing i bildet
figure(6)
subplot(221)
imshow(f,[]);
subplot(222)
imagesc(abs(log(F+1)));
subplot(223)
imshow(g_butter,[]);
subplot(224)
imagesc(log(abs(fftshift(fft2(g_butter)))+1));

%% Gaussisk lavpassfilter
%   Det tredje og siste filteret vi så på var et Gaussisk lavpassfilter.
%   Fordelen med et Gaussfilter er at Fouriertransformasjonen til en Gauss
%   er en Gauss - altså får vi ikke ringen i bildet. Men for et Gauss
%   filter kan vi ikke kontrollerer hvor rask overgang vi har og mister
%   derfor litt kontroll over frekvensresponsen.
%   D0 indikerer her hvor filterverdien er 0.6.
D0 = 0.2;
%   Lager filteret som definert i forelesning foil 24
H_gauss = exp((-D.^2)/(2*D0^2));
g_gauss = real(ifft2(ifftshift(F.*H_gauss)));


figure(5)
subplot(133)
imagesc(H_gauss);
title('Gaussisk lavpassfilter');

%Plotter resultatet. Nok et lavpassfiltrert bilde uten ringing
figure(7)
subplot(221)
imshow(f,[]);
subplot(222)
imagesc(abs(log(F+1)));
subplot(223)
imshow(g_gauss,[]);
subplot(224)
imagesc(log(abs(fftshift(fft2(g_gauss)))+1));


%% Nå vil vi se på de romlige responsene til filterne, altså hvordan de ser
%  ut i bildedomenet. Når vi går fra Fourierdomenet til bildedomenet vil
%  filterne typisk være store, tenk på hvordan vi nullutivder et romlig
%  filter for å bruke dette i Frekvensdomenet. Allikevel er dette stort
%  sett bare verdien 0, og vi kan kutte ut dette i konvolusjonen. Vi ser
%  derfor bare på den sentrale delen av filteret. Vi velger hvor mye av
%  filteret vi vil se på ved hjelp av k_size, og vi ser på 2xk_size.

k_size = 50;
M_cen = floor(M/2+1);
N_cen = floor(N/2+1);
%Romlig representasjon av ideelt lavpass
filter = fftshift(real(ifft2(ifftshift(H))));
kernel_ideal = filter(M_cen-k_size:M_cen+k_size,N_cen-k_size:N_cen+k_size);
%Romplig representasjon av Butterworth
filter = fftshift(real(ifft2(ifftshift(H_butter))));
kernel_butter = filter(M_cen-k_size:M_cen+k_size,N_cen-k_size:N_cen+k_size);
%Romlig representasjon av Gaussisk lavpassfilter
filter = fftshift(real(ifft2(ifftshift(H_gauss))));
kernel_gauss = filter(M_cen-k_size:M_cen+k_size,N_cen-k_size:N_cen+k_size);


figure(9)
subplot(311)
imagesc(kernel_ideal)
title('Ideelt lavpass i bildedomenet');
subplot(312)
imagesc(kernel_butter)
title('Butterworth lavpass i bildedomenet');
subplot(313)
imagesc(kernel_gauss)
title('Gaussisk lavpass i bildedomenet');

%Nedefor plotter vi kun den sentrale dele av aksen i Fourier og i tid
figure(100)
subplot(321)
plot(H(round(end/2),:));
title('Ideelt i Fourier-domenet');
subplot(322)
plot(kernel_ideal(round(end/2),:));
title('Ideelt i bildedomenet');
subplot(323)
plot(H_butter(round(end/2),:))
title('Butterworth i Fourier-domenet');
subplot(324)
plot(kernel_butter(round(end/2),:))
title('Butterworth i bildedomenet');
subplot(325)
plot(H_gauss(round(end/2),:))
title('Gauss i Fourier-domenet');
subplot(326)
plot(kernel_gauss(round(end/2),:))
title('Gauss i bildedomenet');

%Prøv å endre n for Butterworth filteret og hvordan dette endrer responsen
%i bildedomenet. Skjønner du hvorfor ringingen blir borte??
%%  Så skal vi se at å konvolusjonsteoremet faktisk fungerer og
%   konvolusjon i bildedomenet er det samme som punktmultiplikasjon i
%   Frekvensdomenet.

%   Bruker responsen til det ideelle filteret i bildedomenet som
%   konvolusjonsfilter og passer på å bruke 'symmetric' i imfilter. 
%   'symmetric' sørger for at det er sirkulærindeksering som brukes når vi
%   utvider kanten.
g_bildedomenet = imfilter(f,kernel_ideal,'conv','symmetric');

figure(10)
subplot(131)
imshow(f,[])
title('Det orginale bildet')
subplot(132)
imshow(g_bildedomenet,[])
title('Bildet filtrert i bildedomenet');
subplot(133)
imshow(g,[])
title('Bildet filtrert i Fourier-domenet');

%Legg merke til at resultatet er likt.


%% Notch filter example
%   Her skal vi se på et notch filter, samme eksempel som i forelesningen,
%   som fjernet et uønsket mønster (moire pattern) som ligger "over" bildet
%   vi er interessert i
clear all
close all
f = double(imread('car2.png'));
[M,N] = size(f);
F = fft2(f);

u = linspace(-floor(M/2),floor(M/2-1),M)
v = linspace(-floor(N/2),floor(N/2-1),N)

figure(11)
subplot(221);
imshow(f,[]);axis on;
title('Bildet');
subplot(222)
imshow(log(abs(fftshift(F))+1),[]);axis on;
title('Fourier spekteret');

% Lokaliser frekvensene til st�yen (gjort manuelt v.b.a. Fourier-spekteret)
% (frekvensene er her angitt ved array-indeksene etter FFTSHIFT).
% De mest markante frekvensene av st�yen er:
% (166,58) og (86,56) (og deres symmetriske "kopier")
% og de mindre tydelige, men likevel markante frekvensene av st�yen er:
% (41,112) og (45,56) (og deres symmetriske "kopier")
u = [166  86  41  45];
v = [ 58  56 112  56];
% Spesifiser bare de mest markante frekvensparene (er to par):
%u = [166  86];
%v = [ 58  56];
u_sym = M + 2 - rem(M,2) - u;
v_sym = N + 2 - rem(N,2) - v;
u = [u u_sym];
v = [v v_sym];


% Butterworth notch filter.
D0 = 0.1;
n = 4;                                    %Antall hull
D = zeros(M,N,length(u));
H = ones (M,N);
for x = 1:M
    for y = 1:N
        for k = 1:length(u)
            D(x,y,k) = sqrt(((x-u(k))/(M/2))^2 + ((y-v(k))/(N/2))^2);
            H(x,y) = H(x,y) / (1 + (D0/D(x,y,k))^(2*n));
        end
    end
end

% Filtrering.
G = F.*ifftshift(H);

figure(11)
subplot(223); imshow(log( abs(fftshift(G)) + 1), []); axis on;
title('Det filtrerte fourierspekteret');
subplot(224); imshow(real(ifft2(G)),[]); axis on;
title('Det filtrerte bildet');

%% Videre skal vi se at vi kan analysere noen "kjente" convolusjonsfiltere
%   ved å ta fourier transformen og se på frekvensresponsen til filterne

hx = [0 1 0; 0 0 0; 0 -1 0];
hy = [0 0 0; 1 0 -1; 0 0 0];
PrewittHy = [1 0 -1; 1 0 -1; 1 0 -1];

HX = fftshift(abs(fft2(hx,512,512)));
HY = fftshift(abs(fft2(hy,512,512)));
PHY = fftshift(abs(fft2(PrewittHy,512,512)));

figure(100)
subplot(211)
imagesc(log(HX+1));
subplot(212)
imagesc(log(HY+1));

figure(101)
imagesc(log(PHY + 1));
%%  Filtering with zero-padding.
%   Nullutvidelse for å ha en annen måte å takle kantproblematikken på. Se
%   se foilen for nærmere beskrivelse.
fp = [f zeros(M,N) ; zeros(M,N) zeros(M,N)];
Fp = fft2(fp);
Mp = 2*M;
Np = 2*N;
up = u.*2-1;
vp = v.*2-1;
figure(15);subplot(221); imshow(fp,[]); axis on;
subplot(222); imshow(log( abs(fftshift(Fp)) + 1), [0 max(log( abs(Fp(:)) + 1))]); axis on;
%
% Butterworth notch filter.
D0 = 0.1;
n = 4;
Dp = zeros(Mp,Np,length(up));
Hp = ones (Mp,Np);
for x = 1:Mp
    for y = 1:Np
        for k = 1:length(up)
            Dp(x,y,k) = sqrt(((x-up(k))/(Mp/2))^2 + ((y-vp(k))/(Np/2))^2);
            Hp(x,y) = Hp(x,y) / (1 + (D0/Dp(x,y,k))^(2*n));
        end
    end
end
% Filtering.
Gp = Fp.*ifftshift(Hp);
gp = real(ifft2(Gp));
subplot(223); imshow(log( abs(fftshift(Gp)) + 1), [0 max(log( abs(Gp(:)) + 1))]); axis on;
subplot(224); imshow(gp,[]); axis on;


%
figure(14);clf; ax1 = subplot(232); imshow(gp(1:M,1:N),[]); axis on; title('Zero padded');
ax2 = subplot(231); imshow(real(ifft2(G)),[]); axis on; title('Circular index')
%Ser på forskjellen
subplot(233); imshow(real(ifft2(G))-gp(1:M,1:N),[]);
linkaxes([ax1 ax2],'xy');

%%
%   Viser igjen at å convolvere i bildedomenet er det samme som
%   punktmultiplikasjon i Fourier-domenet
%
%   Finner inverstransformasjonen til Notch filteret
filter = fftshift(real(ifft2(ifftshift(H))));
figure(1)
subplot(121)
imagesc(log(H+1));
subplot(122)
imagesc(log(filter+1));

%   Convolverer filterne og bildene i bildedomenet. Her velger vi å bruke
%   sirkulær indeksering som en måte å løse kantproblement på, for det
%   andre bildet brukes det nullutvidelse. Legg merke til at resultatet
%   blir likt som når vi gjør det i Frekvensdomenet!!!
res = imfilter(f,filter,'conv','circular');
res_padded = imfilter(f,filter,'conv');

%   Vi ser at resultatene blir like
figure(14)
subplot(234)
imshow(res,[])
subplot(235)
imshow(res_padded,[])
subplot(236)
%   Ser på forskjellen
imshow(res-res_padded,[])
%% Vindusfunksjon
%   Til slutt skal vi se at å legge en vindusfunksjon på bildet før vi tar
%   Fourier transformasjonen reduserer diskontinuiteten og gir mindre akse
%   bidrag.
%
%   Det pedagogiske poenget her er og nok en gang se dualiteten til
%   konvolusjonsteoremet. Når vi punktmultipliserer med et vindu i
%   bildedomenet vil det si at vi konvolverer med noe som ligner på et
%   middelverdifilter i frekvens.
clear all
close all
f = double(imread('car2.png'));
[M,N] = size(f);
%   Vi lager vinduet
alpha = 0.5; w = (tukeywin(N,alpha)*tukeywin(M,alpha)')';
%   Punktmultipliserer bildet med vinduet
fw = f .* w;
F_window = fftshift(abs(fft2(fw)));
F = fftshift(abs(fft2(f)));

%   Fouriertransformasjonen til vinduet
W = fftshift(abs(fft2(w)));

%   Viser frem bildene og resultatene
figure(15);
subplot(231)
imshow(f,[]);
title('Bildet');
subplot(232)
imshow(w);
title('Vinduet');
subplot(233)
imshow(fw,[]);
title('Bildet punktmultiplisert med vinduet');
subplot(234)
imshow(log( abs(F_window) + 1), [0 max(log( abs(F_window(:)) + 1))]); axis on;
title('Fourier spekteret til bildet');
subplot(235)
imshow(log(W+1))
title('Fourier spekteret til vinduet');
subplot(236); imshow(log( abs(F) + 1), [0 max(log( abs(F(:)) + 1))]); axis on;
title('Fourier spekteret til bildet med ramme');
