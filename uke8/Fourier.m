%   Denne uken handlet om todimensjonal Fourier-transformasjon for bilder og
%   vi skal se på noen eksempler ved hjelp av MATLAB.

clear all
close all

%   Først skal vi se på eksempelet vi gikk gjennom på forelesning med et
%   4x4 bilde, hvor vi fant Fourier transformasjonen ved å bruke
%   definisjonen på transformasjonen (slide 12) og lage alle cosinus og
%   sinus bildene (alle matrisene i basisen). Realdelen til en frekvens (u,v) i
%   Fouriertransformasjonen er da prikkproduktet mellom bildet og cosinus
%   bildet med denne frekvensen (u,v) mens imaginærdelen er gitt av
%   prikkproduktet med sinus bildet for samme frekvens (u,v).
%   Dette gjøres for alle frekvenser u og v.
%   Se foilene for bedre forklaring :)

N = 4;  %Velger oss størrelse 4x4
M = 4;
origImg = sinusImage(N,M,0,1,1);   %Evt bruke sinus bilde som bilde
%origImg = [1 3 2 1; 5 4 5 3; 4 1 1 2; 2 3 2 6;] %Eksempelet fra forelesning
u = [0:N-1];                        %Horisontale frekvenser
v = [0:M-1];                        %Vertikale frekvenser


subIdx = 1;
for i = 1:length(u)                 %For alle horisontale frekvenser
    for j = 1:length(v)             %For alle vertikale frekvenser
        %Lager cosinusbildet for denne frekvensen
        cosImg{i,j} = cosinusImage(N,M,u(i),v(j),1);    
        %Lager sinusbildet for enne frekvensen
        sinImg{i,j} = sinusImage(N,M,u(i),v(j),1);
        
        %Plotter alle cosinus bildene
        figure(1)
        colormap gray
        subplot(N,M,subIdx)
        imagesc(cosImg{i,j});
        title(['u=',num2str(u(i)),', v=',num2str(v(j))]);
        
        %Plotter alle sinus bildene
        figure(2)
        subplot(N,M,subIdx)
        imagesc(sinImg{i,j});
        colormap gray
        title(['u=',num2str(u(i)),', v=',num2str(v(j))]);
        
        %Her gjøres selve "transformasjonen" for denne frekvensen ved å ta
        %prikkproduktet mellom bildet og cosinusbildet for realdel, og
        %sinusbildet for imaginærdel
        F(i,j) = sum(sum(cosImg{i,j}.*origImg)) - 1i*sum(sum(sinImg{i,j}.*origImg));
        
        %Index kun brukt til plotting
        subIdx = subIdx + 1;
    end
end

%Lager aksen, våre frekvenser går fra 0 til N-1
axis = linspace(0,N-1,N)
%Shiftet akse
axis_shifted = linspace(-(N)/2,(N-2)/2,N)
%Hvorfor blir vår akse slik? Tegn opp og se symmetrien, se evt side 238 i
%DIP


figure(3)
subplot(231)
imagesc(origImg)
title('Bildet');
subplot(232)
imagesc(axis,axis,abs((F)))
title('Fourier spectrum');
subplot(233)
imagesc(axis,axis,abs(fft2(origImg)))
title('MATLABs Fourier spectrum');
subplot(234)
imagesc(axis_shifted,axis_shifted,abs(fftshift((fft2(origImg)))))
title('Shifted Fourier spectrum')
subplot(235)
imagesc(real((ifft2(F))))
title('Invers Fourier');
%Ser til slutt at bildet 

%% Fouriertransformen til et cosinusbilde og et sinus bilde med samme frekvens
%   Her skal vi se at magnituden til Fourierspekteret er lik for både sinus
%   og cosinusbildet, men at verdien ligger som den reelle koefisienten for
%   cosinusbildet og som den imaginære coeffisienten for sinusbildet. Dette
%   ser vi ved å se på den relle og imaginære delen hver for seg, samt å se
%   på fasen
clear all
close all
M = 128;                    %Størrelsen
N = 128;
u = 10;                     %Horisontal frekvens
v = 5;                      %Vertikal frekvens
x = linspace(0,M-1,M);
y = linspace(0,N-1,N);
A = 50;                     %Amplituden

sinImg = zeros(N,M);
cosImg = zeros(N,M);
for i = 1:length(x)
    for j = 1:length(y)
        sinImg(i,j) = 50*sin(2*pi*(u.*x(i)/M+v.*y(j)/N));
        cosImg(i,j) = 50*cos(2*pi*(u.*x(i)/M+v.*y(j)/N));
    end
end

%Shiftet akse
axis_shifted = linspace(-(N)/2,(N-2)/2,N)

%Finner Fourier transformasjonen til sinusbildet og cosinusbildet
F_cos = fft2(cosImg);
F_sin = fft2(sinImg);

F_cos_real = real(fft2(cosImg));
F_sin_real = real(fft2(sinImg));

F_cos_imag = imag(fft2(cosImg));
F_sin_imag = imag(fft2(sinImg));

%Sørger for å fjerne alle numeriske unøyaktigheter
F_cos(abs(F_cos)<10^-9) = 0;
F_sin(abs(F_sin)<10^-9) = 0;
F_cos_real(abs(F_cos_real)<10^-9) = 0;
F_sin_real(abs(F_sin_real)<10^-9) = 0;
F_cos_imag(abs(F_cos_imag)<10^-9) = 0;
F_sin_imag(abs(F_sin_imag)<10^-9) = 0;
%Finner også den normaliserte magnituden, skal bruke denne for å finne
%fasen.

h = figure(4)
set(h, 'Position', [100, 100, 1049, 895]);
subplot(521)
imshow(cosImg,[])
title(['Cosinusbilde; u=',num2str(u),', v = ',num2str(v)]);
subplot(522)
imshow(sinImg,[])
title(['Sinusbilde; u=',num2str(u),', v = ',num2str(v)]);
subplot(523)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_cos))+1));
colorbar
xlim([-20 20]);
ylim([-20 20]);
title('Fourier spekteret til cosinusbildet')
subplot(524)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_sin))+1));
colorbar
xlim([-20 20]);
ylim([-20 20]);
title('Fourier spekteret til sinusbildet')

%Kommentar:
%   Vi ser at Fourier spekterne til cosinus bildet og sinus bildet er like.
%   Fourier spekteret er som vi vet absoluttverdien (magnituden) til
%   Fouriertransformasjonen. Og spekterne er som vi forventer likt for både
%   cosinus og sinus bildet.

subplot(525)
imagesc(axis_shifted,axis_shifted,fftshift(angle(F_cos)));
xlim([-20 20]);
ylim([-20 20]);
colorbar
title('Fase spekteret til cosinusbildet')
caxis([-pi/2 pi/2])

%Kommentar: 
%   Fasespekteret til cosinusbildet er, om vi ser bort i fra numeriske
%   unøyaktigheter, lik 0. Legg merke til at størrelsen på verdiene er i
%   størrelsesorden 10^-15, altså ser vi kun numeriske unøyaktigheter.

subplot(526)
imagesc(axis_shifted,axis_shifted,fftshift(angle(F_sin)));
xlim([-20 20]);
ylim([-20 20]);
title('Fase spekteret til sinusbildet')

%Kommentar: 
%   Fasespekteret til sinusbildet har derimot verdiene pi/2 og -pi/2
colorbar
subplot(527)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_cos_real))+1));
xlim([-20 20]);
ylim([-20 20]);
title('Realdelen av Fourier til cos')

colorbar
subplot(528)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_sin_real))+1));
xlim([-20 20]);
ylim([-20 20]);
title('Realdelen av Fourier til sin')
colorbar

colorbar
subplot(529)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_cos_imag))+1));
xlim([-20 20]);
ylim([-20 20]);
title('Imagniærdelen av Fourier til cos')

colorbar
subplot(5,2,10)
imagesc(axis_shifted,axis_shifted,log(fftshift(abs(F_sin_imag))+1));
xlim([-20 20]);
ylim([-20 20]);
title('Imaginærdelen av Fourier til sin')
colorbar

%% For å illustrere poenget bedre ser vi på den endimensjonale fouriertransformarsjonen
%   Her er poenget å bytte mellom et sinus signal og et cosinus signal, og 
%   se at verdien da flytter seg fra den imaginære komponenten til den
%   reelle komponenten.
%   Det er også interessant å 

close all
clear all

f = 5;                      %Center frequency   5 MHz
fs = 512;                   %Sampling frequency 512 MHz
x = 0:1/fs:((5/f)-1/fs)     %Med disse samplene treffer vi perioden akkurat
%x = 0:1/fs:((5/f))     %Med disse samplene bommer vi på perioden. Se hvor mye mer "støy" vi får!

s= cos(2*pi*f.*x); t = 'cos';
s= sin(2*pi*f.*x); t = 'sin';

axis = linspace(-(fs)/2,(fs-2)/2,length(s));
h = figure(5)
set(h, 'Position', [100, 100, 1049, 895]);
subplot(4,2,[1 2])
plot(s,'r','LineWidth',2)
title(['The signal: ',t]);
subplot(4,2,[3 4])
F = abs(fftshift(fft(s)));
F_m = F/max(F);
plot(axis,db(F_m))
title('Magnitude');
xlim([-20 20])
ylim([-400 0])
subplot(425)
F_r = real(fftshift(fft(s)));
plot(axis,db(F_r/max(F)))
title('Real');
xlim([-20 20])
ylim([-400 0])
subplot(426)
F_i = imag(fftshift(fft(s)));
plot(axis,db(F_i/max(F)))
title('Imag')
xlim([-20 20])
ylim([-400 0])
subplot(4,2,[7 8])
F = fft(s);
F(abs(F)<eps*500) = 0;
plot(axis,(fftshift((angle(F))))/pi);
title('Phase');
xlim([-20 20])
ylim([-1 1])
ylabel('Phase in \pi');
title('Phase');
%%

n = 2048;
m = n;

img = zeros(n,m);
img(n/2-n/4:n/2+n/4,m/2-m/16:m/2+m/16) = 1;
axis = linspace(-n/2,n/2,n);
figure(6)
subplot(231)
imshow(img,[]);
subplot(234)
imagesc(axis,axis,log(abs(fftshift(fft2(img)))+1))
colormap jet
ylim([-50 50])
xlim([-50 50])
%axis([-100 100 -100 100])

B = imrotate(img,45);
subplot(232)
imshow(B,[]);
subplot(235)
imagesc(axis,axis,log(abs(fftshift(fft2(B)))+1))
colormap jet
ylim([-50 50])
xlim([-50 50])

B = imrotate(img,90);
subplot(233)
imshow(B,[]);
subplot(236)
imagesc(axis,axis,log(abs(fftshift(fft2(B)))+1))
colormap jet
ylim([-50 50])
xlim([-50 50])
%% Invers Fourier
img = imread('car.png');
img = imread('cameraman.tif');
img = imread('garden.png');


h = figure(198)
set(h, 'Position', [100, 100, 1049, 895]);
clf
subplot(231)
imshow(img)
title('Image');
subplot(232)
F = fftshift(fft2(img));
imagesc(db(abs(F)));
title('Spectrum');
subplot(233)
imagesc(angle(F));
title('Phase');
subplot(234)
imshow((ifft2(ifftshift(abs(F)))),[]);
title('Ifft(Magnitude)');
subplot(235)
F = fft2(img);
phase = F./abs(F);
imshow(abs(ifft2(phase)),[]);
title('Ifft(Phase)');
subplot(236)
imshow(abs(ifft2(ifftshift(F))),[]);
title('Ifft(Phase and magnitude)');
