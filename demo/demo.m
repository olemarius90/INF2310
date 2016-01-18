%% Demoscript
%
%   Et enkelt script som demonstrerer nettsiden for eksempelkode
%   inf2310.olemarius.net til første forelesning.
%


%% Lese inn et bilde i MATLAB
% Først leser vi in ett bilde som er "innebygget" det vil si at bildet
% ligger et sted i MATLABS "path" og vi kan derfor lese det ved å bare
% skrive navnet på filet
img = imread('onion.png');

% Bildet kan vises frem slik
figure(1);
imshow(img);

% Men det kan også vises frem slik
figure(2);
imagesc(img)

% Da legger dere kanskje merke til at bildet er litt strukket, dette kan vi
% enkelt fikse feks slik:
figure(3);
imagesc(img)
axis image

%% Alle bilder
% For å lese inn bilder som ikke er bygget inn i matlab må man sørge for at
% bildet ligger på "path'en", feks samme mappe som filen du kjører.

% For eksempel kan man laste inn et skrytebilde av foreleseren
bilde_ole_m = imread('ole_m_skryt.jpg');
figure(4);
imagesc(bilde_ole_m);
axis image

% Man kan også legge til en tittel
title('Et bilde');

% Og "label" på aksene
xlabel('x-aksen');
ylabel('y-aksen');

%% Gråtonebilder
% Et fargebilde kan enkelt gjøres til en gråtonebilde ved å bruke
% kommandoen "rgb2gray"
bilde_gratone = rgb2gray(bilde_ole_m);

figure(5)
imshow(bilde_gratone,[])

%Men pass på imagesc kommandoen, den kan gjerne velge et annet colormap
figure(6)
imagesc(bilde_gratone);
axis image

% Som kanskje ikke var helt det vi tenkte, men vi kan sørge for at matlab
% bruker gråtoner ved å spesifisere "colormap gray"
figure(7)
imagesc(bilde_gratone);
axis image
colormap gray

% Vi kan også skru på en colorbar for å indikere hvilke farge som har
% hvilke verdi
colorbar

%% Pixel-verdier
% Fra colorbaren kan vi se at bildet inneholder verdier fra 0 til 255,
% dette er relativt vanlig og representasjon av bilder er noe vi skal lære
% mer om i dette kurset. Det er i midlertid en god ting å begynne å tenke
% på bilder som kun en matrise med verdier for hvert pixel i bilde. For
% eksempel kan vi lage oss et enkelt gråtonebilde med verdier fra 0 til 255
% ved å lage en matrise, slik:

bilde_matrise = [0 1 2 3 4 5 6 7 8 9 10 ;
                 0 0 0 0 100 100 100 0 0 0 0 ;
                 0 0 0 100 0 0 0 100 0 0 0 ;
                 0 0 100 0 0 0 0 0 100 0 0 ;
                 0 100 0 0 200 0 200 0 0 100 0 ;
                 0 100 0 0 0 0 0 0 0 100 0 ;
                 0 100 0 0 50 0 50 0 0 100 0 ;
                 0 0 100 0 0 50 0 0 100 0 0 ;
                 0 0 0 100 0 0 0 100 0 0 0 ;
                 0 0 0 0 100 100 100 0 0 0 0 ;
                 255 254 253 252 251 250 249 248 247 246 245]

            
figure(8)
imagesc(bilde_matrise);
colormap gray
axis image
colorbar

% Her ser vi feks at alle verdiene fra 0 til 10 er svarte, og vanskelig å
% se forskjell på, og verdiene fra 245-255 er tilsvarende lyse.

% Dette var alt for nå, og tenkt som en myk intro til MATLAB, bilder og
% hvordan jeg tenker å bruke nettsiden inf2310.olemarius.net