function [p,h,c,c_n] = myHist(img,bits)
%MYHIST Regner ut gråtonehistogrammet til et bilde. Antar her at det er
% 256 (2^8) gråtoner i bildet.
%
%   @input
%   img : Input bildet
%   bits : antall biter
%
%   @output
%   p   : Det normaliserte histogrammet til bildet
%   h   : Histogrammet til bildet
%   c   : Det kumulative histogrammet
%   c_n : Det normalizerte kumulative histogrammet

% Finner størrelsen til bildet
[n,m] = size(img);

%Sjekker hvor mange bits vi bruker
if nargin < 2
    bits = 8;
end

% Preallokerer variable for å holde resultatene
h = zeros(1,2^bits);   % Histogrammet
c = zeros(1,2^bits);   % Det kumulative histogrammet

% Finner histogrammet ved å summere hvor mange piksler som finnes av hver
% gråtone. Dette kan helt klart programmeres mer effektivt, men jeg mener
% denne implementasjonen er relativt intuitiv.
for i = 1:2^bits
    h(i) = sum(sum(uint8(img) == i-1));
end

% Normaliserer histogrammet
p = h/(n*m);

% Legger sammen antallet gråtoner "så langt" for å få det kumulative
% histogrammet som indikerer hvor mange piksler som har gråtone mindre eller
% lik gråtone i.
c(1) = h(1);
for i = 2:2^bits
    c(i) = h(i) + c(i-1);
end

% Normaliserer det kumulative histogrammet
c_n = c/(n*m);
end

