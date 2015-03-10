function [correct] = sjekkKonvolusjonImplementasjon(function_name)
%SJEKKKONVOLUSJONIMPLEMENTASJON Funksjonen sjekker om implementasjonen av
%   din funksjon gir samme resultat som MATLAB's implemetasjon
%   imfilter(I,kernel,'replicate','conv','same')
% 
%   @input:
%   function_name   :   Er navnet på din funksjon. Om din funksjon heter
%                       feks myImfilter skal du kalle denne funksjonen
%                       slik: sjekkKonvolusjonImplementasjon('myImfilter')
%
%   @output:        
%   correct         :   correct = 1 om din funksjon passerer begge testene
%                       correct = 0 om din funksjon ikke passerer begge testene

%   Laster inn et bilde
I = imread('pout.tif');

disp(['Sjekker din funksjon ',function_name,'() mot imfilter med kernel:']);
%   Lager en kernel. Bruker her Frei-Chen operatoren
kernel = [-1 -sqrt(2) -1; 0 0 0; 1 sqrt(2) 1]
m_res = imfilter(double(I),double(kernel),'replicate','conv','same');%Kaller her MATLAB's implementasjon
eval(['res = ',function_name,'(I,kernel);'])%Kaller her deres implementasjon

if size(res) ~= size(I)
    error('Utbildet fra din funksjon har ikke samme størrelse som innbildet');
end

correct1 = sum(sum(res-m_res)) < 1e-10; %Sammenligner resultatene
if correct1
    disp('Dette gikk bra...');
else
    disp('Dette gikk ikke bra...')
end

disp(['Sjekker din funksjon ',function_name,'() mot imfilter med kernel:']);
kernel = [1 2 0 -2 -1; 4 8 0 -8 -4; 6 12 0 -12 -6; 4 8 0 -8 -4; 1 2 0 -2 -1]
m_res = imfilter(double(I),double(kernel),'replicate','conv','same');%Kaller her MATLAB's implementasjon
eval(['res = ',function_name,'(I,kernel);'])         %Kaller her deres implementasjon

if size(res) ~= size(I)
    error('Utbildet fra din funksjon har ikke samme størrelse som innbildet');
end

correct2 = sum(sum(res-m_res)) < 1e-10;          %Sammenligner resultatene
if correct2
    disp('Dette gikk bra...');
else
    disp('Dette gikk ikke bra...')
end

correct = correct1 && correct2;                   %Sjekker om begge ble riktig

if correct
    disp('****************  RIKTIG *************************');
    disp('Din implementasjon av konvolusjon virker korrekt.');
else
    disp('****************  FEIL  ************************');
    disp('Din implementasjon av konvolusjon er ikke korrekt.');
end

end

