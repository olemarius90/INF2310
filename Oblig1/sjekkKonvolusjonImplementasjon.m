function [correct] = sjekkKonvolusjonImplementasjon(function_name)
%SJEKKKONVOLUSJONIMPLEMENTASJON Summary of this function goes here
%   Detailed explanation goes here

I = imread('pout.tif');

disp(['Sjekker din funksjon ',function_name,'() mot imfilter med kernel:\n']);
kernel = [-1 -sqrt(2) -1; 0 0 0; 1 sqrt(2) 1]
m_res = imfilter(I,kernel,'replicate');
eval(['res = ',function_name,'(I,kernel);'])
correct1 = sum(sum(uint8(res)-uint8(m_res))) == 0;
if correct1
    disp('Dette gikk bra...');
else
    disp('Dette gikk ikke bra...')
end

disp(['Sjekker din funksjon ',function_name,'() mot imfilter med kernel:\n']);
kernel = [1 2 0 -2 -1; 4 8 0 -8 -4; 6 12 0 -12 -6; 4 8 0 -8 -4; 1 2 0 -2 -1]
%kernel = [0 0 0; 1 0 1; 0 0 0];
m_res = imfilter(I,kernel,'replicate','conv','same');
eval(['res = ',function_name,'(I,kernel);'])

correct2 = sum(sum(uint8(res)-m_res)) == 0;
if correct2
    disp('Dette gikk bra...');
else
    disp('Dette gikk ikke bra...')
end

correct = correct1 && correct2;

if correct
    disp('****************  RIKTIG *************************');
    disp('Din implementasjon av konvolusjon virker korrekt.');
else
    disp('****************  FEIL  ************************');
    disp('Din implementasjon av konvolusjon er ikke korrekt.');
end

end

