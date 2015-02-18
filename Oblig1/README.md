#	Oblig 1 i INF2310

På denne siden ligger teksten til Oblig1 i INF2310 i filen Oblig1.pdf. 
I filen sjekkKonvolusjonImplementasjon.m har jeg skrevet en funksjon som tester
deres implementasjon mot MATLAB's implementasjon imfilter( ... , ... ,'replicate','conv','same')

Min implementasjon er i funksjonen myImfilter, og jeg kaller da denne filen slik: sjekkKonvolusjonImplementasjon('myImfilter');

>> sjekkKonvolusjonImplementasjon('myImfilter');
Sjekker din funksjon myImfilter() mot imfilter med kernel:\n

kernel =

   -1.0000   -1.4142   -1.0000
         0         0         0
    1.0000    1.4142    1.0000

Dette gikk bra...
Sjekker din funksjon myImfilter() mot imfilter med kernel:\n

kernel =

     1     2     0    -2    -1
     4     8     0    -8    -4
     6    12     0   -12    -6
     4     8     0    -8    -4
     1     2     0    -2    -1

Dette gikk bra...
****************  RIKTIG *************************
Din implementasjon av konvolusjon virker korrekt.


Om din implementasjon ikke er riktig får du beskjed om det.

Gi beskjed til olemarius@olemarius.net om det er feil eller mangler i koden.
