%% Ideelt lavpass filter
clear all
close all

f = double(imread('cameraman.tif'));
[M,N] = size(f);
H = zeros(M,N);
D = zeros(M,N);
D0 = 0.2;
%Ideelt lavpassfilter
for i = 1:M
    for j = 1:N
        D(i,j) = sqrt(((i-floor(M/2 + 1))/(M/2))^2+(((j-floor(N/2+1))/(N/2))^2));
        if D(i,j) <= D0;
            H(i,j) = 1;
        end
    end
end

figure(3)
subplot(211)
imagesc(D)
subplot(212)
imshow(H,[]);

F_window = fftshift(fft2(f));
g = real(ifft2(ifftshift(F_window.*H)));

figure(4)
subplot(221)
imshow(f,[]);
subplot(222)
imagesc(abs(log(F_window+1)));
subplot(223)
imshow(g,[]);
subplot(224)
imagesc(log(abs(fftshift(fft2(g)))+1));

%% Butterworth
n = 5;
H_butter = 1./(1+(D./D0).^(2*n));
g_butter = real(ifft2(ifftshift(F_window.*H_butter)));

figure(5)
subplot(311)
imagesc(H)
subplot(312)
imagesc(H_butter)

figure(6)
subplot(221)
imshow(f,[]);
subplot(222)
imagesc(abs(log(F_window+1)));
subplot(223)
imshow(g_butter,[]);
subplot(224)
imagesc(log(abs(fftshift(fft2(g_butter)))+1));

%% Gaussisk lavpassfilter
H_gauss = exp((-D.^2)/(2*D0^2));
g_gauss = real(ifft2(ifftshift(F_window.*H_gauss)));


figure(5)
subplot(313)
imagesc(H_gauss);

figure(7)
subplot(221)
imshow(f,[]);
subplot(222)
imagesc(abs(log(F_window+1)));
subplot(223)
imshow(g_gauss,[]);
subplot(224)
imagesc(log(abs(fftshift(fft2(g_gauss)))+1));


%% 
k_size = 40;
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


figure(8)
subplot(311)
imagesc(H_gauss)
subplot(312)
imagesc(filter);
subplot(313)
imagesc(kernel_gauss)


figure(9)
subplot(311)
imagesc(kernel_ideal)
subplot(312)
imagesc(kernel_butter)
subplot(313)
imagesc(kernel_gauss)

g = imfilter(f,kernel,'conv','symmetric')

figure(10)
subplot(311)
imshow(f)
subplot(312)
imshow(g,[])
subplot(313)
imshow(g_gauss,[])

%% Notch filter example
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
subplot(222)
imshow(log(abs(fftshift(F))+1),[]);axis on;


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
subplot(224); imshow(real(ifft2(G)),[]); axis on;


%%

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
%% Filtering with zero-padding.
% With zero-padding.
fp = [f zeros(M,N) ; zeros(M,N) zeros(M,N)];
Fp = fft2(fp);
Mp = 2*M;
Np = 2*N;
up = u.*2-1;
vp = v.*2-1;
figure(10); imshow(fp,[]); axis on;
figure(11); imshow(log( abs(fftshift(Fp)) + 1), [0 max(log( abs(Fp(:)) + 1))]); axis on;

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
figure(12); imshow(log( abs(fftshift(Gp)) + 1), [0 max(log( abs(Gp(:)) + 1))]); axis on;
figure(13); imshow(gp,[]); axis on;
figure(14); ax1 = subplot(232); imshow(gp(1:M,1:N),[]); axis on; title('Zero padded');
ax2 = subplot(231); imshow(real(ifft2(G)),[]); axis on; title('Circular index')
subplot(233); imshow(real(ifft2(G))-gp(1:M,1:N),[]);
linkaxes([ax1 ax2],'xy');

%%
filter = fftshift(real(ifft2(ifftshift(H))));

figure(1)
subplot(121)
imagesc(log(H+1));
subplot(122)
imagesc(filter);

res = imfilter(f,filter,'conv','circular');
res_padded = imfilter(f,filter,'conv');

figure(14)
subplot(234)
imshow(res,[])
subplot(235)
imshow(res_padded,[])
subplot(236)
imshow(res-res_padded,[])
%% Vindusfunksjon
% Test: Is it useful to first use a window function in this case, i.e. when
%       we use zero-padding? Here for a Tukey window function:
% Result: No, the result is worse.
f = double(imread('car2.png'));
[M,N] = size(f);
alpha = 0.3; w = (tukeywin(N,alpha)*tukeywin(M,alpha)')';
fw = f .* w;
F_window = fftshift(abs(fft2(fw)));
F = fftshift(abs(fft2(f)));
%fwp = [fw zeros(M,N) ; zeros(M,N) zeros(M,N)];
%Fwp = fft2(fwp);
%figure(15); imshow(fwp); axis on;
figure(15);
subplot(131)
imshow(f,[]);
subplot(132)
imshow(w);
subplot(133)
imshow(fw,[]);


figure(16); subplot(211); imshow(log( abs(F_window) + 1), [0 max(log( abs(F_window(:)) + 1))]); axis on;
subplot(212); imshow(log( abs(F) + 1), [0 max(log( abs(F(:)) + 1))]); axis on;

% Filtering.
% Gp = Fwp.*ifftshift(Hp);
% gp = real(ifft2(Gp));
% figure(17); imshow(log( abs(fftshift(Gp)) + 1), [0 max(log( abs(Gp(:)) + 1))]); axis on;
% figure(18); imshow(gp); axis on;
% figure(19); imshow(gp(1:M,1:N)); axis on;
% figure(20); imshow(gp(1:M,1:N) ./ w); axis on;

%%
