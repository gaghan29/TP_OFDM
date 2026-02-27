clear
clc
close all

% Implémentation modulation/démodulation OFDM

% Chaine de communication OFDM (Wifi, 4G, 4G Edge, 5G, ...)

% Version sans intervalle de garde

%% Paramètres
M = 2; % BPSK
cstl = exp(1j*2*pi*(0:M-1)/M);

N = 128; % Nb de sous porteuses
K = 500; % Nb de symboles par sous-porteuses

nb_classes = 100;

% Paramètres bruit blanc
sigma = 0;

%% Emetteur 
idx = randi([0, M-1], 1, K*N);
S = cstl(idx + 1); %Stream de symboles

% Démultiplexage (Passage d'un signal série à un signal parallèle)
S_reshaped = reshape(S, N, K); 

% Bloc IFFT pour la modulation
S_ifft = ifft(S_reshaped, N, 1)*sqrt(N);

% Multiplexage (Permet de passer à des signaux en série)
x = S_ifft(:).'; 

% Affichage constélation de l'emetteur 
figure;
%plot(real(x));
plot(S);
grid on;

%% Canal 


% Génération d'un bruit blanc Gaussien
bruit_blanc = randn(1, N*K) * sqrt(sigma);

% Signal bruité
x_bruite = x + bruit_blanc;


%% Récepteur
% Signal recu en série

% Démultiplexage (Passage d'un signal série à un signal parallèle)
R_reshaped = reshape(x_bruite, N, K);

% Bloc FFT pour la démodulation
S_fft = fft(R_reshaped, N, 1)/sqrt(N);

% Multiplexage (Permet de passer à des signaux en série)
S_recieved = S_fft(:).'; % Symboles reçus

% Affichage constélation de l'emetteur
figure;
plot(S_recieved);
grid on;

%% Bloc de décision
S_decode = zeros(1, K*N);    % vecteur ligne

for i= 1:(K*N)
    if real(S_recieved(i))> 0 
        S_decode(i) = 0; % BPSK: 1 for positive real part
    else
        S_decode(i) = 1; % BPSK: 0 for non-positive real part
    end
end


%% TEB (Taux d'erreur binaire)
TEB = sum(S_decode ~= idx)/(K*N);

%% Histogrammes

% Histogramme partie réelle trame OFDM
figure;
histogram(real(x),nb_classes);
xlabel('Symbole');
ylabel('Nombre');
title('Histogramme partie réelle trame OFDM');
grid on;

% Histogramme partie imaginaire trame OFDM
figure;
histogram(imag(x),nb_classes);
xlabel('Symbole');
ylabel('Nombre');
title('Histogramme partie imaginaire trame OFDM');
grid on;


%% Fonctions d'autocorélation
trame_OFDM = x;

[R_real, lags_real]= xcorr(real(trame_OFDM), 'biased');
[R_imag, lags_imag] = xcorr(imag(trame_OFDM), 'biased');

[R_inter, lags_inter] = xcorr(real(trame_OFDM), imag(trame_OFDM), 'biased');

% Normalisation et passage en échelle log
R_real = abs(R_real)/max(R_real);
R_realdB = 10*log10(R_real + 1e-10);

R_imag = abs(R_imag)/max(R_imag);
R_imagdB = 10*log10(R_imag + 1e-10);

R_inter = abs(R_inter)/max(R_inter);
R_interdB = 10*log10(R_inter + 1e-10);

figure;
subplot(3,1,1)
grid on;
plot(lags_real, R_realdB, 'r');
xlabel('Lags')
ylabel('Autocorrélation normalisé (dB)')
title('Autocorrélation normalisé de la partie réelle')

subplot(3,1,2)
grid on;
plot(lags_imag, R_imagdB, 'g');
xlabel('Lags')
ylabel('Autocorrélation normalisé (dB)')
title('Autocorrélation normalisé de la partie imaginaire')

subplot(3,1,3)
grid on;
plot(lags_inter, R_interdB, 'b');
xlabel('Lags')
ylabel('Autocorrélation normalisé (dB)')
title('Intercorrélation normalisé de la partie réelle et imaginaire')


