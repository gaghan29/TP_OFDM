clear
clc
close all

% Implémentation modulation/démodulation OFDM

% Chaine de communication OFDM (Wifi, 4G, 4G Edge, 5G, ...)

% Version avec intervalle de garde

%% Paramètres
M = 2; % BPSK
cstl = exp(1j*2*pi*(0:M-1)/M);

N = 128; % Nb de sous porteuses
K = 500; % Nb de symboles par sous-porteuses

nb_classes = 100;
L = 16; % Longueur du préfixe cyclique

% Paramètres bruit blanc
sigma = 0.1;

%% Emetteur 
idx = randi([0, M-1], 1, K*N);
S = cstl(idx + 1); %Stream de symboles

% Démultiplexage (Passage d'un signal série à un signal parallèle)
S_reshaped = reshape(S, N, K); 

% Bloc IFFT pour la modulation
S_ifft = ifft(S_reshaped, N, 1)*sqrt(N);

% Insertion préfixe cylclique
S_Prefixced = S_ifft(N-L+1:end,:); 
Sn_Prefixced = [S_Prefixced; S_ifft];

% Multiplexage (Permet de passer à des signaux en série)
x = Sn_Prefixced(:).';

%% Canal 

% Fonction de transfert du canal (Canal complexe de variance 1/L)
%h=1; %Pour tester que le TEB = 0 couplé à sigma = 0, question (3.2.2)
h = sqrt(1/(2*L)) * (randn(1,L) + 1j*randn(1,L));

% Convolution du signal avec la réponse impulsionnelle du canal
x_canal = conv(x, h);
x_canal = x_canal(1:length(x));

% Génération d'un bruit blanc Gaussien
bruit_blanc = randn(1, (N+L)*K) * sqrt(sigma);

% Signal bruité
x_bruite = x_canal + bruit_blanc;

%% Récepteur
% Signal recu en série

% Démultiplexage (Passage d'un signal série à un signal parallèle)
R_reshaped = reshape(x_bruite, (N+L), K);

% On enleve le préfixe cyclique
R_sansPrefix = R_reshaped(L+1:end, :); % Matrice sans prefixe cyclique

% Bloc FFT pour la démodulation
S_fft = fft(R_sansPrefix, N, 1)/sqrt(N);

% Egalistation/Forcage à zéro
H = fft(h,N).';
S_estime = S_fft ./ H;

% Multiplexage (Permet de passer à des signaux en série)
S_recieved = S_estime(:).'; % Symboles reçus

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

%% Affichages

% Constellations

% Affichage constellation côté émetteur 
figure;
plot(real(S), imag(S), 'o');
xlabel('Partie réelle');
ylabel('Partie imaginaire');
ylim([-1.5;1.5]);
title('Constellation BPSK - Emetteur');
axis equal;
grid on;

% Affichage constellation côté récepteur 
figure;
plot(real(S_recieved), imag(S_recieved), 'o');
xlabel('Partie réelle');
ylabel('Partie imaginaire');
title('Constellation BPSK - Recepteur');
axis equal;
grid on;

% Histogrammes

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
