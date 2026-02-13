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


figure;
plot(S_recieved);
grid on;

%% Bloc de réception
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

