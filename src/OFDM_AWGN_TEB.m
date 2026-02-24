clear
clc
close all

%% Paramètres 

% Constelation
M = 2; % BPSK
cstl = exp(1j*2*pi*(0:M-1)/M);

N = 128; % Nb de sous porteuses
K = 500; % Nb de symboles par sous-porteuses
sigma_symbole = 1;

% Paramètre SNR
SNR_dB = 0:1:10; %Eb/N0 en dB
TEB_sim = zeros(size(SNR_dB)); %Pour collecter les 10 erreurs binaires


%% Main
for i = 1:length(SNR_dB)

    nb_erreur = 0;
    nb_bits_total = 0;

    while nb_erreur <100
        
        %% Emetteur 
        idx = randi([0, M-1], 1, K*N);
        S = cstl(idx + 1); %Stream de symboles
        
        % Démultiplexage (Passage d'un signal série à un signal parallèle)
        S_reshaped = reshape(S, N, K); 
        
        % Bloc IFFT pour la modulation
        S_ifft = ifft(S_reshaped, N, 1)*sqrt(N);
        
        % Multiplexage (Permet de passer à des signaux en série)
        x = S_ifft(:).'; 


        %% Canal 

        % Génération d'un bruit blanc Gaussien
        sigma_bruit = 1 / (10^(SNR_dB(i)/10)); % Valeur de sigma_bruit en fonction de la valeur du SNR (dB)
        bruit_blanc = randn(1, N*K) * sqrt(sigma_bruit);
        
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
        
        %% Bloc de décision
        S_decode = zeros(1, K*N); % vecteur ligne
        
        for k= 1:(K*N)
            if real(S_recieved(k))> 0 
                S_decode(k) = 0; % BPSK: 1 
            else
                S_decode(k) = 1; % BPSK: 0
            end
        end

        % Compter les erreurs
        nb_erreur = nb_erreur + sum(S_decode ~= idx); % Permet de compter les bits d'erreur
        nb_bits_total = nb_bits_total + K*N;

    end
    
    %Mise en place du TEB pour chaque valeur du SNR
    TEB_sim(i) = nb_erreur / nb_bits_total;
end


%% Affichages

%% Courbe théorique (modulation binaire -> Q-function)
SNR_theo = 10.^(SNR_dB/10);
TEB_theo = 0.5 * erfc(sqrt(SNR_theo));

%% Tracé
figure;
semilogy(SNR_dB, TEB_theo, 'r--', 'LineWidth', 2, 'DisplayName', 'Théorie');
hold on;
semilogy(SNR_dB, TEB_sim, 'bo-', 'LineWidth', 2, 'DisplayName', 'Simulation');
hold on;
xlabel('SNR (dB)');
ylabel('TEB');
legend('show');
title('Evolution du TEB en fonction du SNR des sous porteuses');
grid on;