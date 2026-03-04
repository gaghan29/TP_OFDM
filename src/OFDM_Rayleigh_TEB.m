clear
clc
close all

%% Paramètres 

% Constelation
M = 2; % BPSK
cstl = exp(1j*2*pi*(0:M-1)/M);

N = 128; % Nb de sous porteuses
K = 500; % Nb de symboles par sous-porteuses

L = 16; % Longueur du préfixe cyclique

% Paramètre SNR
SNR_dB = 0:1:15; %Eb/N0 en dB
TEB_sim = zeros(size(SNR_dB)); %Pour collecter les 15 erreurs binaires


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

        % Insertion préfixe cylclique
        S_Prefixced = S_ifft(N-L+1:end,:); 
        Sn_Prefixced = [S_Prefixced; S_ifft];

        % Multiplexage (Permet de passer à des signaux en série) 
        x = Sn_Prefixced(:).';

        %% Canal 

        % Fonction de transfert du canal (Canal complexe de variance 1/L)
        h = sqrt(1/(2*L)) * (randn(1,L) + 1j*randn(1,L));
        
        % Convolution du signal avec la réponse impulsionnelle du canal
        x_canal = conv(x, h);
        x_canal = x_canal(1:length(x));

        % Génération d'un bruit blanc Gaussien
        sigma_bruit = 1 / (10^(SNR_dB(i)/10)); % Valeur de sigma_bruit en fonction de la valeur du SNR (dB)
        % On suppose que sigma_symbole^2E[|H(n)|^2] = 1
        bruit_blanc = randn(1, (N+L)*K) * sqrt(sigma_bruit);
        
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

% Courbe théorique (modulation binaire -> Q-function)
SNR_theo = 10.^(SNR_dB/10);
TEB_theo = berfading(SNR_dB,'psk',2,1);

% Tracé
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