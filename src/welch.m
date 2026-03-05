function DSP = welch(x, NFFT, overlap)
    % Algorithme de Welch: découpage en K segments du signal fenêtré puis
    % moyenne. 
    % Recouvrement entre chaque segments de % overlap
    % x : signal 
    % NFFT: Nombre de points de la DSP finale

    N = length(x);
    recouvrement = floor(overlap * NFFT);  % Nombre de points qui recouvrent un morceau du signal
    ecart = NFFT - recouvrement;           % Ecart entre chaque morceau
    K = floor((N - NFFT) / ecart) + 1;
    window = hamming(NFFT)';
    U = sum(window.^2)/NFFT;
    DSP_seg = zeros(K,NFFT);
    for i=1:K
        seg = x((i-1)*ecart+1:(i-1)*ecart+NFFT);
        seg_window = seg.*window;
        DSP_seg(i,:) = (abs(fft(seg_window)).^2) / (NFFT*U);
    end
    DSP = mean(DSP_seg,1);
end