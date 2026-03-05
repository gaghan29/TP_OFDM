function DSP = bartlett(x, NFFT)
    % DSP calculée avec la méthode de Bartlett. 
    % Découpage du signal en K segments, chacun de taille égal.
    N = length(x);
    reste = mod(N, NFFT);
    if reste ~= 0
        x = [x, zeros(1, NFFT - reste)];   % 0-padding si pas le bon nombre de points
    end
    window = ones(1, NFFT);
    U = 1; 
    K = length(x) / NFFT;
    DSP_seg = zeros(K,NFFT);
    for i=1:K
        seg = x((i-1)*NFFT+1:i*NFFT);
        %DSP_seg(i,:) = (abs(fft(seg)).^2)/NFFT;
        DSP_seg(i,:) = (abs(fft(seg.*window)).^2) / (NFFT * U);
    end
    DSP = mean(DSP_seg,1);
end

