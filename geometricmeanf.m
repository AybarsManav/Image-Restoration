function [GME_filter,fft_estimate_gme,GME_image] = geometricmeanf(image,H,alpha,beta,k)
    GME_filter = ((conj(H)./(H.*conj(H))).^alpha).*((conj(H)./(conj(H).*H + beta*k)).^(1-alpha));
    fft_estimate_gme = GME_filter .* image;
    GME_image = ifft2(fftshift(fft_estimate_gme),'symmetric');
    figure;
    subplot(2, 2, 1); imagesc(log(1 + abs(image)));
    title('Frequency Response of the Motion Blurred Image');
    subplot(2, 2, 2); imagesc(log(1 + abs(GME_filter)));
    title('Frequency Response of the Geometric Mean Filter');
    subplot(2, 2, 3); imshow(real(GME_image), []);
    title('Restored Image with Geometric Mean Filter');
    subplot(2, 2, 4); imagesc(log(1 + abs(fft_estimate_gme)));
    title('Frequency Response of the Geometric Mean Filter');
end