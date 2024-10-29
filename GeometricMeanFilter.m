function gme_filtered_image = GeometricMeanFilter(image,H,alpha,beta,k,plot)
    fft_motion_blurred_noisy =  fftshift(fft2(image));
    GME_filter = ((conj(H)./(H.*conj(H))).^alpha).*((conj(H)./(conj(H).*H + beta*k)).^(1-alpha));
    fft_estimate_gme = GME_filter .* fft_motion_blurred_noisy;
    gme_filtered_image = real(ifft2(fftshift(fft_estimate_gme),'symmetric'));
    if plot
        figure;
        subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred_noisy)));
        title('Frequency Response of the Motion Blurred Image');
        subplot(2, 2, 2); imagesc(log(1 + abs(GME_filter)));
        title('Frequency Response of the Geometric Mean Filter');
        subplot(2, 2, 3); imagesc(log(1 + abs(fft_estimate_gme)));
        title('Frequency Response of the Geometric Mean Filter');
        subplot(2, 2, 4); imshow(gme_filtered_image, []);
        title('Restored Image with Geometric Mean Filter');

     end
    
end