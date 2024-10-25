function inverse_filtered = inverseFilter(motion_blurred_image, fft_deformation, plot)

fft_motion_blurred = fftshift(fft2(motion_blurred_image));
fft_filtered = fft_motion_blurred ./ fft_deformation;
inverse_filtered = real(ifft2(fftshift(fft_filtered), "symmetric"));

if plot    
    figure;
    subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred))); title("FFT of Motion Blurred Image");
    subplot(2, 2, 2); imagesc(log(1 + abs(fft_deformation))); title("FFT of Motion Kernel");
    subplot(2, 2, 3); imagesc(log(1 + abs(fft_filtered))); title("FFT of Filtered Image");
    subplot(2, 2, 4); imshow(inverse_filtered, []); title("Restored Image");
end

end
