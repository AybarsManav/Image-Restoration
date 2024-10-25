function wiener_filtered_image = wienerFilter(noisy_motion_blurred_image, deformation_transfer_function, k, plot)

H = deformation_transfer_function;

fft_motion_blurred_noisy =  fftshift(fft2(noisy_motion_blurred_image));
wiener_transfer_func = (H .* conj(H) ./ (H .* conj(H) + k)) ./ H;
fft_estimate = wiener_transfer_func .* fft_motion_blurred_noisy;
wiener_filtered_image = real(ifft2(fftshift(fft_estimate), "symmetric"));

if plot
    figure;
    subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred_noisy)));
    subplot(2, 2, 2); imagesc(log(1 + abs(wiener_transfer_func)));
    subplot(2, 2, 4); imagesc(log(1 + abs(fft_estimate)));
    subplot(2, 2, 3); imshow(wiener_filtered_image, []);
end

end

