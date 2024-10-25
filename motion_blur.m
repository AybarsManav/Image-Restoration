im = imread("cameraman.tif");
im  = double(im);
%%
figure;
imshow(im, []);
%% Visualize Motion Blurred Image
[motion_blurred_image, H] = applyLinearMotionBlur(a, b, T, im, true);

%% Apply Inverse Filtering

filtered = inverseFilter(motion_blurred_image, H, true);

% Evaluate Performance
inverse_filtered_mse = computeMSE(im, real(filtered))
inverse_filtered_snr = computeSNR_db(im, filtered)


%% Add noise as well
mean = 127.5;
varience = 5;
noise = mean + sqrt(varience) * randn(size(im));

noisy_motion_blurred_image = motion_blurred_image + noise;

%% Apply Inverse Filtering

filtered = inverseFilter(noisy_motion_blurred_image, H, true);

inverse_filtered_noisy_mse = computeMSE(im, real(filtered))
inverse_filtered_noisy_snr = computeSNR_db(im, filtered)

%% Wiener Filtering of noisy image
k = 2;
wiener_filtered_image = wienerFilter(noisy_motion_blurred_image, H, k, true);

wiener_filtered_mse = computeMSE(im, real(wiener_filtered_image));
wiener_filtered_snr = computeSNR_db(im, real(wiener_filtered_image));

%% Wiener Filtering of noiseless image
k = 0;
wiener_filtered_image = wienerFilter(motion_blurred_image, H, k, true);

wiener_noisy_filtered_mse = computeMSE(im, real(wiener_filtered_image));
wiener_noisy_filtered_snr = computeSNR_db(im, real(wiener_filtered_image));
%% GEOMETRIC MEAN FILTER with Motion Blur only
alpha = 1;
beta = 1;
k = 10;
[~,~,restored] = geometricmeanf(fft_motion_blurred,H,alpha,beta,k);
gme_filtered_mse = computeMSE(im, real(restored));
gme_filtered_snr = computeSNR_db(im, real(restored));
% GEOMETRIC MEAN FILTER with Motion Blur and Noise
[~,~,restored] = geometricmeanf(fft_motion_blurred_noisy,H,alpha,beta,k);
gme_noisy_filtered_mse = computeMSE(im, real(restored));
gme_noisy_filtered_snr = computeSNR_db(im, real(restored));
%% MSE computation
function MSE = computeMSE(im1, im2)
    n_elements = numel(im1);
    diff2 = (im1 - im2).^2;
    MSE = sum(diff2, "all") / n_elements;
end

function SNR_db = computeSNR_db(original_signal, estimated_signal)
    SNR = sum(original_signal.^2, "all") / sum( (original_signal - estimated_signal).^2, "all");
    SNR_db = 10 * log(SNR) / log(10);
end
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
