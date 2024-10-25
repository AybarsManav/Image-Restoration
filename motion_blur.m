im = imread("cameraman.tif");
im  = double(im);

%% Generate Motion Blur Transfer Function
a = 5;
b = 5;
T = 1;

u = 0:255; u = (u - 128) / 256;
v = 0:255; v = (v - 128) / 256;

[U, V] = meshgrid(u, v);

H = T ./ (pi * (U * a + V * b)) .* ...
    sin(pi * (U * a + V * b)) .* ...
    exp(-1i * pi * (U * a + V * b));

% Before applying the inverse filter it is important to handle values of H
% close to 0.
H(isnan(H)) = T;
H(abs(H) < 1e-4) = T;
%% Visualize the magnitude and the phase
figure;
subplot(2, 1, 1); imagesc((abs(H))); title("Magnitude of H");
subplot(2, 1, 2); imagesc((angle(H))); title("Phase of H");

%% Apply Motion Blur to the image in the frequency domain
fft_im = fftshift(fft2(im));
fft_h = H;
fft_multiplied = fft_im .* fft_h;
ifft_multipled = ifft2(fftshift(fft_multiplied), "symmetric");

%% Visualize Motion Blurred Image
figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_im)));
subplot(2, 2, 2); imagesc(log(1 + abs(fft_h)));
subplot(2, 2, 3); imshow(real(ifft_multipled), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_multiplied)));

%% Apply Inverse Filtering


motion_blurred_image = abs(ifft_multipled); % If we do not take the real of it everything works perfectly.
fft_motion_blurred = fftshift(fft2(motion_blurred_image));
fft_filtered = fft_motion_blurred ./ H;
filtered = ifft2(fftshift(fft_filtered), "symmetric");


%%
figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred)));
subplot(2, 2, 2); imagesc(log(1 + abs(H)));
subplot(2, 2, 3); imshow(real(filtered), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_filtered)));

%%
inverse_filtered_mse = computeMSE(im, real(filtered))
inverse_filtered_snr = computeSNR_db(im, filtered)


%% Add noise as well
mean = 127.5;
varience = 1;
noise = mean + sqrt(varience) * randn(size(im));

noisy_motion_blurred_image = motion_blurred_image + noise;

%% Apply Inverse Filtering

fft_motion_blurred = fftshift(fft2(noisy_motion_blurred_image));
fft_filtered = fft_motion_blurred ./ H;
filtered = ifft2(fftshift(fft_filtered), "symmetric");

%% Visualize

figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred)));
subplot(2, 2, 2); imagesc(log(1 + abs(H)));
subplot(2, 2, 3); imshow(real(filtered), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_filtered)));

inverse_filtered_noisy_mse = computeMSE(im, real(filtered))
inverse_filtered_noisy_snr = computeSNR_db(im, filtered)

%% Wiener Filtering of noisy image
k = 2;
fft_motion_blurred_noisy =  fftshift(fft2(noisy_motion_blurred_image));
wiener_transfer_func = (H .* conj(H) ./ (H .* conj(H) + k)) ./ H;
fft_estimate = wiener_transfer_func .* fft_motion_blurred_noisy;
wiener_filtered_image = ifft2(fftshift(fft_estimate), "symmetric");

figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred_noisy)));
subplot(2, 2, 2); imagesc(log(1 + abs(wiener_transfer_func)));
subplot(2, 2, 3); imshow(real(wiener_filtered_image), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_estimate)));

wiener_filtered_mse = computeMSE(im, real(wiener_filtered_image));
wiener_filtered_snr = computeSNR_db(im, real(wiener_filtered_image));


%% Wiener Filtering of noisy image
k = 1;
fft_motion_blurred =  fftshift(fft2(motion_blurred_image));
wiener_transfer_func = (H .* conj(H) ./ (H .* conj(H) + k)) ./ H;
fft_estimate = wiener_transfer_func .* fft_motion_blurred;
wiener_filtered_image = ifft2(fftshift(fft_estimate), "symmetric");

figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred)));
subplot(2, 2, 2); imagesc(log(1 + abs(wiener_transfer_func)));
subplot(2, 2, 3); imshow(real(wiener_filtered_image), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_estimate)));

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