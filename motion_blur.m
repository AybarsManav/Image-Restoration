im = imread("cameraman.tif");
im  = double(im);

%% Generate Motion Blur Transfer Function
a = 15;
b = 15;
T = 1;

u = 0:255; u = (u - 128) / 256;
v = 0:255; v = (v - 128) / 256;

[U, V] = meshgrid(u, v);

H = T ./ (pi * (U * a + V * b) + 1e-12) .* ...
    sin(pi * (U * a + V * b)) .* ...
    exp(-1i * pi * (U * a + V * b));

% Before applying the inverse filter it is important to handle values of H
% close to 0.

H(abs(H) < 1e-6) = 1;
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

%% MSE computation

function MSE = computeMSE(im1, im2)
    n_elements = numel(im1);
    diff2 = (im1 - im2).^2;
    MSE = sum(diff2, "all") / n_elements;
end
