function [motion_blurred_image, motion_blur_kernel_spectrum] = applyLinearMotionBlur(a, b, T, image, plot)
% Parameters:
% a : Horizontal velocity
% b : Vertical velocity
% T : Exposure time
% image : Input image
% plot : Plot the images
% Output parameters:
% motion_blurred_image : Motion blur applied image.
% motion_blur_kernel_spectrum : Spectrum of motion blur.

image = double(image);
width = size(image, 2);
height = size(image, 1);

u = 0:width-1; u = (u - ceil(width / 2)) / width;
v = 0:height-1; v = (v - ceil(height / 2)) / height;

[U, V] = meshgrid(u, v);

H = T ./ (pi * (U * a + V * b)) .* ...
    sin(pi * (U * a + V * b)) .* ...
    exp(-1i * pi * (U * a + V * b));

% Take the limit of H as (U * a + V * b) goes to 0:
H(isnan(H)) = T; 
H(abs(H) < 1e-4) = T;

motion_blur_kernel_spectrum = H;


% Apply Motion Blur to the image in the frequency domain
fft_im = fftshift(fft2(image));
fft_h = H;
fft_multiplied = fft_im .* fft_h;
ifft_multipled = ifft2(fftshift(fft_multiplied), "symmetric");

motion_blurred_image = real(ifft_multipled);

if plot
% Visualize the magnitude and the phase
figure;
subplot(2, 1, 1); imagesc((abs(H))); title("Magnitude of H");
subplot(2, 1, 2); imagesc((angle(H))); title("Phase of H");

figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_im))); title("FFT of Image");
subplot(2, 2, 2); imagesc(log(1 + abs(fft_h))); title("FFT of Motion Blur");
subplot(2, 2, 3); imagesc(log(1 + abs(fft_multiplied))); title("FFT of Multiplied Spectrums");
subplot(2, 2, 4); imshow(real(ifft_multipled), []); title("Motion Blurred Image");
end

end


