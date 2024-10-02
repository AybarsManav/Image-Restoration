im = imread("cameraman.tif");

im = padarray(im, [1 1],"replicate" , 'post');
im  = double(im);
%%

a = 6;
b = 6;
T = 1;

u = 0:256; u = (u - 128) / 257;
v = 0:256; v = (v - 128) / 257;

[U, V] = meshgrid(u, v);

H = T ./ (pi * (U * a + V * b) + 1e-12) .* ...
    sin(pi * (U * a + V * b)) .* ...
    exp(-1i * pi * (U * a + V * b));
H(H == 0) = 1;
%%
figure;
subplot(2, 2, 1); imagesc((real(H)));
subplot(2, 2, 2); imagesc((imag(H)));
x = (abs(H));
subplot(2, 2, [3 4]); plot(x(127, :));
%%
fft_im = fftshift(fft2(im));
fft_h = H;
% fft_im = fft2(im);
% fft_h = H;
fft_multiplied = fft_im .* fft_h;
ifft_multipled = ifft2(fftshift(fft_multiplied));
% ifft_multipled = ifft2((fft_multiplied));

%%
figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_im)));
subplot(2, 2, 2); imagesc(log(1 + abs(fft_h)));
subplot(2, 2, 3); imshow(real(ifft_multipled), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_multiplied)));

%%
motion_blurred_image = real(ifft_multipled);
fft_motion_blurred = fftshift(fft2(motion_blurred_image));
fft_filtered = fft_motion_blurred ./ H;
filtered = ifft2(fftshift(fft_filtered));


%%
figure;
subplot(2, 2, 1); imagesc(log(1 + abs(fft_motion_blurred)));
subplot(2, 2, 2); imagesc(log(1 + abs(H)));
subplot(2, 2, 3); imshow(real(filtered), []);
subplot(2, 2, 4); imagesc(log(1 + abs(fft_filtered)));