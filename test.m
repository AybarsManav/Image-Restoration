im = imread("cameraman.tif");
im  = double(im);
im = padarray(im, [1 1], "replicate", "pre");

fft_im = fftshift(fft2((im)));

figure;
imagesc(log(1 + abs(fft_im)));

%%

% Check for symmetry in the real part
real_part = real(H);
symmetric_real = isequal(real_part, rot90(real_part, 2));

% Check for antisymmetry in the imaginary part
imag_part = imag(H);
antisymmetric_imag = isequal(imag_part, -rot90(imag_part, 2));

disp(['Real part symmetric: ', num2str(symmetric_real)]);
disp(['Imaginary part antisymmetric: ', num2str(antisymmetric_imag)]);

% Compare positive and negative frequency components
diff_real = real_part - fliplr(flipud(real_part));
diff_imag = imag_part + fliplr(flipud(imag_part));

max_diff_real = max(abs(diff_real(:)));
max_diff_imag = max(abs(diff_imag(:)));

disp(['Max real part difference: ', num2str(max_diff_real)]);
disp(['Max imaginary part difference: ', num2str(max_diff_imag)]);

%%
% Check for symmetry in the real part
real_part = real(fft_im);
symmetric_real = isequal(real_part, rot90(real_part, 2));

% Check for antisymmetry in the imaginary part
imag_part = imag(fft_im);
antisymmetric_imag = isequal(imag_part, -rot90(imag_part, 2));

disp(['Real part symmetric: ', num2str(symmetric_real)]);
disp(['Imaginary part antisymmetric: ', num2str(antisymmetric_imag)]);

% Compare positive and negative frequency components
diff_real = real_part - fliplr(flipud(real_part));
diff_imag = imag_part + fliplr(flipud(imag_part));

max_diff_real = max(abs(diff_real(:)));
max_diff_imag = max(abs(diff_imag(:)));

disp(['Max real part difference: ', num2str(max_diff_real)]);
disp(['Max imaginary part difference: ', num2str(max_diff_imag)]);