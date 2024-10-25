function MSE = computeMSE(im1, im2)
    n_elements = numel(im1);
    diff2 = (im1 - im2).^2;
    MSE = sum(diff2, "all") / n_elements;
end
