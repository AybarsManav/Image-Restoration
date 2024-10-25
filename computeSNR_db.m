function SNR_db = computeSNR_db(original_signal, estimated_signal)
    SNR = sum(original_signal.^2, "all") / sum( (original_signal - estimated_signal).^2, "all");
    SNR_db = 10 * log(SNR) / log(10);
end
