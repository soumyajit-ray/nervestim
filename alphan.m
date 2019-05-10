function a = alphan(v)
    a = 0.01 .* (v + 55) ./ (1 - exp(-(v + 55) / 10));
end