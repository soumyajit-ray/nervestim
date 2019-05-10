function a = alpham(v)
    a = 0.1 * (v + 40) ./ (1 - exp(-(v + 40) / 10));
end