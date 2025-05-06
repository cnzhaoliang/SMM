function result = wv2epsr(w, wp, ve)

    result = 1 - wp .^ 2 ./ (w .^ 2 + ve .^ 2) - ...
        1j .* ve ./ w .* wp .^ 2 ./ (w .^ 2 + ve .^ 2);

end
