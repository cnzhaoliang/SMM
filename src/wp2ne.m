function result = wp2ne(wp)

    EPS0 = 8.854187817e-12; % F/m
    E = 1.60217653e-19; % C
    ME = 9.10938215e-31; % kg

    result = wp .^ 2 .* EPS0 .* ME ./ E .^ 2; % m^-3

end
