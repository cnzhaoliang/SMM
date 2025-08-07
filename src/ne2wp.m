function result = ne2wp(ne)

    EPS0 = 8.854187817e-12; % F/m
    E = 1.60217653e-19; % C
    ME = 9.10938215e-31; % kg

    result = sqrt(ne .* E .^ 2 ./ EPS0 ./ ME); % rad/s

end
