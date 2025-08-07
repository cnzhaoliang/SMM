function [R, T] = SMM(n, d, f)
    % Calculate reflection and transmission coefficients for multilayer media by the use of the Scattering Matrix Method (SMM).
    % TEM waves incident perpendicularly to the boundary and all media are non-magnetized(μ = μ₀).
    %
    % E1 = E1f * exp(-1j * k * z) + E1b * exp(1j * k * z), H1 = 1 / eta1 * (E1f * exp(-1j * k * z) - E1b * exp(1j * k * z))
    % E2 = E2f * exp(-1j * k * z) + E2b * exp(1j * k * z), H2 = 1 / eta2 * (E2f * exp(-1j * k * z) - E2b * exp(1j * k * z))
    %
    % Boundary condition: E1(0) = E2(0), H1(0) = H2(0)
    % Scattering Matrix:  [E1b; E2f] = [S11, S12; S21, S22] * [E1f; E2b]
    %
    %   Inputs:
    %       n - Refractive index vector     [n₁, n₂, ..., nₙ]
    %       d - Layer thickness vector      [d₁, d₂, ..., dₙ₋₂] (m)
    %       f - Frequency (Hz)
    %
    %   Outputs:
    %       R - Reflection coefficient
    %       T - Transmission coefficient
    %
    %   Example:
    %       n = [1, 1.1 - 1j * 0.02, 1];    % Air-Dielectric-Air
    %       d = 0.01;                       % 1cm
    %       f = 10e9;                       % 10GHz
    %       [R, T] = SMM(n, d, f);

    arguments
        n {mustBeNumeric, mustBeVector}
        d {mustBeNumeric, mustBeVector, mustBeReal, mustBePositive}
        f {mustBeNumeric, mustBeReal, mustBePositive}
    end

    assert(length(n) == length(d) + 2, ...
    'Invalid parameter vector lengths.');

    % k = beta - 1j * alpha = omega * n
    k = 2 * pi * f / 2.99792458e8 .* n; % rad/m

    % Find elements in k vector where abs(k) < eps and replace them with eps
    k(abs(k) < eps) = eps;

    % Calculate the scattering matrix
    S = S4Boundary(n(1), n(2));

    for i = 1:length(d)

        S = RedhefferStarProduct(S, S4Area(k(i + 1), d(i)));
        S = RedhefferStarProduct(S, S4Boundary(n(i + 1), n(i + 2)));

    end

    % Calculate the reflection and transmission coefficients
    [R, T] = deal(S(1, 1), S(2, 1));

end

function S = S4Boundary(n1, n2)

    S = 1 / (n1 + n2) * [n1 - n2, 2 * n2; 2 * n1, n2 - n1];

end

function S = S4Area(k, d)

    expTerm = exp(-1j * k * d);
    S = [0, expTerm; expTerm, 0];

end

function SAB = RedhefferStarProduct(SA, SB)

    % Redheffer star product for 2x2 matrices.
    % Combines two scattering matrices to form a overall scattering matrix.
    % SA and SB are two different scattering matrices.
    % SAB is the combined scaterring matrix.

    temp1 = 1 / (1 - SB(1, 1) * SA(2, 2));
    temp2 = 1 / (1 - SA(2, 2) * SB(1, 1));
    SAB = zeros(2, 2);
    SAB(1, 1) = SA(1, 1) + SA(1, 2) * temp1 * SB(1, 1) * SA(2, 1);
    SAB(1, 2) = SA(1, 2) * temp1 * SB(1, 2);
    SAB(2, 1) = SB(2, 1) * temp2 * SA(2, 1);
    SAB(2, 2) = SB(2, 2) + SB(2, 1) * temp2 * SA(2, 2) * SB(1, 2);

end
