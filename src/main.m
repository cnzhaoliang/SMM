clear;
clc;
close all;

% 启用并行计算池
if isempty(gcp('nocreate'))
    parpool('local');
end

% ---------------------常量--------------------- %

nemax = 1.2e19; la = 0.05; sigma1 = 0.03; sigma2 = 0.465; % 峰值密度(m^-3), 峰值位置(m), 高斯展宽1(m), 高斯展宽2(m)
ve = 0; % 平均碰撞频率(rad/s)
zmax = 2; % 鞘套厚度(m)
k = 100; % 细分比例

% ---------------------变量--------------------- %

wFunc = @(f) (2 * pi * f); lambdaFunc = @(f) (2.99792458e8 ./ f);
neFunc = @(z) arrayfun(@(zi) DGNE(zi, nemax, la, sigma1, sigma2), z);
nFunc = @(f, z) sqrt(wv2epsr(wFunc(f), ne2wp(neFunc(z)), ve));

% ---------------------主计算--------------------- %

figure("Name", "Plasma Density Profile");
z = linspace(0, zmax, 1000);
ne = neFunc(z);
semilogy(z, ne, 'b', 'LineWidth', 2);
xlabel('$z(\rm{m})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);
ylabel('$n_{\rm{e}}(\rm{m^{-3}})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16, 'Rotation', 0);
grid on;

f = (20:0.01:30).' * 1e9; % Hz
R = zeros(length(f), 1);

waiter = ParWaiter(length(f), 0.01);
D = parallel.pool.DataQueue;
afterEach(D, @(~) waiter.update());

parfor i = 1:length(f)
    [points, sizes] = getGrid(0, zmax, lambdaFunc(f(i)) / k);
    n = [1, nFunc(f(i), points), 1];
    [R(i), ~] = SMM(n, sizes, f(i));
    send(D, []);
end

waiter.close();

figure("Name", "Reflection Coefficient");
plot(f / 1e9, real(R), 'r', 'LineWidth', 2); hold on;
plot(f / 1e9, imag(R), 'b', 'LineWidth', 2); hold on;
legend('Real Part', 'Imaginary Part', 'Location', 'Best', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('$f (\rm{GHz})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('$\rm{R}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12, 'Rotation', 0);
ylim([-1, 1]);
yticks((-1:0.2:1));
grid on;

% ---------------------函数--------------------- %

function [points, sizes] = getGrid(xmin, xmax, min_size)
    num = ceil((xmax - xmin) / min_size);
    size = (xmax - xmin) / num;
    points = linspace(xmin + size / 2, xmax - size / 2, num); % 网格点
    sizes = size * ones(1, num); % 网格间距
end

function result = DGNE(x, nemax, la, sigma1, sigma2)
    result = ...
        (x < la) * (nemax * exp(- (x - la) ^ 2 / (2 * sigma1 ^ 2))) + ...
        (x > la) * (nemax * exp(- (x - la) ^ 2 / (2 * sigma2 ^ 2))) + ...
        (x == la) * nemax;
end
