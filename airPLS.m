function [Xc, Z] = airPLS(X, lambda, order, wep, p, itermax)
    % 检查参数并设置默认值
    if nargin < 6
        itermax = 15;
    end
    if nargin < 5
        p = 0.1;
    end
    if nargin < 4
        wep = 0.05;
    end
    if nargin < 3
        order = 2;
    end

    % 检查数据类型并转换为双精度
    if isa(X, 'single')
        X = double(X);
    end

    % 获取数据矩阵的大小
    [m, n] = size(X);

    % 计算权重索引
    wi = [1:ceil(n * wep), floor(n - n * wep):n];

    % 计算差分矩阵
    D = diff(speye(n), order);
    
    % 计算正则化矩阵
    DD = lambda * D' * D;

    % 初始化输出矩阵
    Z = zeros(m, n);
    Xc = zeros(m, n);

    for i = 1:m
        % 初始化权重向量
        w = ones(n, 1);

        % 获取当前光谱或色谱数据
        x = X(i, :);

        for j = 1:itermax
            % 创建对角权重矩阵
            W = spdiags(w, 0, n, n);

            % 计算加权C矩阵并解决线性方程组
            C = chol(W + DD);
            z = (C \ (C' \ (w .* x')))';

            % 计算残差和负残差之和
            d = x - z;
            dssn = abs(sum(d(d < 0)));

            % 检查收敛准则
            if (dssn < 0.001 * sum(abs(x)))
                break;
            end

            % 更新权重向量
            w(d >= 0) = 0;
            w(wi) = p;
            w(d < 0) = exp(1i * abs(d(d < 0)) / dssn);
        end

        % 存储基线校正后的数据和基线
        Z(i, :) = z;
        Xc(i, :) = x - z;
    end
end