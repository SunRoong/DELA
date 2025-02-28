%% Six Masses–Six Springs Simulation with Sine Wave Springs
%
% 描述：
%  本程序模拟一个垂直悬挂的弹簧–质量系统：
%    - 系统有 n 个质量（每个质量 m，默认为 1 kg），由 n 根弹簧（弹簧常数 k，
%      自然长度 L0）依次连接，并悬挂于天花板（ y = 0 ）下。
%    - 重力加速度为 g（向下），系统初始处于静态平衡状态，其各质量位置根据 L0 重新计算。
%    - 从 t = 0 开始，对最下方的质量施加一个额外向下的力 F_ext，使系统偏离平衡。


clear; clc; close all;

%% Adjustable Parameters
n = 6;              % 质量个数
m = 1;              % 每个质量（kg），本模型默认 m=1
k = 3;              % 弹簧常数（N/m）
L0 = 0.5;           % 弹簧自然长度（m）  % 修改 L0 后，初始状态会相应更新
g = 10;             % 重力加速度（m/s²）
F_ext = -6;          % 施加在最下方质量上的额外向下力（N）
animSpeed = 10000;  % 动画速度因子（1：真实时间；>1 加快；<1 减慢）

%% 计算静力平衡状态并设置初始条件
% 调用 init_state 函数，根据当前 L0、g、k 重新计算初始平衡位置
z0 = init_state(n, L0, g, k);
y_eq = z0(1:2:end);  % 提取各质量平衡位置（用于后续箭头绘制等）

%% 模拟时间设置
T = 80;              % 模拟总时长（秒）
tspan = [0 T];

%% 将参数打包到结构体中传递给 ODE 函数
params.n = n;
params.k = k;
params.L0 = L0;
params.g = g;
params.F_ext = F_ext;

%% 使用 ode45 求解系统运动
[t, sol] = ode45(@(t,z) mass_spring_ode(t,z,params), tspan, z0);

%% 动画参数设置
pauseTime = 0.05 / animSpeed;  % 每帧暂停时间

%% 动画演示
fig = figure;
set(fig, 'Name', 'Six Masses–Six Springs Simulation (Sine Wave Springs)');
axis manual;

i = 1;
while ishandle(fig)
    clf;
    hold on;
    
    % 绘制天花板（y = 0 的水平粗线）
    plot([-1.5, 1.5], [0, 0], 'k-', 'LineWidth', 3);
    
    % 计算各质量位置（x 坐标固定为 0）
    massPositions = zeros(n,2);
    for j = 1:n
        massPositions(j,:) = [0, sol(i, 2*j - 1)];
    end
    
    % 绘制从天花板到第 1 个质量的弹簧（正弦波形）
    sp = getSpringPoints(0, 0, massPositions(1,1), massPositions(1,2));
    plot(sp(1,:), sp(2,:), 'b-', 'LineWidth', 2);
    
    % 绘制相邻质量之间的弹簧
    for j = 1:(n-1)
        sp = getSpringPoints(massPositions(j,1), massPositions(j,2), ...
            massPositions(j+1,1), massPositions(j+1,2));
        plot(sp(1,:), sp(2,:), 'b-', 'LineWidth', 2);
    end
    
    % 绘制质量（红色圆点）
    for j = 1:n
        plot(massPositions(j,1), massPositions(j,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    end
    
    % 绘制施加在最下方质量上的额外力箭头
    arrowLength = (y_eq(end)/100) * F_ext;  % 按比例计算箭头长度
    arrowStart = massPositions(end,:);
    quiver(arrowStart(1), arrowStart(2), 0, arrowLength, 0, 'LineWidth', 2, 'MaxHeadSize', 0.01);
    
    % 设置坐标轴：固定 x 轴范围，y 轴从 0 到足够显示最低质量
    xlim([-1.5 1.5]);
    ylim([0, max(sol(:, 2*n - 1)) + 20]);
    set(gca, 'YDir', 'reverse');  % 使 y=0 显示在顶部
    title(sprintf('Time = %.2f s', t(i)));
    ylabel('Distance from Ceiling (m)');
    grid on;
    
    drawnow;
    pause(pauseTime);
    
    i = i + 1;
    if i > length(t)
        i = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 局部函数：计算初始状态（静力平衡位置）
function z0 = init_state(n, L0, g, k)
    % 对于第 i 个质量，其静力平衡位置：
    %   y_eq(i) = i*L0 + (g/k) * [ i*(n+1) - (i*(i+1))/2 ]
    y_eq = zeros(n, 1);
    for i = 1:n
        y_eq(i) = i * L0 + (g/k) * ( i*(n+1) - (i*(i+1))/2 );
    end
    % 初始状态：位置取静力平衡位置，速度均为 0
    z0 = zeros(2*n, 1);
    for i = 1:n
        z0(2*i - 1) = y_eq(i);
        z0(2*i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 局部函数：ODE 模型
function dzdt = mass_spring_ode(~, z, params)
    % 状态向量 z = [y1; v1; y2; v2; ...; yn; vn]
    n = params.n;
    k = params.k;
    L0 = params.L0;
    g = params.g;
    F_ext = params.F_ext;
    
    dzdt = zeros(2*n, 1);
    
    % 质量 1（连接天花板和质量 2）
    dzdt(1) = z(2);
    dzdt(2) = - k * ((z(1) - 0) - L0) + k * ((z(3) - z(1)) - L0) + g;
    
    % 质量 2 至 n-1
    for i = 2:(n-1)
        dzdt(2*i - 1) = z(2*i);
        dzdt(2*i) = - k * ((z(2*i - 1) - z(2*i - 3)) - L0) ...
            + k * ((z(2*i + 1) - z(2*i - 1)) - L0) + g;
    end
    
    % 质量 n（最下方，受额外 F_ext）
    dzdt(2*n - 1) = z(2*n);
    dzdt(2*n) = - k * ((z(2*n - 1) - z(2*n - 3)) - L0) + g + F_ext;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 局部函数：生成正弦波形弹簧轨迹
function sp = getSpringPoints(x1, y1, x2, y2)
    % 生成连接 (x1,y1) 与 (x2,y2) 的弹簧轨迹，采用正弦函数模拟
    % 要求：
    %   - 弹簧在竖直方向上从 y1 到 y2 绘制
    %   - 绘制固定个周期（本例设定 nCycles = 3，可根据需要调整）
    %   - 振幅固定为 amplitude
    %
    % 返回：
    %   sp：2×N 矩阵，第一行为 x 坐标，第二行为 y 坐标

    nCycles = 3;      % 绘制 3 个完整周期
    amplitude = 0.06;   % 固定振幅
    N = 100;           % 采样点数
    t = linspace(0, 1, N);
    
    % 竖直方向线性插值
    y = y1 + t * (y2 - y1);
    % 水平方向正弦波形
    x = amplitude * sin(2*pi*nCycles*t);
    
    sp = [x; y];
end
