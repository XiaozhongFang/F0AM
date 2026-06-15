%% F0AM ODE 组装与求解 — 逐步可视化教程
% 本脚本用 5 个物种、6 个反应的微型机理，展示 F0AM 核心的
% ODE 组装 → 稀疏矩阵 → 雅可比矩阵 → ode15s 求解 全过程。
%
% 可以在 MATLAB 中用"实时编辑器 (Live Editor)"打开本文件:
%   编辑器选项卡 → 另存为 → 保存类型: Live Code Files (.mlx)
% 或者直接在当前编辑器中使用"运行节"功能逐段执行。

clear; close all; clc;
fprintf('===== F0AM ODE 求解器可视化教程 =====\n\n')

%% 步骤 1：定义微型化学机理

% 物种列表 (nSp = 5)
Cnames = {'ONE'    % 占位物种，浓度恒为 1 (用于一级反应)
           'RO2'   % 有机过氧自由基总量 (占位)
           'A'     % 反应物 (类似 VOC)
           'B'     % 自由基 (类似 OH)
           'C'};   % 稳定产物

nSp = length(Cnames);
fprintf('物种 (%d 个):\n', nSp);
for i = 1:nSp
    fprintf('  %2d: %s\n', i, Cnames{i});
end

% 反应列表 (nRx = 6)
Rnames = {'A + B  → C + B'   % ① 催化氧化: A 被 B 氧化为 C，B 再生
           'B      → (loss)'  % ② 自由基终止
           'A      → (loss)'  % ③ A 光解损失
           'B + B  → (loss)'  % ④ B 自反应
           '→ A'              % ⑤ A 的排放 (零级)
           'C      → (loss)'};% ⑥ C 沉降损失

nRx = length(Rnames);
fprintf('\n反应 (%d 个):\n', nRx);
for i = 1:nRx
    fprintf('  %2d: %s\n', i, Rnames{i});
end

% 速率常数 (1/s 或 1/(molec·cm³·s) 取决于反应级数)
% 这里简化，直接给数值
k_values = [1e-3   % k₁: A + B → C + B
            0.01   % k₂: B → loss
            1e-4   % k₃: A → loss
            0.1    % k₄: B + B → loss
            0.5    % k₅: → A  (排放速率)
            1e-4]; % k₆: C → loss

fprintf('\n速率常数:\n');
for i = 1:nRx
    fprintf('  k%d = %.4g\n', i, k_values(i));
end

%% 步骤 2：构建核心数据结构 — iG 和稀疏矩阵 f

% 这是 F0AM 最关键的一步：将化学反应翻译成数值矩阵。
%
% iG (nRx × 3): 每个反应的 1~3 个反应物在 Cnames 中的索引
%   G(irxn) = conc(iG(irxn,1)) × conc(iG(irxn,2)) × conc(iG(irxn,3))
%   对于一级反应，用 ONE 填充多余的索引 (conc(ONE)=1)
%
% f (nRx × nSp): 稀疏化学计量系数矩阵
%   f(irxn, isp) = 反应 irxn 中物种 isp 的净化学计量系数
%   (负 = 消耗，正 = 生成)

% ---- 构建 iG ----
% 思路：每种反应需要 1~3 个反应物的乘积
% ONE 的索引 = 1，浓度恒为 1

iONE = 1;  % ONE 总是列表第一个物种

iG = zeros(nRx, 3);

% 反应 1: A + B → C + B   (二级反应: [A][B])
iG(1,:) = [3, 4, iONE];   % = [A, B, ONE] → G = [A]·[B]·1

% 反应 2: B → loss        (一级反应: [B])
iG(2,:) = [4, iONE, iONE]; % = [B, ONE, ONE] → G = [B]·1·1

% 反应 3: A → loss        (一级反应: [A])
iG(3,:) = [3, iONE, iONE]; % = [A, ONE, ONE]

% 反应 4: B + B → loss    (二级反应: [B]²)
iG(4,:) = [4, 4, iONE];    % = [B, B, ONE] → G = [B]·[B]·1

% 反应 5: → A             (零级反应: 1)
iG(5,:) = [iONE, iONE, iONE]; % G = 1·1·1 = 1

% 反应 6: C → loss        (一级反应: [C])
iG(6,:) = [5, iONE, iONE]; % = [C, ONE, ONE]

fprintf('\niG 矩阵 (nRx × 3): 每个反应的 1~3 个反应物索引\n');
fprintf('  %-20s | %-8s %-8s %-8s\n', '反应', 'idx₁', 'idx₂', 'idx₃');
for i = 1:nRx
    fprintf('  %-20s | %-8d %-8d %-8d\n', Rnames{i}, iG(i,1), iG(i,2), iG(i,3));
end

% ---- 构建稀疏矩阵 f ----
% 用 full 矩阵先演示，再转为 sparse

f_full = zeros(nRx, nSp);

% 反应 1: A + B → C + B  (A:-1, C:+1, B: net 0)
f_full(1, 3) = -1;  % A consumed
f_full(1, 5) = +1;  % C produced
% B: net 0 (consumed 1, produced 1)

% 反应 2: B → loss
f_full(2, 4) = -1;

% 反应 3: A → loss
f_full(3, 3) = -1;

% 反应 4: B + B → loss
f_full(4, 4) = -2;  % 2 B consumed

% 反应 5: → A
f_full(5, 3) = +1;  % A produced

% 反应 6: C → loss
f_full(6, 5) = -1;

% 转换为稀疏矩阵 (F0AM 实际使用的格式)
f_sparse = sparse(f_full);

fprintf('\n\n化学计量矩阵 f (nRx × nSp):\n');
fprintf('  行=反应, 列=物种 | ONE  RO2   A    B    C\n');
fprintf('%s\n', repmat('-', 1, 44));
for i = 1:nRx
    fprintf('  %-20s |', Rnames{i});
    for j = 1:nSp
        if f_full(i,j) ~= 0
            fprintf(' %+3d ', f_full(i,j));
        else
            fprintf('  .  ');
        end
    end
    fprintf('\n');
end

% ---- 可视化 f 矩阵 ----
figure('Position', [100, 100, 700, 400]);

subplot(1,2,1);
imagesc(f_full);
colormap(redbluecm);
colorbar;
title('f 矩阵 (full)');
xlabel('物种索引'); ylabel('反应索引');
set(gca, 'XTick', 1:nSp, 'XTickLabel', Cnames, ...
    'YTick', 1:nRx, 'YTickLabel', Rnames, 'YTickLabelRotation', 0);
xtickangle(45);

subplot(1,2,2);
spy(f_sparse, 'b', 15);
title('f 矩阵 (sparse 可视化)');
xlabel('物种列'); ylabel('反应行');
set(gca, 'XTick', 1:nSp, 'XTickLabel', Cnames, ...
    'YTick', 1:nRx, 'YTickLabel', Rnames);
xtickangle(45);
sgtitle('化学计量系数矩阵', 'FontSize', 14);

%% 步骤 3：dydt 求导 — 核心计算

% 给定当前浓度向量 conc (molec/cm³)，计算 dydt
% 这是 ODE 求解器每步都调用的核心函数
%
% 三个步骤:
%   ① G  = conc(iG₁) · conc(iG₂) · conc(iG₃)  (每个反应的"浓度积")
%   ② rates = k · G                             (每个反应的反应速率)
%   ③ dydt = rates · f                          (所有反应对每个物种的净变化之和)
%
% ★ f 是 sparse → 矩阵乘法极快，即使 20000 反应 × 7000 物种

fprintf('\n\n===== dydt 求导计算演示 =====\n');

% 假设当前时刻的浓度 (ppb)
conc_ppb = [1;      % ONE: 恒为 1
            0;      % RO2: 暂为 0
            100;    % A: 100 ppb
            0.05;   % B: 0.05 ppb
            0];     % C: 0 ppb

fprintf('\n当前浓度 (ppb):\n');
for i = 1:nSp
    fprintf('  %-6s = %.4g ppb\n', Cnames{i}, conc_ppb(i));
end

% F0AM 内部使用 molec/cm³，但这里为了简单直接使用 ppb 展示计算过程
% 因为 G 的乘积性质在两种单位下原理相同

% 步骤 ①: 计算 G (每个反应的浓度积)
G = zeros(nRx, 1);
for i = 1:nRx
    G(i) = conc_ppb(iG(i,1)) * conc_ppb(iG(i,2)) * conc_ppb(iG(i,3));
end

fprintf('\n步骤①: G(irxn) = conc(iG₁) × conc(iG₂) × conc(iG₃)\n');
for i = 1:nRx
    fprintf('  G%d (%s) = %s × %s × %s = %.4g\n', ...
        i, Rnames{i}, ...
        Cnames{iG(i,1)}, Cnames{iG(i,2)}, Cnames{iG(i,3)}, G(i));
end

% 步骤 ②: 计算反应速率
rates = k_values(:) .* G;

fprintf('\n步骤②: rates = k × G\n');
for i = 1:nRx
    fprintf('  rate%d (%s) = %.4g × %.4g = %.4g\n', ...
        i, Rnames{i}, k_values(i), G(i), rates(i));
end

% 步骤 ③: dydt = rates · f  (矩阵乘法)
dydt = rates' * f_full;  % 1 × nSp

fprintf('\n步骤③: dydt = rates · f   (关键一步!)\n');
fprintf('  对每个物种: dydt(isp) = Σ(rate(irxn) × f(irxn, isp))\n\n');
for j = 1:nSp
    contrib = '';
    total = 0;
    for i = 1:nRx
        if f_full(i,j) ~= 0
            val = rates(i) * f_full(i,j);
            total = total + val;
            if val > 0
                contrib = [contrib sprintf(' + %.4g × (%+d)', rates(i), f_full(i,j))];
            else
                contrib = [contrib sprintf(' - %.4g × (%+d)', abs(rates(i)), f_full(i,j))];
            end
        end
    end
    fprintf('  d[%s]/dt = %s = %.4g ppb/s\n', Cnames{j}, contrib, total);
end

%% 步骤 4：可视化 Jacobian 矩阵

% F0AM 提供解析 Jacobian → 求解速度快 10~100 倍
%
% 对于反应 k·A·B:
%   ∂(k·A·B)/∂A = k·B
%   ∂(k·A·B)/∂B = k·A
%
% 矩阵形式: J = f' · (∂rates/∂conc)
%   其中 ∂rates/∂conc 是稀疏的 nRx × nSp 矩阵
%
% 最终 J 是 nSp × nSp 的方阵

fprintf('\n\n===== Jacobian 矩阵 (解析) =====\n');

% 计算 ∂rates/∂conc (每个反应对每个物种的偏导)
DrDy = zeros(nRx, nSp);

for i = 1:nRx
    % 对第 1 个反应物求偏导:  ∂(k·c₁·c₂·c₃)/∂c₁ = k·c₂·c₃
    DrDy(i, iG(i,1)) = k_values(i) * conc_ppb(iG(i,2)) * conc_ppb(iG(i,3));
    
    % 对第 2 个反应物求偏导:  ∂(k·c₁·c₂·c₃)/∂c₂ = k·c₁·c₃
    DrDy(i, iG(i,2)) = DrDy(i, iG(i,2)) + k_values(i) * conc_ppb(iG(i,1)) * conc_ppb(iG(i,3));
    
    % 对第 3 个反应物求偏导:  ∂(k·c₁·c₂·c₃)/∂c₃ = k·c₁·c₂
    DrDy(i, iG(i,3)) = DrDy(i, iG(i,3)) + k_values(i) * conc_ppb(iG(i,1)) * conc_ppb(iG(i,2));
    
    % 处理自反应 (B+B): ∂(k·B²)/∂B = 2k·B
    if iG(i,1) == iG(i,2) && iG(i,1) ~= iONE
        DrDy(i, iG(i,1)) = 2 * k_values(i) * conc_ppb(iG(i,1));
    end
end

fprintf('\n∂rates/∂conc 矩阵 (nRx × nSp):\n');
fprintf('  行=反应, 列=物种 | ONE    RO2     A       B       C\n');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:nRx
    fprintf('  %-20s |', Rnames{i});
    for j = 1:nSp
        if abs(DrDy(i,j)) > 1e-10
            fprintf(' %6.3f', DrDy(i,j));
        else
            fprintf('    .  ');
        end
    end
    fprintf('\n');
end

% 最终 Jacobian: J = f' · DrDy
J_full = f_full' * DrDy;  % nSp × nSp

fprintf('\nJacobian J = f^T · DrDy (nSp × nSp):\n');
fprintf('  列=对谁的偏导, 行=哪个物种的方程\n');
fprintf('  %-6s', 'J(i,j)');
for j = 1:nSp
    fprintf(' | ∂/d%s    ', Cnames{j});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 6 + 13*nSp));
for i = 1:nSp
    fprintf('  d[%s]/dt', Cnames{i});
    for j = 1:nSp
        if abs(J_full(i,j)) > 1e-10
            fprintf(' | %10.4f', J_full(i,j));
        else
            fprintf(' |    .     ');
        end
    end
    fprintf('\n');
end

% ---- 可视化 Jacobian ----
figure('Position', [150, 150, 600, 500]);

imagesc(J_full);
colorbar;
colormap(redbluecm);
title(sprintf('Jacobean 矩阵 J (nSp×nSp)\n当前 [B]=%.3f ppb', conc_ppb(4)));
xlabel('偏导对象 (∂/∂cⱼ)');
ylabel('物种方程 (d[cᵢ]/dt)');
set(gca, 'XTick', 1:nSp, 'XTickLabel', Cnames, ...
    'YTick', 1:nSp, 'YTickLabel', Cnames);
xtickangle(45);

% 标注对角线 (自反馈)
hold on;
for i = 1:nSp
    if J_full(i,i) < 0
        text(i, i, '← 负反馈(稳定)', 'HorizontalAlignment', 'right', ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10);
    elseif J_full(i,i) > 0
        text(i, i, '→ 正反馈(不稳定)', 'HorizontalAlignment', 'left', ...
            'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10);
    end
end

%% 步骤 5：调用 ode15s 求解 ODE 系统

% 建立 ODE 函数 (与 F0AM 的 dydt_eval 完全一致)
fprintf('\n\n===== ODE 求解 =====\n');

% 定义 ODE 方程函数
dydt_func = @(t, y) ode_system(t, y, k_values, iG, f_full, nSp);

% 初始条件 (ppb)
y0 = [1; 0; 100; 0.05; 0];  % ONE, RO2, A, B, C

% 时间范围 (1 小时 = 3600 秒)
tspan = [0, 3600];

fprintf('初始条件: [A]=%.1f ppb, [B]=%.3f ppb, [C]=%.1f ppb\n', y0(3), y0(4), y0(5));
fprintf('模拟时长: %.0f 秒 (%.1f 小时)\n', tspan(2), tspan(2)/3600);

% 设置 ODE 选项: 提供 Jacobian (加速!)
fprintf('\n设置 odeset:\n');
fprintf('  ✔ Jacobian = 解析 Jacobian 函数\n');
fprintf('  ✔ 使用 ode15s (刚性求解器)\n');

J_func = @(t, y) jacobian_system(t, y, k_values, iG, f_full, nSp);

options = odeset('Jacobian', J_func);

% 求解
fprintf('\n正在调用 ode15s 求解...\n');
tic;
[t, y] = ode15s(dydt_func, tspan, y0, options);
elapsed = toc;
fprintf('完成! 用时 %.3f 秒, %d 个时间步\n', elapsed, length(t));

% ---- 绘图: 浓度时间序列 ----
figure('Position', [200, 200, 800, 500]);

subplot(2,2,1);
plot(t, y(:,3), 'b-', 'LineWidth', 2); hold on;
plot(t, y(:,5), 'r-', 'LineWidth', 2);
plot(t, y(:,4), 'g-', 'LineWidth', 2);
xlabel('时间 (秒)'); ylabel('浓度 (ppb)');
legend('A (反应物)', 'C (产物)', 'B (自由基)', 'Location', 'best');
title('浓度随时间变化');
grid on; box on;

subplot(2,2,2);
semilogy(t, y(:,4), 'g-', 'LineWidth', 2);
xlabel('时间 (秒)'); ylabel('[B] (ppb, 对数)');
title('B (自由基) 浓度 — 对数尺度');
grid on; box on;

% 反应速率随时间变化
subplot(2,2,3);
G_all = zeros(length(t), nRx);
for j = 1:nRx
    G_all(:,j) = y(:, iG(j,1)) .* y(:, iG(j,2)) .* y(:, iG(j,3));
end
rates_all = k_values(:)' .* G_all;
plot(t, rates_all, 'LineWidth', 2);
xlabel('时间 (秒)'); ylabel('反应速率 (ppb/s)');
legend(Rnames, 'Location', 'best', 'FontSize', 8);
title('各反应速率变化');
grid on; box on;

% 浓度和 (质量守恒检查)
subplot(2,2,4);
conc_sum = sum(y(:,3:5), 2);
plot(t, conc_sum, 'k-', 'LineWidth', 2);
xlabel('时间 (秒)'); ylabel('A+B+C 总和 (ppb)');
title('质量守恒检查 (A+B+C)');
grid on; box on;
ylim([min(conc_sum)*0.95, max(conc_sum)*1.05]);

sgtitle('ode15s 求解结果 — 微型化学机理', 'FontSize', 14);

%% 步骤 6：Jacobian 随时间演化 — 热力图动画

fprintf('\n\n===== Jacobian 热力图动画 =====\n');

% 选取几个时间点展示 Jacobian 如何随浓度变化
n_frames = 6;
idx_frames = round(linspace(1, length(t), n_frames));
time_points = t(idx_frames);

figure('Position', [250, 250, 900, 600]);

for k = 1:n_frames
    yy = y(idx_frames(k), :)';
    J_now = f_full' * drdy_at_conc(yy, k_values, iG, nSp);
    
    subplot(2, 3, k);
    imagesc(J_now, [-max(abs(J_now(:)))*0.8, max(abs(J_now(:)))*0.8]);
    colormap(redbluecm);
    colorbar;
    title(sprintf('t = %.0f s\n[B] = %.4f ppb', t(idx_frames(k)), yy(4)));
    xlabel('∂/∂cⱼ'); ylabel('dcᵢ/dt');
    set(gca, 'XTick', 1:nSp, 'XTickLabel', Cnames, ...
        'YTick', 1:nSp, 'YTickLabel', Cnames);
    xtickangle(45);
end
sgtitle('Jacobian 矩阵随时间的演化', 'FontSize', 14);

%% 步骤 7：完整 F0AM 风格求解 (含 RO2)

% 在 F0AM 中，RO2 = 所有有机过氧自由基之和
% 这里模拟 RO2 由 A 氧化产生
%
% 真正 F0AM: RO2 在每步 dydt 中更新:
%   conc(:,2) = sum(conc(:, iRO2_list), 2);
%
% 这里我们展示完整流程: A+B → B+C 生成的 C 的一部分是 RO2
% 简化: 设 RO2 = C × 0.3 (30% 产物是自由基)

fprintf('\n\n===== 完整 F0AM 流程模拟 =====\n');

% 重新定义 ODE 函数: 包含 RO2 更新
dydt_f0am = @(t, y) ode_system_f0am(t, y, k_values, iG, f_full, nSp, 3);

% 初始条件 (包含 RO2)
y0_f0am = [1; 0.01; 100; 0.05; 0];  % RO2 初始 0.01 ppb

% 求解
options2 = odeset('Jacobian', @(t,y) jacobian_system_f0am(t,y,k_values,iG,f_full,nSp,3));
[t2, y2] = ode15s(dydt_f0am, [0, 3600], y0_f0am, options2);

% 绘图: 对比有无 RO2 的差异
figure('Position', [300, 300, 900, 400]);

subplot(1,3,1);
plot(t2, y2(:,3), 'b-', 'LineWidth', 2); hold on;
plot(t2, y2(:,5), 'r-', 'LineWidth', 2);
plot(t2, y2(:,2), 'm-', 'LineWidth', 2);
plot(t2, y2(:,4), 'g-', 'LineWidth', 2);
xlabel('时间 (秒)'); ylabel('浓度 (ppb)');
legend('A', 'C', 'RO2', 'B', 'Location', 'best');
title('含 RO₂ 的完整模拟');
grid on; box on;

subplot(1,3,2);
% dydt 各组分可视化 (在 t=100s 时)
[~, idx_t] = min(abs(t2-100));
y_sample = y2(idx_t,:)';
dydt_sample = dydt_f0am(t2(idx_t), y_sample);
bar(dydt_sample);
set(gca, 'XTick', 1:nSp, 'XTickLabel', Cnames);
title(sprintf('dydt @ t=100s (ppb/s)'));
ylabel('dC/dt (ppb/s)');
grid on; box on;

subplot(1,3,3);
spy(f_full' * drdy_at_conc(y_sample, k_values, iG, nSp), 'b', 12);
title(sprintf('Jacobian 稀疏模式 @ t=100s'));
xlabel('物种'); ylabel('物种');

sgtitle('F0AM 风格完整求解', 'FontSize', 14);

%% 总结

fprintf('\n\n===== 总结 =====\n');
fprintf('F0AM ODE 求解的核心步骤:\n\n');
fprintf('  ① 化学机制 → 稀疏矩阵 f (nRx × nSp)\n');
fprintf('  ② iG 记录每个反应的 1~3 个反应物索引\n');
fprintf('  ③ 每步 ODE 求导:\n');
fprintf('       G = conc(iG₁) .* conc(iG₂) .* conc(iG₃)\n');
fprintf('       rates = k .* G\n');
fprintf('       dydt = rates * f    ← 一行代码!\n');
fprintf('  ④ 解析 Jacobian: J = f^T · DrDy\n');
fprintf('  ⑤ ode15s 利用 J 高速求解\n');
fprintf('  ⑥ 后处理: ppb 转换 + 速率分析\n');
fprintf('\n关键: 稀疏矩阵让 20000 反应 × 7000 物种\n');
fprintf('      的矩阵乘法在毫秒级完成!\n');

%% ========== 辅助函数 ==========

function dydt = ode_system(~, y, k, iG, f, nSp)
    % ODE 导数函数 (与 F0AM dydt_eval 原理一致)
    y = y(:);
    
    % 计算反应物浓度积
    G = y(iG(:,1)) .* y(iG(:,2)) .* y(iG(:,3));
    
    % 反应速率
    rates = k(:) .* G;
    
    % 核心: dydt = rates · f
    dydt = rates' * f;
    dydt = dydt(:);
    
    % ONE 不动, RO2 暂为 0
    dydt(1) = 0;
    dydt(2) = 0;
end

function J = jacobian_system(~, y, k, iG, f, nSp)
    % 解析 Jacobian (与 F0AM Jac_eval 原理一致)
    y = y(:);
    
    % 计算 ∂rates/∂conc
    DrDy = drdy_at_conc(y, k, iG, nSp);
    
    % J = f^T · DrDy
    J = f' * DrDy;
    
    % ONE 和 RO2 的行置零
    J(1,:) = 0;
    J(2,:) = 0;
end

function dydt = ode_system_f0am(~, y, k, iG, f, nSp, iRO2_species)
    % 含 RO2 更新的 ODE 导数函数
    y = y(:);
    
    % 更新 RO2 = 指定物种的浓度
    y(2) = y(iRO2_species) * 0.3;  % 30% 产物作为 RO2
    
    % 计算反应物浓度积
    G = y(iG(:,1)) .* y(iG(:,2)) .* y(iG(:,3));
    
    % 反应速率
    rates = k(:) .* G;
    
    % 核心: dydt = rates · f
    dydt = rates' * f;
    dydt = dydt(:);
    
    dydt(1) = 0;  % ONE 不动
end

function J = jacobian_system_f0am(~, y, k, iG, f, nSp, iRO2_species)
    % 含 RO2 更新的 Jacobian
    y = y(:);
    y(2) = y(iRO2_species) * 0.3;
    
    DrDy = drdy_at_conc(y, k, iG, nSp);
    J = f' * DrDy;
    J(1,:) = 0;
    J(2,:) = 0;
end

function DrDy = drdy_at_conc(y, k, iG, nSp)
    % 计算 ∂rates/∂conc 在当前浓度 y 下的值
    nRx = length(k);
    DrDy = zeros(nRx, nSp);
    
    for i = 1:nRx
        % 对每个反应物的偏导
        DrDy(i, iG(i,1)) = k(i) * y(iG(i,2)) * y(iG(i,3));
        DrDy(i, iG(i,2)) = DrDy(i, iG(i,2)) + k(i) * y(iG(i,1)) * y(iG(i,3));
        DrDy(i, iG(i,3)) = DrDy(i, iG(i,3)) + k(i) * y(iG(i,1)) * y(iG(i,2));
        
        % 自反应修正
        if (iG(i,1) == iG(i,2) && iG(i,1) ~= 1)
            DrDy(i, iG(i,1)) = 2 * k(i) * y(iG(i,1));
        end
        if (iG(i,1) == iG(i,3) && iG(i,1) ~= 1)
            DrDy(i, iG(i,1)) = 2 * k(i) * y(iG(i,1));
        end
        if (iG(i,3) == iG(i,2) && iG(i,3) ~= 1)
            DrDy(i, iG(i,3)) = 2 * k(i) * y(iG(i,3));
        end
    end
end

function cmap = redbluecm
    % 红蓝分色图 (正=红, 负=蓝, 0=白)
    n = 64;
    cmap = zeros(n, 3);
    for i = 1:n
        t = (i-1)/(n-1);
        if t < 0.5
            % 蓝到白
            cmap(i,:) = [t*2, t*2, 1];
        else
            % 白到红
            cmap(i,:) = [1, (1-t)*2, (1-t)*2];
        end
    end
end
