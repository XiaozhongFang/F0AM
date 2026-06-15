%% 烟雾箱示例 — 机制解析与雅可比矩阵组装
% 本脚本基于真实的 MCMv3.3.1 机制 (610 物种, 1974 反应)，
% 逐步展示 F0AM 如何将化学反应翻译为数值矩阵，并构建解析雅可比。
%
% 前置条件: 运行本脚本前确保 F0AM 路径已添加:
%   addpath(genpath('D:\Code\F0AM\Core'));
%   addpath(genpath('D:\Code\F0AM\Chem'));

clear; close all; clc;
fprintf('===== F0AM Chamber 示例: 机制解析与雅可比组装 =====\n\n');

%% ============================
%% 第一步：设置气象与初始条件
%% ============================

% 与 ExampleSetup_Chamber.m 完全一致
Met = {...
    'P'         1013                       ; % mbar
    'T'         298                        ; % K
    'RH'        10                         ; % %
    'LFlux'     'ExampleLightFlux.txt'     ; % 光谱文件
    'jcorr'     1                          ; % 光强因子
    };

InitConc = {...
    'C5H8'      10                    0;   % 异戊二烯 10 ppb
    'NO2'       [.1; 1; 10]           0;   % 三组 NOx
    'H2O2'      200                   0;   % OH 源
    };

ChemFiles = {...
    'MCMv331_K(Met)';
    'MCMv331_J(Met,1)'; % Jmethod=1: BottomUp
    'MCMv331_Inorg_Isoprene';
    'CH3ONO_hv';
    };

BkgdConc = {...
    'DEFAULT'       0;
    };

ModelOptions.Verbose       = 0;  % 安静模式
ModelOptions.EndPointsOnly = 0;
ModelOptions.LinkSteps     = 0;
ModelOptions.IntTime       = 3*3600;
ModelOptions.SavePath      = 'DoNotSave';  % 不保存
ModelOptions.GoParallel    = 0;

%% ============================
%% 第二步：解析化学机制 — 核心
%% ============================

fprintf('正在解析化学机制...\n');

% 模拟 ModelCore 的内部初始化
% 先转换输入格式
InitConc_raw = InitConc;
[Met_str, InitConc_str, BkgdConc_str] = deal(Met, InitConc, BkgdConc);

holdFlag = logical(cell2mat(InitConc(:,3)));
Met      = breakout(Met(:,2),Met(:,1));
InitConc = breakout(InitConc(:,2),InitConc(:,1));
BkgdConc = breakout(BkgdConc(:,2),BkgdConc(:,1));

% 计算数密度
Met.M = NumberDensity(Met.P,Met.T);
Met.H2O = ConvertHumidity(Met.T,Met.P,Met.RH,'RH','NumberDensity');

% ===== 关键: 调用 InitializeChemistry =====
tic;
[Cnames, Rnames, k, f_sparse, iG, iRO2, jcorr, jcorr_all, iLR] = ...
    InitializeChemistry(Met, ChemFiles, ModelOptions, 1);
t_parse = toc;

nSp = length(Cnames);
nRx = length(Rnames);
nIc = length(Met.T);

fprintf('完成! 用时 %.2f 秒\n', t_parse);
fprintf('  物种数: %d\n', nSp);
fprintf('  反应数: %d\n', nRx);
fprintf('  步骤数: %d (NO2 的三个浓度)\n', nIc);

%% ---- 可视化: 机制全貌 ----
figure('Position', [50, 50, 1200, 800]);

% 稀疏矩阵 f 的全貌
subplot(2,3,1);
spy(f_sparse, 'k', 1);
title(sprintf('稀疏矩阵 f (%d×%d)\n非零元素: %d (%.4f%%)', ...
    nRx, nSp, nnz(f_sparse), 100*nnz(f_sparse)/(nRx*nSp)));
xlabel('物种'); ylabel('反应');

% 每行(反应)的非零元素数分布
subplot(2,3,2);
row_nz = sum(f_sparse ~= 0, 2);
histogram(row_nz, 0:max(row_nz), 'FaceColor', [0.3 0.5 0.8]);
title('每个反应涉及的物种数分布');
xlabel('非零系数个数/反应'); ylabel('反应数');
grid on;

% 每列(物种)的非零元素数分布
subplot(2,3,3);
col_nz = sum(f_sparse ~= 0, 1);
histogram(col_nz, 0:2:max(col_nz), 'FaceColor', [0.8 0.4 0.3]);
title('每个物种参与的反应数分布');
xlabel('参与反应数/物种'); ylabel('物种数');
grid on;

% 反应级数分布
subplot(2,3,4);
rxn_order = sum(iG ~= 1, 2); % 非 ONE 的反应物数 = 反应级数
rxn_order(rxn_order > 2) = 2; % 三级反应归为 2+
counts = histcounts(rxn_order, [0 1 2 3]);
pie(counts, {'一级(光解/热解)', '二级(碰撞)', '三级+'});
title(sprintf('反应级数分布 (共 %d 个反应)', nRx));

% 反应类型: 光解 vs 热反应
subplot(2,3,5);
is_jval = contains(Rnames, '+ hv');
pie([sum(is_jval), sum(~is_jval)], ...
    {sprintf('光解反应 (%d)', sum(is_jval)), ...
     sprintf('热反应 (%d)', sum(~is_jval))});
title('反应类型分布');

% k 值分布
subplot(2,3,6);
k_vals = k(1,:); % 第一组步的 k 值
k_plot = k_vals(k_vals > 0 & k_vals < 1e6); % 去掉极端值
histogram(log10(k_plot), 50, 'FaceColor', [0.2 0.7 0.4]);
title('速率常数 k 分布 (log₁₀)');
xlabel('log₁₀(k)'); ylabel('反应数');
grid on;

sgtitle(sprintf('MCMv3.3.1 机制全景 (%d 物种 × %d 反应)', nSp, nRx), ...
    'FontSize', 16);

%% ---- 局部放大: 无机组分 ----

% 找出无机物种索引 (常见无机物)
inorg_list = {'O','O1D','O3','OH','HO2','H2O2','NO','NO2','NO3',...
    'N2O5','HONO','HNO3','HO2NO2','CO','H2','SO2','SO3','HSO3','NA','SA'};
[~, inorg_idx] = ismember(inorg_list, Cnames);
inorg_idx(inorg_idx == 0) = [];
inorg_names = Cnames(inorg_idx);

% 找出仅涉及无机物的反应
rxn_inorg = all(ismember(iG, [1; inorg_idx(:)]), 2);
rxn_inorg_idx = find(rxn_inorg);
fprintf('\n无机组分: %d 个物种\n', length(inorg_idx));
fprintf('无机反应: %d 个\n', sum(rxn_inorg));

figure('Position', [100, 100, 800, 600]);

% 无机子矩阵
subplot(2,2,1);
f_inorg = f_sparse(rxn_inorg, inorg_idx);
spy(f_inorg, 'b', 8);
title(sprintf('无机组分 f 子矩阵 (%d×%d)', size(f_inorg,1), size(f_inorg,2)));
xlabel('物种'); ylabel('反应');
set(gca, 'XTick', 1:length(inorg_names), 'XTickLabel', inorg_names, 'FontSize', 7);
xtickangle(90);

% 无机反应举例
subplot(2,2,2);
example_n = min(20, sum(rxn_inorg));
example_idx = rxn_inorg_idx(1:example_n);
T_example = table(Rnames(example_idx), ...
    'VariableNames', {'反应式'});
disp(T_example);
% 在图中用 text 显示
cla;
text(0, 0.5, sprintf('无机反应示例 (前 %d 个):\n\n', example_n), ...
    'FontSize', 12, 'FontWeight', 'bold');
for i = 1:example_n
    text(0, 0.5 - i*0.045, sprintf('  %s', Rnames{example_idx(i)}), ...
        'FontSize', 9, 'FontName', 'FixedWidth');
end
axis off;
title('前 20 个无机反应');

% iG 中无机部分的反应物分布
subplot(2,2,3);
iG_inorg = iG(rxn_inorg, :);
reactant1 = iG_inorg(:,1); reactant1(reactant1==1) = [];
reactant2 = iG_inorg(:,2); reactant2(reactant2==1) = [];
all_reactants = [reactant1; reactant2];
[unique_r, ~, idx_r] = unique(all_reactants);
counts_r = accumarray(idx_r, 1);
[sorted_counts, sort_i] = sort(counts_r, 'descend');
top_n = min(10, length(sorted_counts));
bar(categorical(Cnames(unique_r(sort_i(1:top_n)))), sorted_counts(1:top_n));
title('无机反应中最常出现的反应物');
ylabel('出现次数');
grid on;

% 无机反应的 k 值
subplot(2,2,4);
k_inorg = k(1, rxn_inorg);
k_inorg_plot = k_inorg(k_inorg > 0);
histogram(log10(k_inorg_plot), 20, 'FaceColor', [0.2 0.6 0.3]);
title('无机反应 k 值分布');
xlabel('log₁₀(k)'); ylabel('反应数');
grid on;

sgtitle('无机子机制详解', 'FontSize', 14);

%% ============================
%% 第三步: 初始化浓度并转换单位
%% ============================

fprintf('\n===== 初始化浓度 =====\n');

% 建立初始浓度向量 (ppb → molec/cm³)
conc_init = zeros(1, nSp);
[isSp, iC] = ismember(fieldnames(InitConc), Cnames);
fn = fieldnames(InitConc);
for i = 1:length(iC)
    if isSp(i)
        conc_init(iC(i)) = InitConc.(fn{i})(1); % 第一组 NOx
    end
end

% ppb → molec/cm³
conc_init = conc_init .* Met.M(1) ./ 1e9;
conc_init(1) = 1; % ONE = 1

% 输出初始浓度 (前 20 个非零)
init_ppb = conc_init ./ Met.M(1) .* 1e9;
non_zero = find(init_ppb > 1e-6);
fprintf('初始非零浓度物种 (%d 个):\n', length(non_zero));
for i = 1:min(10, length(non_zero))
    idx = non_zero(i);
    if idx <= length(Cnames)
        fprintf('  %s = %.2f ppb\n', Cnames{idx}, init_ppb(idx));
    end
end

%% ============================
%% 第四步: 组装 dydt 和 Jacobian
%% ============================

fprintf('\n===== dydt 与 Jacobian 组装 =====\n');

% 计算 dydt
[G_vals, rates, dydt] = compute_dydt(conc_init, k(1,:)', iG, f_sparse, iRO2, nSp);

% 找出最重要的反应 (速率最大的前 10)
[rates_sorted, rates_idx] = sort(rates, 'descend');
top_rates = rates_idx(1:min(10, length(rates_idx)));

fprintf('\n速率最大的前 10 个反应:\n');
fprintf('  %-50s %12s %12s\n', '反应', 'k', '速率(ppb/s)');
for i = 1:length(top_rates)
    ri = top_rates(i);
    fprintf('  %-50s %10.4g  %10.4g\n', Rnames{ri}, k(1,ri), rates_sorted(i));
end

% 找出 dydt 最大的物种
[dydt_sorted, dydt_idx] = sort(abs(dydt), 'descend');
top_dydt = dydt_idx(1:min(10, length(dydt_idx)));
fprintf('\nd|C|/dt 最大的前 10 个物种:\n');
for i = 1:length(top_dydt)
    si = top_dydt(i);
    fprintf('  %-20s dC/dt = %+10.4g ppb/s (C = %.4g ppb)\n', ...
        Cnames{si}, dydt(si), init_ppb(si));
end

%% ---- 可视化重要反应和通量 ----

figure('Position', [150, 150, 900, 500]);

subplot(1,3,1);
barh(dydt_sorted(min(10,length(dydt_sorted)):-1:1), 'FaceColor', [0.3 0.6 0.8]);
set(gca, 'YTickLabel', Cnames(top_dydt(min(10,length(top_dydt)):-1:1)));
xlabel('dC/dt (ppb/s)');
title('前 10 大物种通量');
grid on;

subplot(1,3,2);
% 参与 C5H8 反应的通量分析
iC5H8 = find(strcmp(Cnames, 'C5H8'));
if ~isempty(iC5H8)
    c5_rxns = find(f_sparse(:, iC5H8) ~= 0);
    c5_rates = rates(c5_rxns);
    [c5_sorted, c5_si] = sort(abs(c5_rates), 'descend');
    top_c5 = min(10, length(c5_sorted));
    barh(c5_sorted(top_c5:-1:1), 'FaceColor', [0.8 0.4 0.3]);
    set(gca, 'YTickLabel', Rnames(c5_rxns(c5_si(top_c5:-1:1))), 'FontSize', 7);
    xlabel('反应速率 (ppb/s)');
    title('C₅H₈ 相关反应通量');
    grid on;
end

subplot(1,3,3);
% 稀疏 Jacobian 可视化
DrDy = compute_drdy(conc_init, k(1,:)', iG, nSp, iRO2);
J_sparse = f_sparse' * DrDy;
spy(J_sparse, 'b', 1);
title(sprintf('Jacobian 稀疏模式\n(%d×%d, %d 非零)', nSp, nSp, nnz(J_sparse)));
xlabel('物种'); ylabel('物种');

sgtitle('初始时刻的化学通量与 Jacobian', 'FontSize', 14);

%% ============================
%% 第五步: 逐步组装 Jacobian 的细节
%% ============================

fprintf('\n===== Jacobian 细节解析 =====\n');
fprintf('采用 F0AM 的解析方法 J = f^T · DrDy\n');
fprintf('DrDy: %d 反应 × %d 物种\n', nRx, nSp);

% 选取一个代表性区域: 围绕 C5H8 的 Jacobian 子矩阵
iC5H8 = find(strcmp(Cnames, 'C5H8'));
iOH   = find(strcmp(Cnames, 'OH'));
iO3   = find(strcmp(Cnames, 'O3'));
iNO3  = find(strcmp(Cnames, 'NO3'));

key_species = [iC5H8, iOH, iO3, iNO3];
key_names = Cnames(key_species);

if all(key_species > 0)
    % 提取子 Jacobian
    J_sub = J_sparse(key_species, key_species);
    J_sub_full = full(J_sub);
    
    figure('Position', [200, 200, 700, 500]);
    
    subplot(1,2,1);
    imagesc(J_sub_full);
    colorbar;
    colormap(redbluecm_custom);
    set(gca, 'XTick', 1:length(key_names), 'XTickLabel', key_names);
    set(gca, 'YTick', 1:length(key_names), 'YTickLabel', key_names);
    title(sprintf('关键物种 Jacobian 子矩阵\n(初始时刻, 单位: 1/s)'));
    % 标注数值
    for i = 1:length(key_species)
        for j = 1:length(key_species)
            text(j, i, sprintf('%.2g', J_sub_full(i,j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', ...
                'Color', 0.5*(1-sign(J_sub_full(i,j)))*[1 1 1]);
        end
    end
    
    subplot(1,2,2);
    % 追踪一个 Jacobian 元素的贡献来源
    % 例如 d(d[C5H8]/dt)/d[OH] = fᵀ · DrDy
    target_sp = iC5H8; % C5H8 的方程
    wrt_sp    = iOH;    % 对 OH 求偏导
    
    % 找到所有涉及 C5H8 和 OH 的反应
    both_rxns = find(f_sparse(:, target_sp) ~= 0 & ...
                     any(iG == wrt_sp, 2));
    
    if ~isempty(both_rxns)
        contribution = zeros(length(both_rxns), 1);
        for ri = 1:length(both_rxns)
            rxn = both_rxns(ri);
            % J_ij = Σ f_ir · DrDy_rj
            contribution(ri) = full(f_sparse(rxn, target_sp)) * DrDy(rxn, wrt_sp);
        end
        
        [contrib_sorted, ci] = sort(abs(contribution), 'descend');
        top_contrib = min(15, length(contrib_sorted));
        
        barh(contrib_sorted(top_contrib:-1:1), 'FaceColor', [0.3 0.5 0.8]);
        set(gca, 'YTickLabel', Rnames(both_rxns(ci(top_contrib:-1:1))), 'FontSize', 8);
        xlabel('贡献值');
        title(sprintf('J(%s, %s) 各反应贡献分解', ...
            Cnames{target_sp}, Cnames{wrt_sp}));
        grid on;
    end
    
    sgtitle('Jacobian 关键元素分解', 'FontSize', 14);
end

%% ============================
%% 第六步: 运行 ODE 并追踪变化
%% ============================

fprintf('\n===== ODE 求解 =====\n');

% 构建 ODE 函数句柄 (使用完整的 f_sparse 和 Jacobian)
ode_fun = @(t, y) dydt_full(t, y, k(1,:)', iG, f_sparse, iRO2, nSp);
ode_jac = @(t, y) jacobian_full(t, y, k(1,:)', iG, f_sparse, iRO2, nSp);

options = odeset('Jacobian', ode_jac, ...
    'RelTol', 1e-6, 'AbsTol', 1e-12);

fprintf('求解中 (3 小时 = 10800 秒)...\n');
tic;
[t, y] = ode15s(ode_fun, [0 10800], conc_init', options);
t_ode = toc;
fprintf('完成! %.1f 秒, %d 个时间步\n', t_ode, length(t));

% 转换回 ppb
y_ppb = y ./ Met.M(1) .* 1e9;

%% ---- 结果可视化 ----

fprintf('\n绘制结果...\n');

% 关键物种跟踪
key_track = {'C5H8','OH','HO2','NO','NO2','O3','HCHO','H2O2','MVK','MACR'};
[~, key_idx] = ismember(key_track, Cnames);
key_idx(key_idx == 0) = [];

figure('Position', [50, 50, 1200, 700]);

% 主反应物
subplot(2,3,1);
plot(t, y_ppb(:, key_idx(1)), 'b-', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('ppb');
title('C₅H₈ (异戊二烯)');
grid on;

% 自由基
subplot(2,3,2);
semilogy(t, max(y_ppb(:, key_idx(2)), 1e-10), 'r-', 'LineWidth', 1.5); hold on;
semilogy(t, max(y_ppb(:, key_idx(3)), 1e-10), 'g-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('ppb (log)');
legend('OH', 'HO₂', 'Location', 'best');
title('自由基');
grid on;

% NOx
subplot(2,3,3);
plot(t, y_ppb(:, key_idx(4)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_ppb(:, key_idx(5)), 'r-', 'LineWidth', 1.5);
plot(t, y_ppb(:, key_idx(4)) + y_ppb(:, key_idx(5)), 'k--', 'LineWidth', 1);
xlabel('时间 (s)'); ylabel('ppb');
legend('NO', 'NO₂', 'NOₓ', 'Location', 'best');
title('氮氧化物');
grid on;

% O3
subplot(2,3,4);
plot(t, y_ppb(:, key_idx(6)), 'Color', [0 0.6 0.8], 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('ppb');
title('O₃');
grid on;

% 产物
subplot(2,3,5);
plot(t, y_ppb(:, key_idx(7)), 'm-', 'LineWidth', 1.5); hold on;
plot(t, y_ppb(:, key_idx(8)), 'Color', [0.8 0.5 0.2], 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('ppb');
legend('HCHO', 'H₂O₂', 'Location', 'best');
title('产物');
grid on;

% MVK + MACR (异戊二烯氧化标志物)
subplot(2,3,6);
plot(t, y_ppb(:, key_idx(9)), 'LineWidth', 1.5); hold on;
plot(t, y_ppb(:, key_idx(10)), 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('ppb');
legend('MVK', 'MACR', 'Location', 'best');
title('异戊二烯氧化产物');
grid on;

sgtitle('Chamber 示例 — MCMv3.3.1 模拟结果', 'FontSize', 16);

%% ---- Jacobian 随时间的演化 ----

fprintf('\n生成 Jacobian 时间演化...\n');

% 选取 n 个时间点
n_snapshots = 12;
snap_idx = round(linspace(1, length(t), n_snapshots));

figure('Position', [50, 50, 1400, 800]);

% 创建热点图: 关键物种子 Jacobian 随时间变化
key_idx2 = [iC5H8, iOH, ...
    find(strcmp(Cnames, 'HO2')), ...
    find(strcmp(Cnames, 'NO')), ...
    find(strcmp(Cnames, 'NO2')), ...
    find(strcmp(Cnames, 'O3'))];
key_idx2(key_idx2 == 0) = [];
key_names2 = Cnames(key_idx2);

k_vals = k(1,:)';  % 保存速率常数

for snap = 1:n_snapshots
    tk = t(snap_idx(snap));
    yk = y(snap_idx(snap), :)';
    
    Dk = compute_drdy(yk, k_vals, iG, nSp, iRO2);
    Jk = f_sparse' * Dk;
    Jk_sub = full(Jk(key_idx2, key_idx2));
    
    subplot(3, 4, snap);
    imagesc(Jk_sub, [-max(abs(Jk_sub(:)))*0.5, max(abs(Jk_sub(:)))*0.5]);
    colormap(redbluecm_custom);
    if snap == 1 || snap == 5 || snap == 9
        set(gca, 'YTick', 1:length(key_names2), 'YTickLabel', key_names2, 'FontSize', 7);
    end
    if snap >= 9
        set(gca, 'XTick', 1:length(key_names2), 'XTickLabel', key_names2, 'FontSize', 7);
        xtickangle(45);
    end
    title(sprintf('t = %.0f s', tk), 'FontSize', 9);
end
sgtitle('关键物种 Jacobian 子矩阵随时间演化', 'FontSize', 16);
colorbar('Position', [0.92 0.1 0.02 0.8]);

%% ============================
%% 总结
%% ============================

fprintf('\n\n===== 总结 =====\n');
fprintf('MCMv3.3.1 机制解析过程:\n');
fprintf('  ① ChemFiles 中的 K/J 函数 → 计算 k 向量 (%d×1)\n', nRx);
fprintf('  ② SpeciesToAdd/ReactionToAdd → 构建 iG (%d×3) + f (%d×%d, 稀疏)\n', nRx, nRx, nSp);
fprintf('  ③ G = conc(iG₁)·conc(iG₂)·conc(iG₃) → 浓度积向量\n');
fprintf('  ④ rates = k·G → 反应速率\n');
fprintf('  ⑤ dydt = rates·f → 物种浓度变化率\n\n');
fprintf('Jacobian 组装:\n');
fprintf('  ⑥ DrDy (∂rates/∂conc) 解析计算:\n');
fprintf('     - 二级反应: ∂(k·A·B)/∂A = k·B, ∂(k·A·B)/∂B = k·A\n');
fprintf('     - 自反应: ∂(k·B²)/∂B = 2k·B\n');
fprintf('  ⑦ J = f^T · DrDy → %d×%d 稀疏 Jacobian (%d 非零)\n', nSp, nSp, nnz(J_sparse));
fprintf('  ⑧ ode15s 利用解析 Jacobian 加速求解\n');

%% ========== 辅助函数 ==========

function [G, rates, dydt] = compute_dydt(conc, k, iG, f, iRO2, nSp)
    conc = conc(:);
    conc(2) = sum(conc(iRO2)); % RO2 更新
    
    G = conc(iG(:,1)) .* conc(iG(:,2)) .* conc(iG(:,3));
    rates = k .* G;
    dydt = (rates' * f)';
    dydt(1) = 0; % ONE
    % 转换为 ppb/s 用于展示
    M_ref = 2.46e19; % 标准大气数密度
    dydt = dydt ./ M_ref .* 1e9;
end

function DrDy = compute_drdy(conc, k, iG, nSp, iRO2)
    nRx = length(k);
    conc = conc(:);
    conc(2) = sum(conc(iRO2));
    
    DrDy = zeros(nRx, nSp);
    for i = 1:nRx
        DrDy(i, iG(i,1)) = k(i) * conc(iG(i,2)) * conc(iG(i,3));
        DrDy(i, iG(i,2)) = DrDy(i, iG(i,2)) + k(i) * conc(iG(i,1)) * conc(iG(i,3));
        DrDy(i, iG(i,3)) = DrDy(i, iG(i,3)) + k(i) * conc(iG(i,1)) * conc(iG(i,2));
        
        % 自反应修正
        if (iG(i,1) == iG(i,2) && iG(i,1) ~= 1)
            DrDy(i, iG(i,1)) = 2 * k(i) * conc(iG(i,1));
        end
        if (iG(i,1) == iG(i,3) && iG(i,1) ~= 1)
            DrDy(i, iG(i,1)) = 2 * k(i) * conc(iG(i,1));
        end
        if (iG(i,3) == iG(i,2) && iG(i,3) ~= 1)
            DrDy(i, iG(i,3)) = 2 * k(i) * conc(iG(i,3));
        end
    end
end

function dydt = dydt_full(~, y, k, iG, f, iRO2, ~)
    y = y(:);
    y(2) = sum(y(iRO2));
    G = y(iG(:,1)) .* y(iG(:,2)) .* y(iG(:,3));
    rates = k .* G;
    dydt = rates' * f;
    dydt = dydt(:);
    dydt(1) = 0;
end

function J = jacobian_full(~, y, k, iG, f, iRO2, nSp)
    y = y(:);
    y(2) = sum(y(iRO2));
    DrDy = compute_drdy(y, k, iG, nSp, iRO2);
    J = f' * DrDy;
    J(1,:) = 0;
end

function cmap = redbluecm_custom
    n = 64;
    cmap = zeros(n, 3);
    for i = 1:n
        t = (i-1)/(n-1);
        if t < 0.5
            cmap(i,:) = [t*2, t*2, 1];
        else
            cmap(i,:) = [1, (1-t)*2, (1-t)*2];
        end
    end
end
