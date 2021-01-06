%% 优化算法
% 使用Lagrange乘子方法
clear all;
clc;
%% 
% 输入权重向量
data1 = xlsread('D:\大学数学资源\优化理论与算法\优化大作业\LibraryData.xlsx', '工作表 1 - 表格 3', 'B3: B24');
data2 = xlsread('D:\大学数学资源\优化理论与算法\优化大作业\LibraryData.xlsx', '工作表 1 - 表格 3', 'C3: C24');
w1 = data1;
w1 = w1 ./ norm(w1);
w2 = data2;
w2 = w2 ./ norm(w2);
%% 
% 初始化最优权重因子
alpha = [1, 2]';
% alpha_inf = 0;
% alpha_sup = 1;
%% 求值区域
sigma = 1e+1;
rho = 2;
epsilon = 1e-6;
yeta = 1e-4;
step_max = 100;
stl1 = 0.6;  % 0 < stl1 <= stl2 <= 1
stl2 = 0.9;
[best_alpha, sigma, t, k1, k2, step] = SUMT_outer(alpha, sigma, rho, epsilon, yeta, step_max, w1, w2, stl1, stl2);
%% 后处理
disp('算法结束, 打印最优值：');
disp(best_alpha);
fprintf('最优值对应的函数值：%f \n', opt_function(best_alpha, w1, w2));
fprintf('最优值对应的乘子：t = %f, k1 = %f, k2 = %f \n', t, k1, k2);
fprintf('最优值对应的惩罚因子：%f \n', sigma);
fprintf('当前迭代步数：%f \n', step);
disp('----------------------------');

%% 函数求值区域
function f = opt_function(alpha, w1, w2)
% 本函数实现函数求值
f = 0.4*norm(alpha(1).*w1 + alpha(2).*w2 - w1) + 0.6*norm(alpha(1).*w1 + alpha(2).*w2 - w2);
% f = norm(alpha(1).*w1 + alpha(2).*w2 - w1) + norm(alpha(1).*w1 + alpha(2).*w2 - w2);
end

% function [st1, st2, st3] = st_function(alpha)
% % 本函数实现三种外罚约束适应度求值
% st1 = - alpha(1);
% st2 = - alpha(2);
% st3 = sum(alpha) - 1;
% end

% function g = st_function_new(alpha, sigma)
% % 本函数实现约束适应度函数求值
% [st1, st2, st3] = st_function(alpha);
% g = sigma*(st3^2 + max([0, st1])^2 + max([0, st2])^2) / 2;
% end

% function f = opt_function_new(alpha, w1, w2, sigma)
% % 本函数实现惩罚函数求值
% f = opt_function(alpha, w1, w2);
% g = st_function_new(alpha, sigma);
% f = f + g;
% end

function g = st_function(alpha, sigma, gamma)
% 本函数实现罚函数（约束适应度函数）求值
% gamma: 乘子
g = ((sum(alpha) - 1)^2 + (-alpha(1) + gamma(1))^2 + (-alpha(2) + gamma(2))^2) * sigma / 2;
end

function f = lagrange_function(alpha, w1, w2, sigma, gamma, t, k1, k2)
% 本函数实现lagrange乘子函数求值
% gamma: 松弛变量
% t: 等式约束乘子
% k1, k2: 含松弛变量的等式约束乘子
f = opt_function(alpha, w1, w2) + t*(sum(alpha) - 1) + k1*(-alpha(1) + gamma(1)) + k2*(-alpha(2) + gamma(2));
f = f + st_function(alpha, sigma, gamma);
end

function g = lagrange_st_function(alpha, w1, w2, sigma, t, k1, k2)
% 本函数实现lagrange约束函数罚函数求值(约减形式)
% t: 等式约束乘子
% k1, k2: 含松弛变量的等式约束乘子
g = opt_function(alpha, w1, w2) + t*(sum(alpha) - 1) + sigma * (sum(alpha) - 1)^2 / 2;
g = g + sigma * (max([0, k1 / sigma - alpha(1)])^2 + max([0, k2 / sigma - alpha(2)])^2 - (k1^2 + k2^2) / (sigma^2)) / 2;
end

function f = lagrange_function_exs(alpha, w1, w2, sigma, t, k1, k2)
% 本函数实现lagrange乘子函数求值(约减形式)
% t: 等式约束乘子
% k1, k2: 含松弛变量的等式约束乘子
f = opt_function(alpha, w1, w2) + lagrange_st_function(alpha, w1, w2, sigma, t, k1, k2);
end

function grad_f = grad(alpha, w1, w2)
% 本函数实现求梯度值. 对象：alpha
% 对alpha1求导的结果
grad_f1_alpha1 = 0.4 * ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
grad_f2_alpha1 = 0.6 * ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% grad_f1_alpha1 = ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
% grad_f2_alpha1 = ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% 对alpha2求导的结果
grad_f1_alpha2 = 0.4 * ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
grad_f2_alpha2 = 0.6 * ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% grad_f1_alpha2 = ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
% grad_f2_alpha2 = ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% 组合
grad_f = [grad_f1_alpha1 + grad_f2_alpha1; grad_f1_alpha2 + grad_f2_alpha2];
end

function grad_f = grad_lagrange(alpha, w1, w2, sigma, t, k1, k2)
% 本函数实现求拉格朗日函数的梯度值
% 对象：alpha
grad_f = grad(alpha, w1, w2);
grad_g1_alpha = t .* [1; 1];
grad_g2_alpha = sigma*(sum(alpha) - 1) .* [1; 1];
if (k1 / sigma - alpha(1)) > 0
    if (k2 / sigma - alpha(2)) > 0
        grad_g3_alpha = [-sigma*(k1 / sigma - alpha(1)); -sigma*(k2 / sigma - alpha(2))];
    else
        grad_g3_alpha = [-sigma*(k1 / sigma - alpha(1)); 0];
    end
else
    if (k2 / sigma - alpha(2)) > 0
        grad_g3_alpha = [0; -sigma*(k2 / sigma - alpha(2))];
    else
        grad_g3_alpha = [0; 0];
    end
end
grad_f = grad_f + grad_g1_alpha + grad_g2_alpha + grad_g3_alpha;
end

% function grad_f = grad_new(alpha, w1, w2, sigma)
% % 本函数实现求梯度值. 对象：alpha
% grad_f = grad(alpha, w1, w2);
% % 判断此时的st1与st2是否满足非负
% [st1, st2, ~] = st_function(alpha);
% if st1 >= 0
%     if st2 >= 0
%         grad_g_alpha1 = sigma*(sum(alpha) - 1 + alpha(1));
%         grad_g_alpha2 = sigma*(sum(alpha) - 1 + alpha(2));
%     else
%         grad_g_alpha1 = sigma*(sum(alpha) - 1 + alpha(1));
%         grad_g_alpha2 = sigma*(sum(alpha) - 1);
%     end
% else
%     if st2 >= 0
%         grad_g_alpha1 = sigma*(sum(alpha) - 1);
%         grad_g_alpha2 = sigma*(sum(alpha) - 1 + alpha(2));
%     else
%         grad_g_alpha1 = sigma*(sum(alpha) - 1);
%         grad_g_alpha2 = sigma*(sum(alpha) - 1);
%     end
% end
% grad_g = [grad_g_alpha1; grad_g_alpha2];
% grad_f = grad_f + grad_g;
% end

function v = lagrange_st_eq(alpha, k1, k2, sigma)
% 本函数求解对应Lagrange问题的约束违反度指标
v = sqrt((sum(alpha) - 1)^2 + max([-alpha(1), -k1 / sigma])^2 + max([-alpha(2), -k2 / sigma])^2);
end


%% 一维线搜索模型
% function lamda_f = lamda_function(lamda, alpha, d, w1, w2, sigma)
% % 本函数实现一维线搜索求值. 对象：lamda. 
% lamda_f = opt_function_new(alpha + lamda.*d, w1, w2, sigma);
% end

function lamda_f = lamda_lagrange_function(lamda, alpha, d, w1, w2, sigma, t, k1, k2)
% 本函数实现lagrange乘子函数求值(约减形式)
% gamma: 松弛变量
% t: 等式约束乘子
% k1, k2: 含松弛变量的等式约束乘子
lamda_f = lagrange_function_exs(alpha + lamda.*d, w1, w2, sigma, t, k1, k2);
end
% ----------------------- %
% function lamda_grad = lamda_grad_new(lamda, alpha, d, w1, w2, sigma)
% % 本函数实现求导数值. 对象：lamda. 
% q = d(1).*w1 + d(2).*w2;
% p1 = alpha(1).*w1 + alpha(2).*w2 - w1;
% p2 = alpha(1).*w1 + alpha(2).*w2 - w2;
% grad1 = 0.4 * (p1 + lamda.*q)'*q / norm(p1 + lamda.*q) + 0.6 * (p2 + lamda.*q)'*q / norm(p2 + lamda.*q);
% if alpha(1) + lamda*d(1) < 0
%     if alpha(2) + lamda*d(2) < 0
%         grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(1) + lamda*d(1))*d(1) + (alpha(2) + lamda*d(2))*d(2);
%     else
%         grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(1) + lamda*d(1))*d(1);
%     end
% else
%     if alpha(2) + lamda*d(2) < 0
%         grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(2) + lamda*d(2))*d(2);
%     else
%         grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d);
%     end
% end
% lamda_grad = grad1 + grad2;
% end
function lamda_grad = lamda_lagrange_grad(lamda, alpha, d, w1, w2, sigma, t, k1, k2)
% 本函数实现求导数值.
% 目标对象：lamda
q = d(1).*w1 + d(2).*w2;
p1 = alpha(1).*w1 + alpha(2).*w2 - w1;
p2 = alpha(1).*w1 + alpha(2).*w2 - w2;
grad1 = 0.4 * (p1 + lamda.*q)'*q / norm(p1 + lamda.*q) + 0.6 * (p2 + lamda.*q)'*q / norm(p2 + lamda.*q);
% grad1 = (p1 + lamda.*q)'*q / norm(p1 + lamda.*q) + (p2 + lamda.*q)'*q / norm(p2 + lamda.*q);
grad21 = t * (sum(d));
grad22 = sigma * (sum(alpha) + sum(d)*lamda - 1)*sum(d);
if (k1 / sigma - alpha(1)) > 0
    if (k2 / sigma - alpha(2)) > 0
        grad23 = -sigma*(k1/sigma - alpha(1) - lamda*d(1))*d(1) - sigma*(k2/sigma - alpha(2) - lamda*d(2))*d(2);
    else
        grad23 = -sigma*(k1/sigma - alpha(1) - lamda*d(1))*d(1);
    end
else
    if (k2 / sigma - alpha(2)) > 0
        grad23 = -sigma*(k2/sigma - alpha(2) - lamda*d(2))*d(2);
    else
        grad23 = 0;
    end
end
lamda_grad = grad1 + grad21 + grad22 + grad23;
end
% ----------------------- %
% function lamda_best = Wolfe(lamda, alpha, d, w1, w2, sigma, max_time)
% % 对lamda进行线搜索（非精确）而确定一个合适的lamda
% % 采用：Wolfe准则
% % 默认参数：c1 = 0.5, c2 = 0.9
% % 初始lamda: 1
% c1 = 0.2;
% c2 = 0.9;
% time = 0;
% while true
%     choose1 = lamda_function(lamda, alpha, d, w1, w2, sigma);
%     choose1 = choose1 - opt_function_new(alpha, w1, w2, sigma) - (c1*lamda).*grad_new(alpha, w1, w2, sigma)'*d;
%     choose2 = c2.*grad_new(alpha, w1, w2, sigma)'*d - grad_new(alpha + lamda.*d, w1, w2, sigma)'*d;
%     if (choose1 <= 0 && choose2 <= 0) || abs(lamda) < 1e-6 || time > max_time
%         break
%     end
%     lamda = lamda / 2;
%     % 执行插值操作。找到一个二次函数，使得满足：g
%     % g(0) = opt_function_new(alpha, w1, w2, sigma)
%     % g'(0) = lamda_grad_new(0, alpha, d, w1, w2, sigma)
%     % g(lamda) = lamda_function(lamda, alpha, d, w1, w2, sigma)
%     % 取得g(x)在(0, lamda)的极小值点作为lamda的更新点继续
% %     b = lamda_grad_new(0, alpha, d, w1, w2, sigma);
% %     c = opt_function_new(alpha, w1, w2, sigma);
% %     a = (choose1 - b*lamda - c) / lamda^2;
% %     if abs(a) < 1e-6
% %         if b > 0
% %             lamda = lamda/2;
% %         end
% %     else
% %         be_func = - b/(2*a);
% %         if be_func > 0 && be_func < abs(lamda)
% % %             lamda = min([0, lamda, a*be_func^2 + b*be_func + c]);
% %             lamda = min([lamda/2, lamda, a*be_func^2 + b*be_func + c]);
% %         else
% % %             lamda = min([0, lamda]);
% %             lamda = min([lamda/2, lamda]);
% %         end
% %     end
%     time = time + 1;
% end
% lamda_best = lamda;
% end
% ----------------------- %
function lamda_best = Newton_spline(lamda, alpha, d, w1, w2, sigma, max_time, t, k1, k2)
% 对lamda进行线搜索（非精确）而确定一个合适的lamda
% 采用：割线Newton准则
% 默认参数
% 初始lamda: 1
time = 0;
lamda_sub = max([1e+4, 100*sigma]);
lamda_old = 2 * lamda;
if lamda_lagrange_function(lamda, alpha, d, w1, w2, sigma, t, k1, k2) < lamda_lagrange_function(lamda_old, alpha, d, w1, w2, sigma, t, k1, k2)
    lamda_history = lamda;
else
    lamda_history = lamda_old;
end
while lamda_sub < 1e-6 || time < max_time
    delta_grad = lamda_lagrange_grad(lamda, alpha, d, w1, w2, sigma, t, k1, k2) - lamda_lagrange_grad(lamda_old, alpha, d, w1, w2, sigma, t, k1, k2);
    if abs(delta_grad) < 1e-6
        break
    end
    lamda_new = lamda - lamda_lagrange_grad(lamda, alpha, d, w1, w2, sigma, t, k1, k2) * (lamda - lamda_old) / (delta_grad);
    if lamda_lagrange_grad(lamda_new, alpha, d, w1, w2, sigma, t, k1, k2) < lamda_lagrange_grad(lamda_history, alpha, d, w1, w2, sigma, t, k1, k2)
        lamda_history = lamda;
    end
    lamda_sub = abs(lamda_new - lamda);
    if 0.8*lamda_lagrange_grad(lamda_new, alpha, d, w1, w2, sigma, t, k1, k2) > lamda_lagrange_grad(lamda, alpha, d, w1, w2, sigma, t, k1, k2)
        break
    end
    lamda_old = lamda;
    lamda = lamda_new;
    time = time + 1;
end
lamda_best = lamda_history;
end
%% 基于线搜索确定合适的搜索方向
% 使用DFP算法进行
function [min_x, lamda, end_point] = DFP_BFGS(alpha, w1, w2, sigma, yeta0, max_time, t, k1, k2)
% epsilon: 允许误差
% OPTION: 默认参数为1('DFS'),也可设置为2('BFGS')
OPTION = 2;
% 此时alpha是初始迭代点，整个迭代基于alpha进行
time = 0;
tip = 0;
x_size = length(alpha);
H = eye(x_size);
lamda = 1;
end_point = 0;
while true
    if time > max_time
        break
    end
    G = grad_lagrange(alpha, w1, w2, sigma, t, k1, k2);
    if norm(G) < 1e-6 && time == 0
        end_point = 1;
        break
    end
%     disp('梯度：');
%     disp(G);
    d = - H*G;
    d = d ./ norm(d);
%     disp('方向：');
%     disp(d);
%     lamda = Wolfe(lamda, alpha, d, w1, w2, sigma, max_time);
    lamda = Newton_spline(lamda, alpha, d, w1, w2, sigma, max_time, t, k1, k2);
%     fprintf('当前步长：%f \n', lamda);
    if lamda < 1e-6
        break
    end
    fprintf('当前解的模：%f \n', norm(alpha));
    fprintf('当前解对应的函数值：%f \n', lagrange_function_exs(alpha, w1, w2, sigma, t, k1, k2));
    alpha = alpha + lamda.*d;
%     disp('当前解：');
%     disp(alpha);
    if norm(grad_lagrange(alpha, w1, w2, sigma, t, k1, k2)) < yeta0
        break
    end
    time = time + 1;
    tip = tip + 1;
    if tip == x_size
        H = eye(x_size);
        tip = 0;
        continue
    end
    G_new = grad_lagrange(alpha, w1, w2, sigma, t, k1, k2);
    P = lamda.*d;
    Q = G_new - G;
    if OPTION == 1
        H = H + (P*P') ./ (P'*Q) - (H*(Q*Q')*H) ./ (Q'*H*Q);
    elseif OPTION == 2
        H = H + (1 + (Q'*H*Q) ./ (P'*Q)) .* (P*P') ./ (P'*Q) - (P*Q'*H + H*Q*P') ./ (P'*Q);
    end
end
min_x = alpha;
end
%% 序列优化模型
function [min_x, sigma, t, k1, k2, step] = SUMT_outer(alpha, sigma, rho, epsilon, yeta, step_max, w1, w2, stl1, stl2)
% alpha0: 初始迭代点
% sigma: 初始罚因子
% c: 放大系数
% epsilon: 允许误差
% step_max: 允许迭代的最大次数
% 0 < stl1 <= stl2 < 1
% rho = 2;
t = 0.8;
k1 = 0.8;
k2 = 0.8;
yeta0 = 1 / sigma;
epsilon0 = yeta^stl1;
disp('正在使用 增广Lagrange方法 进行SUMT优化...');
disp('----------------------------');
step = 1;
max_time = step_max;
% 存储历史迭代点(figure1: x - y; figure2: step - f(alpha))
history_best_x1 = [];
history_best_x2 = [];
history_best_f_st = [];
history_best_f = [];
% 
while true
    if step > step_max
        break
    end
    fprintf('进行第%d次优化过程... \n', step);
    if step > 1
        st_value_old = st_value;
    end
    f_st_old = lagrange_function_exs(alpha, w1, w2, sigma, t, k1, k2);
    [alpha, lamda, end_point] = DFP_BFGS(alpha, w1, w2, sigma, yeta0, max_time, t, k1, k2);
%     f_st_new = lagrange_function_exs(alpha, w1, w2, sigma, t, k1, k2);
    f_st_new = lagrange_st_eq(alpha, k1, k2, sigma);
    % 存储历史时间点
    history_best_x1 = [history_best_x1, alpha(1)];
    history_best_x2 = [history_best_x2, alpha(2)];
    history_best_f_st = [history_best_f_st, f_st_new];
    history_best_f = [history_best_f, opt_function(alpha, w1, w2)];
    % 停机或更新乘子继续
    if lagrange_st_eq(alpha, k1, k2, sigma) < epsilon0
        if lagrange_st_eq(alpha, k1, k2, sigma) < epsilon && norm(grad_lagrange(alpha, w1, w2, sigma, t, k1, k2)) < yeta
            disp('停机准则：函数已经取得逼近解...');
            fprintf('迭代终止时, 函数梯度模为%f \n', norm(grad_lagrange(alpha, w1, w2, sigma, t, k1, k2)));
            disp('----------------------------');
            break
        else
            % 更新乘子
            t = t + sigma*(sum(alpha) - 1);
            k1 = max([k1 - sigma*alpha(1), 0]);
            k2 = max([k2 - sigma*alpha(2), 0]);
            yeta0 = yeta0 / sigma;
            epsilon0 = epsilon0 / (sigma^stl2);
        end
        
    else
        sigma = sigma * rho;
        yeta0 = 1 / (sigma^stl1);
        epsilon0 = 1 / sigma;
    end
    % 约束函数模
    st_value = abs(lagrange_st_function(alpha, w1, w2, sigma, t, k1, k2));
    if step == 1
        st_value_old = st_value;
    end
    % 停机准则
%     bool1 = (abs(st_function_new(alpha, sigma)) < epsilon2 && end_point == 1);
%     bool2 = (abs(lamda) + norm(grad_new(alpha, w1, w2, sigma)) < 1e-6);
    bool3 = (0.8*abs(st_value) > abs(st_value_old));
    bool4 = (0.8*f_st_new > f_st_old);
%     if bool1 || bool2 || bool3 || bool4
%         disp('----------------------------');
%         if bool1 || bool2
%             disp('停机准则：函数已取得相对最优值...');
%         elseif bool4 || bool3
%             disp('停机准则：当前迭代病态...');
%         end
%         fprintf('迭代终止时, 函数梯度模为%f \n', norm(grad(alpha, w1, w2)));
%         break
%         disp('----------------------------');
%     end
    if bool3 || bool4
        disp('停机准则：当前迭代病态...');
        fprintf('迭代终止时, 函数梯度模为%f \n', norm(grad(alpha, w1, w2)));
        disp('----------------------------');
        break
    end
    disp('----------------------------');
    fprintf('第%d次优化过程取得的结果如下： \n', step);
    disp('当前解：');
    disp(alpha);
    fprintf('当前解对应的约束函数模：%f \n', st_value);
    fprintf('当前解对应的原函数值：%f \n', opt_function(alpha, w1, w2));
    disp('----------------------------');
    step = step + 1;
end
min_x = alpha;
% 绘图选项
figure
plot(history_best_x1, history_best_x2, '-*', 'LineWidth', 2, 'MarkerSize', 5, 'color', [0.5, 0.5, 0]);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Alpha - (1)', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('Alpha - (2)', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Alpha Distribution Map', 'Fontname', 'Times New Roman','FontSize', 18);
figure
plot(1:step, history_best_f_st, '-*', 'LineWidth', 2, 'MarkerSize', 5, 'color', [0.5, 0, 0.5]);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('Function (with S.T.)', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Function (with S.T.)', 'Fontname', 'Times New Roman','FontSize', 18);
figure
plot(1:step, history_best_f, '-*', 'LineWidth', 2, 'MarkerSize', 5);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('Function', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Function', 'Fontname', 'Times New Roman','FontSize', 18);
end
