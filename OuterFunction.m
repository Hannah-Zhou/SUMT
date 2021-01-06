%% �Ż��㷨
% ʹ����㷣���������Ķ��η���
clear all;
clc;
%% 
% ����Ȩ������
data1 = xlsread('D:\��ѧ��ѧ��Դ\�Ż��������㷨\�Ż�����ҵ\LibraryData.xlsx', '������ 1 - ��� 3', 'B3: B24');
data2 = xlsread('D:\��ѧ��ѧ��Դ\�Ż��������㷨\�Ż�����ҵ\LibraryData.xlsx', '������ 1 - ��� 3', 'C3: C24');
w1 = data1;
w1 = w1 ./ norm(w1);
w2 = data2;
w2 = w2 ./ norm(w2);
%% 
% ��ʼ������Ȩ������
alpha = [1, 2]';
% alpha_inf = 0;
% alpha_sup = 1;
%% ��ֵ����
sigma = 1e+1;
c = 2;
epsilon1 = 1e-6;
epsilon2 = 1e-8;
step_max = 100;
% ƫִ��
% c1, c2: c1 + c2 -> sum(c1, c2) = 1
[best_alpha, sigma, step] = SUMT_outer(alpha, sigma, c, epsilon1, epsilon2, step_max, w1, w2);
%% ����
disp('�㷨����, ��ӡ����ֵ��');
disp(best_alpha);
fprintf('����ֵ��Ӧ�ĺ���ֵ��%f \n', opt_function(best_alpha, w1, w2));
fprintf('����ֵ��Ӧ�ĳͷ��ɳ�ֵ��%f \n', st_function_new(best_alpha, sigma));
fprintf('����ֵ��Ӧ�ĳͷ����ӣ�%f \n', sigma);
fprintf('��ǰ����������%f \n', step);
disp('----------------------------');

%% ������ֵ����
function f = opt_function(alpha, w1, w2)
% ������ʵ�ֺ�����ֵ
f = 0.4*norm(alpha(1).*w1 + alpha(2).*w2 - w1) + 0.6*norm(alpha(1).*w1 + alpha(2).*w2 - w2);
% f = norm(alpha(1).*w1 + alpha(2).*w2 - w1) + norm(alpha(1).*w1 + alpha(2).*w2 - w2);
end

function [st1, st2, st3] = st_function(alpha)
% ������ʵ�������ⷣԼ����Ӧ����ֵ

st1 = - alpha(1);
st2 = - alpha(2);
st3 = sum(alpha) - 1;
end

function g = st_function_new(alpha, sigma)
% ������ʵ��Լ����Ӧ�Ⱥ�����ֵ
[st1, st2, st3] = st_function(alpha);
g = sigma*(st3^2 + max([0, st1])^2 + max([0, st2])^2) / 2;
end

function f = opt_function_new(alpha, w1, w2, sigma)
% ������ʵ�ֳͷ�������ֵ
f = opt_function(alpha, w1, w2);
g = st_function_new(alpha, sigma);
f = f + g;
end

function grad_f = grad(alpha, w1, w2)
% ������ʵ�����ݶ�ֵ. ����alpha
% ��alpha1�󵼵Ľ��
grad_f1_alpha1 = 0.4 * ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
grad_f2_alpha1 = 0.6 * ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% grad_f1_alpha1 = ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
% grad_f2_alpha1 = ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w1) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% ��alpha2�󵼵Ľ��
grad_f1_alpha2 = 0.4 * ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
grad_f2_alpha2 = 0.6 * ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% grad_f1_alpha2 = ((alpha(1).*w1 + alpha(2).*w2 - w1)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w1));
% grad_f2_alpha2 = ((alpha(1).*w1 + alpha(2).*w2 - w2)'*w2) / (norm(alpha(1).*w1 + alpha(2).*w2 - w2));
% ���
grad_f = [grad_f1_alpha1 + grad_f2_alpha1; grad_f1_alpha2 + grad_f2_alpha2];
end

function grad_f = grad_new(alpha, w1, w2, sigma)
% ������ʵ�����ݶ�ֵ. ����alpha
grad_f = grad(alpha, w1, w2);
% �жϴ�ʱ��st1��st2�Ƿ�����Ǹ�
[st1, st2, ~] = st_function(alpha);
if st1 >= 0
    if st2 >= 0
        grad_g_alpha1 = sigma*(sum(alpha) - 1 + alpha(1));
        grad_g_alpha2 = sigma*(sum(alpha) - 1 + alpha(2));
    else
        grad_g_alpha1 = sigma*(sum(alpha) - 1 + alpha(1));
        grad_g_alpha2 = sigma*(sum(alpha) - 1);
    end
else
    if st2 >= 0
        grad_g_alpha1 = sigma*(sum(alpha) - 1);
        grad_g_alpha2 = sigma*(sum(alpha) - 1 + alpha(2));
    else
        grad_g_alpha1 = sigma*(sum(alpha) - 1);
        grad_g_alpha2 = sigma*(sum(alpha) - 1);
    end
end
grad_g = [grad_g_alpha1; grad_g_alpha2];
grad_f = grad_f + grad_g;
end
%% һά������ģ��
function lamda_f = lamda_function(lamda, alpha, d, w1, w2, sigma)
% ������ʵ��һά��������ֵ. ����lamda. 
lamda_f = opt_function_new(alpha + lamda.*d, w1, w2, sigma);
end
% ----------------------- %
function lamda_grad = lamda_grad_new(lamda, alpha, d, w1, w2, sigma)
% ������ʵ������ֵ. ����lamda. 
q = d(1).*w1 + d(2).*w2;
p1 = alpha(1).*w1 + alpha(2).*w2 - w1;
p2 = alpha(1).*w1 + alpha(2).*w2 - w2;
grad1 = 0.4 * (p1 + lamda.*q)'*q / norm(p1 + lamda.*q) + 0.6 * (p2 + lamda.*q)'*q / norm(p2 + lamda.*q);
% grad1 = (p1 + lamda.*q)'*q / norm(p1 + lamda.*q) + (p2 + lamda.*q)'*q / norm(p2 + lamda.*q);
if alpha(1) + lamda*d(1) < 0
    if alpha(2) + lamda*d(2) < 0
        grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(1) + lamda*d(1))*d(1) + (alpha(2) + lamda*d(2))*d(2);
    else
        grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(1) + lamda*d(1))*d(1);
    end
else
    if alpha(2) + lamda*d(2) < 0
        grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d) + (alpha(2) + lamda*d(2))*d(2);
    else
        grad2 = sigma*(sum(alpha) + lamda*sum(d) - 1)*sum(d);
    end
end
lamda_grad = grad1 + grad2;
end
% ----------------------- %
% function lamda_best = Wolfe(lamda, alpha, d, w1, w2, sigma, max_time)
% % ��lamda�������������Ǿ�ȷ����ȷ��һ�����ʵ�lamda
% % ���ã�Wolfe׼��
% % Ĭ�ϲ�����c1 = 0.5, c2 = 0.9
% % ��ʼlamda: 1
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
%     % ִ�в�ֵ�������ҵ�һ�����κ�����ʹ�����㣺g
%     % g(0) = opt_function_new(alpha, w1, w2, sigma)
%     % g'(0) = lamda_grad_new(0, alpha, d, w1, w2, sigma)
%     % g(lamda) = lamda_function(lamda, alpha, d, w1, w2, sigma)
%     % ȡ��g(x)��(0, lamda)�ļ�Сֵ����Ϊlamda�ĸ��µ����
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
function lamda_best = Newton_spline(lamda, alpha, d, w1, w2, sigma, max_time)
% ��lamda�������������Ǿ�ȷ����ȷ��һ�����ʵ�lamda
% ���ã�����Newton׼��
% Ĭ�ϲ���
% ��ʼlamda: 1
time = 0;
lamda_sub = max([1e+4, 100*sigma]);
lamda_old = 2 * lamda;
if lamda_function(lamda, alpha, d, w1, w2, sigma) < lamda_function(lamda_old, alpha, d, w1, w2, sigma)
    lamda_history = lamda;
else
    lamda_history = lamda_old;
end
while lamda_sub < 1e-6 || time < max_time
    delta_grad = lamda_grad_new(lamda, alpha, d, w1, w2, sigma) - lamda_grad_new(lamda_old, alpha, d, w1, w2, sigma);
    if abs(delta_grad) < 1e-6
        break
    end
    lamda_new = lamda - lamda_grad_new(lamda, alpha, d, w1, w2, sigma) * (lamda - lamda_old) / (delta_grad);
    if lamda_function(lamda_new, alpha, d, w1, w2, sigma) < lamda_function(lamda_history, alpha, d, w1, w2, sigma)
        lamda_history = lamda;
    end
    lamda_sub = abs(lamda_new - lamda);
    if 0.8*lamda_function(lamda_new, alpha, d, w1, w2, sigma) > lamda_function(lamda, alpha, d, w1, w2, sigma)
        break
    end
    lamda_old = lamda;
    lamda = lamda_new;
    time = time + 1;
end
lamda_best = lamda_history;
end
%% ����������ȷ�����ʵ���������
% ʹ��DFP�㷨����
function [min_x, lamda, end_point] = DFP_BFGS(alpha, w1, w2, sigma, epsilon, max_time)
% epsilon: �������
% OPTION: Ĭ�ϲ���Ϊ1('DFS'),Ҳ������Ϊ2('BFGS')
OPTION = 2;
% ��ʱalpha�ǳ�ʼ�����㣬������������alpha����
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
    G = grad_new(alpha, w1, w2, sigma);
    if norm(G) < 1e-6 && time == 0
        end_point = 1;
        break
    end
%     disp('�ݶȣ�');
%     disp(G);
    d = - H*G;
    d = d ./ norm(d);
%     disp('����');
%     disp(d);
%     lamda = Wolfe(lamda, alpha, d, w1, w2, sigma, max_time);
    lamda = Newton_spline(lamda, alpha, d, w1, w2, sigma, max_time);
%     fprintf('��ǰ������%f \n', lamda);
    if lamda < 1e-6
        break
    end
    fprintf('��ǰ���ģ��%f \n', norm(alpha));
    fprintf('��ǰ���Ӧ�ĺ���ֵ��%f \n', opt_function_new(alpha, w1, w2, sigma));
    alpha = alpha + lamda.*d;
%     disp('��ǰ�⣺');
%     disp(alpha);
    if norm(grad_new(alpha, w1, w2, sigma)) < epsilon
        break
    end
    time = time + 1;
    tip = tip + 1;
    if tip == x_size
        H = eye(x_size);
        tip = 0;
        continue
    end
    G_new = grad_new(alpha, w1, w2, sigma);
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
%% �����Ż�ģ��
function [min_x, sigma, step] = SUMT_outer(alpha, sigma, c, epsilon1, epsilon2, step_max, w1, w2)
% alpha0: ��ʼ������
% sigma: ��ʼ������
% c: �Ŵ�ϵ��
% epsilon: �������
% step_max: ���������������
disp('����ʹ�� ��㷣�������� ����SUMT�Ż�...');
disp('----------------------------');
step = 1;
max_time = step_max;
% �洢��ʷ������(figure1: x - y; figure2: step - f(alpha))
history_best_x1 = [];
history_best_x2 = [];
history_best_f_st = [];
history_best_f = [];
%
while true
    if step > step_max
        break
    end
    fprintf('���е�%d���Ż�����... \n', step);
    if step > 1
        st_value_old = st_value;
    end
    f_st_old = opt_function_new(alpha, w1, w2, sigma);
    [alpha, lamda, end_point] = DFP_BFGS(alpha, w1, w2, sigma, epsilon1, max_time);
    f_st_new = opt_function_new(alpha, w1, w2, sigma);
    % �洢��ʷʱ���
    history_best_x1 = [history_best_x1, alpha(1)];
    history_best_x2 = [history_best_x2, alpha(2)];
    history_best_f_st = [history_best_f_st, f_st_new];
    history_best_f = [history_best_f, opt_function(alpha, w1, w2)];
    % ����Լ������ģ
    st_value = abs(st_function_new(alpha, sigma));
    if step == 1
        st_value_old = st_value;
    end
    % ͣ��׼��
    bool1 = (abs(st_function_new(alpha, sigma)) < epsilon2 && end_point == 1);
    bool2 = (abs(lamda) + norm(grad_new(alpha, w1, w2, sigma)) < 1e-6);
    bool3 = (0.8*abs(st_value) > abs(st_value_old));
    bool4 = (0.8*f_st_new > f_st_old);
    if bool1 || bool2 || bool3 || bool4
        disp('----------------------------');
        if bool1 || bool2
            disp('ͣ��׼�򣺺�����ȡ���������ֵ...');
        elseif bool4 || bool3
            disp('ͣ��׼�򣺵�ǰ������̬...');
        end
        fprintf('������ֹʱ, �����ݶȷ���Ϊ%f \n', norm(grad(alpha, w1, w2)));
        break
        disp('----------------------------');
    end
    disp('----------------------------');
    fprintf('��%d���Ż�����ȡ�õĽ�����£� \n', step);
    disp('��ǰ�⣺');
    disp(alpha);
    fprintf('��ǰ���Ӧ��Լ������ģ��%f \n', st_value);
    disp('----------------------------');
    sigma = sigma*c;
    step = step + 1;
end
min_x = alpha;
% ��ͼѡ��
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



   





