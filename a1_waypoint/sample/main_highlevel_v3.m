clc
close all
clear all

PaperPosition = [-0.25 -0.1 8 6];
PaperSize     = [7.25 5.8];
Fontsize      = 12;
print_pdf     = 0;
matlabpath = 'figures/';
Linewidth     = 2;

% ---- Obstacle definitions (keep in sync with cbf_circle_obstacles.c) --------
obs = [ -2.0,  1.0, 0.35;
        -4.0, -1.0, 0.35;
        -2.0, -1.0, 0.35;
        -4.0,  1.0, 0.35 ];

% =============================================================================
% FIGURE 1 — original tracking plots
% =============================================================================
data

i = 1;
step        = all_data(:,i); i = i+1;
x_ref       = all_data(:,i); i = i+1;
y_ref       = all_data(:,i); i = i+1;
theta_ref   = all_data(:,i); i = i+1;
x_p         = all_data(:,i); i = i+1;
y_p         = all_data(:,i); i = i+1;
theta       = all_data(:,i); i = i+1;
error_x     = all_data(:,i); i = i+1;
error_y     = all_data(:,i); i = i+1;
error_theta = all_data(:,i); i = i+1;
cmd1        = all_data(:,i); i = i+1;
cmd2        = all_data(:,i); i = i+1;
cmd3        = all_data(:,i); i = i+1;

hh = figure(1);
subplot(2,2,[1 2]);
plot(x_p, y_p, 'r', 'Linewidth', 2); hold on
plot(x_ref, y_ref, 'k-.', 'Linewidth', 1);
th = linspace(0, 2*pi, 200);
for k = 1:size(obs,1)
    cx = obs(k,1);  cy = obs(k,2);  r = obs(k,3);
    fill(cx + r*cos(th), cy + r*sin(th), 'k', ...
         'EdgeColor','k','FaceAlpha',0.25,'LineWidth',1.5);
    text(cx, cy, sprintf('O%d',k), ...
         'HorizontalAlignment','center', ...
         'FontSize',Fontsize-2,'Color','k','FontWeight','bold');
end
ylabel('y','Fontsize',Fontsize);
xlabel('x','Fontsize',18);
legend('act','ref','Fontsize',Fontsize,'Location','Best');
ax = gca; ax.FontSize = Fontsize;
grid on
ylim([-0.8 0.8])
xlim([-4.2 0.2])

subplot(2,2,4);
plot(step, cmd1, 'b:', 'Linewidth',2); hold on
plot(step, cmd2, 'r',  'Linewidth',2);
plot(step, cmd3, 'k',  'Linewidth',2);
ylabel('command','Fontsize',Fontsize);
legend('v_x','v_y','\omega','Fontsize',Fontsize);
ax = gca; ax.FontSize = Fontsize;
xlabel('step no','Fontsize',Fontsize);

subplot(2,2,3);
plot(step, error_x,     'b:', 'Linewidth',2); hold on
plot(step, error_y,     'r',  'Linewidth',2);
plot(step, error_theta, 'k-.','Linewidth',2);
ylabel('error','Fontsize',Fontsize);
legend('x','y','theta','Fontsize',Fontsize);
ax = gca; ax.FontSize = Fontsize;
xlabel('step no','Fontsize',Fontsize);

if print_pdf
    print(hh, '-dpdf', [path 'all_tracking']);
end

% =============================================================================
% FIGURE 2 — PCBF logged data (h, V_T, gradients)
% reads pcbf_log.csv written by cbf_pcbf_circles_v3.c
% =============================================================================
log_file = 'pcbf_log.csv';

if ~isfile(log_file)
    warning('pcbf_log.csv not found — skipping Figure 2. Run the simulation first.');
    return
end

T = readtable(log_file);

% Expected columns:
%   step, x, y, theta,
%   h0,h1,h2,h3,
%   V0,V1,V2,V3,
%   dVdx0,dVdy0,dVdt0,
%   dVdx1,dVdy1,dVdt1,
%   dVdx2,dVdy2,dVdt2,
%   dVdx3,dVdy3,dVdt3

log_step = T.step;
obs_colors = {'b','r',[0 0.6 0],[0.8 0.4 0]};  % one color per obstacle

hh2 = figure(2);
set(hh2,'PaperPosition',PaperPosition,'PaperSize',PaperSize);

% ---- Subplot 1: Phase portrait with h=0 and V=0 contours -------------------
subplot(2,3,[1 2]);
hold on; grid on;

% Draw obstacles
th = linspace(0, 2*pi, 200);
for k = 1:size(obs,1)
    cx = obs(k,1);  cy = obs(k,2);  r = obs(k,3);
    fill(cx + r*cos(th), cy + r*sin(th), 'k', ...
         'EdgeColor','k','FaceAlpha',0.25,'LineWidth',1.5);
    text(cx, cy, sprintf('O%d',k), ...
         'HorizontalAlignment','center', ...
         'FontSize',Fontsize-2,'Color','k','FontWeight','bold');
end

% Actual trajectory from log
plot(T.x, T.y, 'r-', 'LineWidth', 2);
ylabel('y','FontSize',Fontsize);
xlabel('x','FontSize',Fontsize);
title('Robot trajectory (from pcbf\_log)', 'FontSize', Fontsize);
legend('obstacles','','','','trajectory','FontSize',Fontsize-1,'Location','best');
ax = gca; ax.FontSize = Fontsize;
ylim([-1.5 1.5]);
xlim([-4.5 0.5]);

% ---- Subplot 2: h(x) current state per obstacle ----------------------------
subplot(2,3,3);
hold on; grid on;
for k = 1:4
    col_name = sprintf('h%d', k-1);
    if ismember(col_name, T.Properties.VariableNames)
        plot(log_step, T.(col_name), 'Color', obs_colors{k}, 'LineWidth', Linewidth);
    end
end
yline(0, 'k--', 'LineWidth', 1.5, 'Label', 'h=0');
ylabel('h(x_0)', 'FontSize', Fontsize);
xlabel('step no', 'FontSize', Fontsize);
title('Barrier function h at current state', 'FontSize', Fontsize);
legend('h_0','h_1','h_2','h_3','FontSize',Fontsize-1,'Location','best');
ax = gca; ax.FontSize = Fontsize;

% ---- Subplot 3: V_T per obstacle -------------------------------------------
subplot(2,3,4);
hold on; grid on;
for k = 1:4
    col_name = sprintf('V%d', k-1);
    if ismember(col_name, T.Properties.VariableNames)
        plot(log_step, T.(col_name), 'Color', obs_colors{k}, 'LineWidth', Linewidth);
    end
end
yline(0, 'k--', 'LineWidth', 1.5, 'Label', 'V=0');
ylabel('V_T(x_0)', 'FontSize', Fontsize);
xlabel('step no', 'FontSize', Fontsize);
title('PCBF value function V_T', 'FontSize', Fontsize);
legend('V_0','V_1','V_2','V_3','FontSize',Fontsize-1,'Location','best');
ax = gca; ax.FontSize = Fontsize;

% ---- Subplot 4: h vs V comparison (worst obstacle each step) ---------------
subplot(2,3,5);
hold on; grid on;
% Worst case h and V across all obstacles at each step
h_mat = [T.h0, T.h1, T.h2, T.h3];
V_mat = [T.V0, T.V1, T.V2, T.V3];
h_worst = max(h_mat, [], 2);
V_worst = max(V_mat, [], 2);
plot(log_step, h_worst, 'b-', 'LineWidth', Linewidth);
plot(log_step, V_worst, 'r--','LineWidth', Linewidth);
yline(0, 'k--', 'LineWidth', 1.5);
ylabel('value', 'FontSize', Fontsize);
xlabel('step no', 'FontSize', Fontsize);
title('Worst-case h vs V_T', 'FontSize', Fontsize);
legend('h worst', 'V_T worst', 'FontSize', Fontsize-1, 'Location', 'best');
ax = gca; ax.FontSize = Fontsize;

% ---- Subplot 5: Gradient magnitudes ----------------------------------------
subplot(2,3,6);
hold on; grid on;
for k = 1:4
    dx_name = sprintf('dVdx%d', k-1);
    dy_name = sprintf('dVdy%d', k-1);
    if ismember(dx_name, T.Properties.VariableNames) && ...
       ismember(dy_name, T.Properties.VariableNames)
        grad_mag = sqrt(T.(dx_name).^2 + T.(dy_name).^2);
        plot(log_step, grad_mag, 'Color', obs_colors{k}, 'LineWidth', Linewidth);
    end
end
ylabel('||\nablaV||', 'FontSize', Fontsize);
xlabel('step no', 'FontSize', Fontsize);
title('Gradient magnitude ||\nablaV_T||', 'FontSize', Fontsize);
legend('\nablaV_0','\nablaV_1','\nablaV_2','\nablaV_3', ...
       'FontSize',Fontsize-1,'Location','best');
ax = gca; ax.FontSize = Fontsize;

sgtitle('PCBF Logged Data — h, V_T, Gradients', 'FontSize', Fontsize+1);

if print_pdf
    print(hh2, '-dpdf', [path 'pcbf_log_plots']);
end
