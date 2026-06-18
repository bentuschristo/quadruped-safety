clc
close all
clear all

PaperPosition = [-0.25 -0.1 8 6]; %location on printed page. rect = [left, bottom, width, height]
PaperSize = [7.25 5.8]; %[width height]
Fontsize = 12;
print_pdf = 1;
path = 'figures/';
Linewidth = 2;

% ---- Obstacle definitions (keep in sync with cbf_circle_obstacles.c) --
obs = [ -2.0,  1, 0.35;   % [cx, cy, r]
        -4.0, -1, 0.35;
        -2.0, -1, 0.35; 
        -4.0,  1, 0.35 ];
% -----------------------------------------------------------------------

data
%data_turn
%data_yaw

i = 1;
step      = all_data(:,i); i = i+1;
x_ref     = all_data(:,i); i = i+1;
y_ref     = all_data(:,i); i = i+1;
theta_ref = all_data(:,i); i = i+1;
x_p       = all_data(:,i); i = i+1;
y_p       = all_data(:,i); i = i+1;
theta     = all_data(:,i); i = i+1;
error_x   = all_data(:,i); i = i+1;
error_y   = all_data(:,i); i = i+1;
error_theta = all_data(:,i); i = i+1;
cmd1      = all_data(:,i); i = i+1;
cmd2      = all_data(:,i); i = i+1;
cmd3      = all_data(:,i); i = i+1;

%%%%%%%%%%%%%%%% alternate figures %%%%%%%%%
hh = figure(1);

subplot(2,2,[1 2]);
plot(x_p, y_p, 'r', 'Linewidth', 2); hold on
plot(x_ref, y_ref, 'k-.', 'Linewidth', 1);

% Draw obstacles as filled black circles
th = linspace(0, 2*pi, 200);
for k = 1:size(obs,1)
    cx = obs(k,1);
    cy = obs(k,2);
    r  = obs(k,3);
    fill(cx + r*cos(th), cy + r*sin(th), 'k', ...
         'EdgeColor', 'k', 'FaceAlpha', 0.25, 'LineWidth', 1.5);
    % label each obstacle
    text(cx, cy, sprintf('O%d', k), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', Fontsize-2, 'Color', 'k', 'FontWeight', 'bold');
end

ylabel('y', 'Fontsize', Fontsize);
xlabel('x', 'Fontsize', 18);
legend('act', 'ref', 'Fontsize', Fontsize, 'Location', 'Best');
ax = gca; ax.FontSize = Fontsize;
grid on
ylim([-0.8 0.8])
xlim([-4.2 0.2])

subplot(2,2,4);
plot(step, cmd1, 'b:', 'Linewidth', 2); hold on
plot(step, cmd2, 'r',  'Linewidth', 2);
plot(step, cmd3, 'k',  'Linewidth', 2);
ylabel('command', 'Fontsize', Fontsize);
legend('v_x', 'v_y', '\omega', 'Fontsize', Fontsize);
ax = gca; ax.FontSize = Fontsize;
xlabel('step no', 'Fontsize', Fontsize);

subplot(2,2,3);
plot(step, error_x,     'b:', 'Linewidth', 2); hold on
plot(step, error_y,     'r',  'Linewidth', 2);
plot(step, error_theta, 'k-.','Linewidth', 2);
ylabel('error', 'Fontsize', Fontsize);
legend('x', 'y', 'theta', 'Fontsize', Fontsize);
ax = gca; ax.FontSize = Fontsize;
xlabel('step no', 'Fontsize', Fontsize);

string = [path, 'all_tracking'];

if (print_pdf == 1)
    print(hh, '-dpdf', string);
end


%% ================= Figure 2: nominal, backup, and filtered controls =================
% The C controller log is expected to have 12 columns:
% step,
% vx_safe, vy_safe, wz_safe,
% vx_nom,  vy_nom,  wz_nom,
% vx_backup, vy_backup, wz_backup,
% mu, status
%
% This file is usually generated as data.txt by my_controller.c.
% If your filename is different, change ctrl_file below.

ctrl_file = 'data.txt';
if exist(ctrl_file, 'file') == 2
    ctrl_data = readmatrix(ctrl_file);

    if size(ctrl_data, 2) >= 12
        cstep = ctrl_data(:,1);

        vx_safe = ctrl_data(:,2);
        vy_safe = ctrl_data(:,3);
        wz_safe = ctrl_data(:,4);

        vx_nom = ctrl_data(:,5);
        vy_nom = ctrl_data(:,6);
        wz_nom = ctrl_data(:,7);

        vx_backup = ctrl_data(:,8);
        vy_backup = ctrl_data(:,9);
        wz_backup = ctrl_data(:,10);

        mu_ctrl = ctrl_data(:,11);
        status_ctrl = ctrl_data(:,12);
        qp_time_sec = ctrl(:,13);
        
        hh2 = figure(2);
        set(hh2, 'PaperUnits', 'inches');
        set(hh2, 'PaperPosition', PaperPosition);
        set(hh2, 'PaperSize', PaperSize);

        % ---- v_x ----
        subplot(3,1,1);
        plot(cstep, vx_nom, 'k--', 'LineWidth', 1.6); hold on
        plot(cstep, vx_backup, 'b:', 'LineWidth', 2.0);
        plot(cstep, vx_safe, 'r', 'LineWidth', 2.0);
        ylabel('v_x', 'Fontsize', Fontsize);
        legend('nominal', 'backup', 'filtered', 'Fontsize', Fontsize-2, 'Location', 'Best');
        ax = gca; ax.FontSize = Fontsize;
        grid on

        % ---- v_y ----
        subplot(3,1,2);
        plot(cstep, vy_nom, 'k--', 'LineWidth', 1.6); hold on
        plot(cstep, vy_backup, 'b:', 'LineWidth', 2.0);
        plot(cstep, vy_safe, 'r', 'LineWidth', 2.0);
        ylabel('v_y', 'Fontsize', Fontsize);
        legend('nominal', 'backup', 'filtered', 'Fontsize', Fontsize-2, 'Location', 'Best');
        ax = gca; ax.FontSize = Fontsize;
        grid on

        % ---- omega ----
        subplot(3,1,3);
        plot(cstep, wz_nom, 'k--', 'LineWidth', 1.6); hold on
        plot(cstep, wz_backup, 'b:', 'LineWidth', 2.0);
        plot(cstep, wz_safe, 'r', 'LineWidth', 2.0);
        ylabel('\omega', 'Fontsize', Fontsize);
        xlabel('step no', 'Fontsize', Fontsize);
        legend('nominal', 'backup', 'filtered', 'Fontsize', Fontsize-2, 'Location', 'Best');
        ax = gca; ax.FontSize = Fontsize;
        grid on

        sgtitle('Nominal, Backup, and Filtered Control Inputs', 'Fontsize', Fontsize+1);

        if ~exist(path, 'dir')
            mkdir(path);
        end

        string2 = [path, 'control_inputs_nominal_backup_filtered'];
        if (print_pdf == 1)
            print(hh2, '-dpdf', string2);
        end
    else
        warning('Control log %s has %d columns. Expected at least 12 columns.', ctrl_file, size(ctrl_data,2));
    end
else
    warning('Control log file %s not found. Figure 2 was not generated.', ctrl_file);
end
