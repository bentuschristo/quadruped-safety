%% compare_backup_methods_update.m
% Compare Backup CBF, Blended Inputs, and Optimal Interpolation
% for two backup-controller scenarios:
%   Figure 1: Backup in Y
%   Figure 2: Zero-Control
%
% Expected directory structure relative to this script:
%
% results/Backup CBFs/update/
%   Backup CBFs - Baseline/
%       Backup in Y/data.csv
%       Backup in Y/data.txt
%       Zero-Control/data.csv
%       Zero-Control/data.txt
%   Blended Inputs/
%       Backup in Y/data.csv
%       Backup in Y/data.txt
%       Zero-Control/data.csv
%       Zero-Control/data.txt
%   Optimal Interpolation/
%       Backup in Y/data.csv
%       Backup in Y/data.txt
%       Zero-Control/data.csv
%       Zero-Control/data.txt
%
% data.csv is assumed to contain the trajectory/state log.
% data.txt is assumed to contain the controller log:
%   col 1:  step
%   col 2:  vx_safe
%   col 3:  vy_safe
%   col 4:  wz_safe
%   col 5:  vx_nom
%   col 6:  vy_nom
%   col 7:  wz_nom
%   col 8:  vx_backup
%   col 9:  vy_backup
%   col 10: wz_backup
%   col 11: mu
%   col 12: status
%   col 13: qp_time_sec, optional

clear; clc; close all;

%% ========================= USER SETTINGS =========================

% If this script is placed inside /results/Backup CBFs/update/, leave this:
root_dir = pwd;

% Otherwise, uncomment and edit this absolute path:
% root_dir = '/sample/results/Backup CBFs/update';
% root_dir = '/home/bgerodac/quadruped/a1_waypoint/sample/results/Backup CBFs/update';

save_figures = true;
fig_dir = fullfile(root_dir, 'figures_compare_methods');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

methods = { ...
    'Backup CBF - Baseline', ...
    'Blended Inputs', ...
    'Optimal Interpolation' ...
};

method_labels = { ...
    'bCBF-QP', ...
    'Blended Inputs', ...
    'Optimal Interpolation' ...
};

scenarios = { ...
    'Backup in Y', ...
    'Zero-Control' ...
};

scenario_titles = { ...
    'Backup in Y', ...
    'Zero-Control Backup' ...
};

% Plot style
Linewidth = 2.0;
Fontsize = 13;

colors = { ...
    [0.00 0.00 0.00], ...   % bCBF black
    [0.00 0.45 0.74], ...   % blended blue
    [0.47 0.67 0.19] ...    % OI green
};

linestyles = {'-', '--', '-.'};

% Obstacles used in the current examples. Edit if your obstacle layout differs.
obs = [ ...
    -2.0,  1.0, 0.35;
    -4.0, -1.0, 0.35;
    -2.0, -1.0, 0.35;
    -4.0,  1.0, 0.35
];

% If you still include the center obstacles, uncomment these:
% obs = [obs;
%     -2.8,  0.0, 0.30;
%     -3.3,  0.0, 0.30];

draw_obstacles = true;

%% ========================= LOAD DATA =========================

data = struct();

for s = 1:numel(scenarios)
    scenario = scenarios{s};

    for m = 1:numel(methods)
        method = methods{m};

        folder = fullfile(root_dir, method, scenario);
        csv_file = fullfile(folder, 'data.csv');
        txt_file = fullfile(folder, 'data.txt');

        entry = struct();
        entry.method = method;
        entry.method_label = method_labels{m};
        entry.scenario = scenario;
        entry.folder = folder;

        if exist(csv_file, 'file')
            entry.traj = readmatrix(csv_file);
        else
            warning('Missing trajectory file: %s', csv_file);
            entry.traj = [];
        end

        if exist(txt_file, 'file')
            entry.ctrl = readmatrix(txt_file);
        else
            warning('Missing controller file: %s', txt_file);
            entry.ctrl = [];
        end

        data(s,m).entry = entry;
    end
end

%% ========================= FIGURE 1 and FIGURE 2 =========================
% Each figure compares the three methods for one backup-controller scenario.
% Layout:
%   (1) x-y trajectory
%   (2) vx filtered control
%   (3) vy filtered control
%   (4) omega filtered control
%   (5) mu / interpolation parameter
%   (6) solve time, if available

for s = 1:numel(scenarios)

    fig = figure(s);
    clf;
    set(fig, 'Color', 'w', 'Position', [100 100 1250 850]);

    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    %% ----- Trajectory subplot -----
    nexttile(1);
    hold on; grid on; box on;

    for m = 1:numel(methods)
        entry = data(s,m).entry;
        traj = entry.traj;

        if isempty(traj)
            continue;
        end

        [x_col, y_col] = detect_xy_columns(traj);

        plot(traj(:,x_col), traj(:,y_col), ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end

    if draw_obstacles
        plot_obstacles(obs);
    end

    axis equal;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title(sprintf('%s: trajectory comparison', scenario_titles{s}), ...
        'Interpreter', 'none', 'FontSize', Fontsize);
    legend('Location', 'best', 'FontSize', 10);

    %% ----- vx subplot -----
    nexttile(2);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        entry = data(s,m).entry;
        ctrl = entry.ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        step = ctrl(:,1);
        vx_safe = ctrl(:,2);
        plot(step, vx_safe, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$v_x^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $v_x$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    legend('Location', 'best', 'FontSize', 10);

    %% ----- vy subplot -----
    nexttile(3);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        entry = data(s,m).entry;
        ctrl = entry.ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        step = ctrl(:,1);
        vy_safe = ctrl(:,3);
        plot(step, vy_safe, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$v_y^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $v_y$', 'Interpreter', 'latex', 'FontSize', Fontsize);

    %% ----- omega subplot -----
    nexttile(4);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        entry = data(s,m).entry;
        ctrl = entry.ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        step = ctrl(:,1);
        wz_safe = ctrl(:,4);
        plot(step, wz_safe, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\omega^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $\omega$', 'Interpreter', 'latex', 'FontSize', Fontsize);

    %% ----- mu subplot -----
    nexttile(5);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        entry = data(s,m).entry;
        ctrl = entry.ctrl;
        if isempty(ctrl) || size(ctrl,2) < 11
            continue;
        end
        step = ctrl(:,1);
        mu = ctrl(:,11);

        % bCBF uses mu = -1 by convention; hide it from this plot.
        mu(mu < 0) = NaN;

        plot(step, mu, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\mu$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Interpolation / blending parameter', 'Interpreter', 'latex', 'FontSize', Fontsize);
    ylim([-0.05, 1.05]);
    legend('Location', 'best', 'FontSize', 10);

    %% ----- computation time subplot -----
    nexttile(6);
    hold on; grid on; box on;
    has_time = false;

    for m = 1:numel(methods)
        entry = data(s,m).entry;
        ctrl = entry.ctrl;
        if isempty(ctrl) || size(ctrl,2) < 13
            continue;
        end
        step = ctrl(:,1);
        qp_time_ms = 1000.0 * ctrl(:,13);

        plot(step, qp_time_ms, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});

        has_time = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('solve time (ms)', 'FontSize', Fontsize);
    title('Computation / solve time', 'FontSize', Fontsize);
    if has_time
        legend('Location', 'best', 'FontSize', 10);
    else
        text(0.1, 0.5, 'No column 13 found in data.txt', ...
            'Units', 'normalized', 'FontSize', Fontsize);
    end

    sgtitle(sprintf('Method Comparison: %s', scenario_titles{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');

    if save_figures
        out_png = fullfile(fig_dir, sprintf('fig%d_%s_comparison.png', s, sanitize_name(scenarios{s})));
        out_pdf = fullfile(fig_dir, sprintf('fig%d_%s_comparison.pdf', s, sanitize_name(scenarios{s})));
        exportgraphics(fig, out_png, 'Resolution', 300);
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
        fprintf('Saved:\n  %s\n  %s\n', out_png, out_pdf);
    end
end

%% ========================= OPTIONAL: CONTROL DECOMPOSITION FIGURES =========================
% These are useful if you want to inspect nominal vs backup vs filtered
% for each method/scenario separately. Set this to true if needed.

make_control_decomposition_figs = false;

if make_control_decomposition_figs
    fig_id = 10;

    for s = 1:numel(scenarios)
        for m = 1:numel(methods)
            entry = data(s,m).entry;
            ctrl = entry.ctrl;

            if isempty(ctrl) || size(ctrl,2) < 10
                continue;
            end

            step = ctrl(:,1);

            vx_safe = ctrl(:,2);  vy_safe = ctrl(:,3);  wz_safe = ctrl(:,4);
            vx_nom  = ctrl(:,5);  vy_nom  = ctrl(:,6);  wz_nom  = ctrl(:,7);
            vx_b    = ctrl(:,8);  vy_b    = ctrl(:,9);  wz_b    = ctrl(:,10);

            fig = figure(fig_id); fig_id = fig_id + 1;
            clf;
            set(fig, 'Color', 'w', 'Position', [150 150 1000 750]);

            tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');

            nexttile;
            plot(step, vx_safe, 'r-', 'LineWidth', 1.7); hold on;
            plot(step, vx_nom, 'k--', 'LineWidth', 1.5);
            plot(step, vx_b, 'b:', 'LineWidth', 2.5);
            grid on; ylabel('$v_x$', 'Interpreter', 'latex');
            title(sprintf('%s / %s', scenario_titles{s}, method_labels{m}), 'Interpreter', 'none');
            legend('filtered', 'nominal', 'backup', 'Location', 'best');

            nexttile;
            plot(step, vy_safe, 'r-', 'LineWidth', 1.7); hold on;
            plot(step, vy_nom, 'k--', 'LineWidth', 1.5);
            plot(step, vy_b, 'b:', 'LineWidth', 2.5);
            grid on; ylabel('$v_y$', 'Interpreter', 'latex');

            nexttile;
            plot(step, wz_safe, 'r-', 'LineWidth', 1.7); hold on;
            plot(step, wz_nom, 'k--', 'LineWidth', 1.5);
            plot(step, wz_b, 'b:', 'LineWidth', 2.5);
            grid on; ylabel('$\omega$', 'Interpreter', 'latex');
            xlabel('step');

            if save_figures
                out_png = fullfile(fig_dir, sprintf('controls_%s_%s.png', ...
                    sanitize_name(scenarios{s}), sanitize_name(method_labels{m})));
                exportgraphics(fig, out_png, 'Resolution', 300);
                fprintf('Saved:\n  %s\n', out_png);
            end
        end
    end
end

%% ========================= HELPER FUNCTIONS =========================

function [x_col, y_col] = detect_xy_columns(traj)
% Try to infer x/y columns from the trajectory CSV.
% This is intentionally defensive because the existing logs may vary.

    ncol = size(traj,2);

    % Common case from previous high-level plotting:
    % all_data(:,2) = x, all_data(:,3) = y
    if ncol >= 3
        x_col = 2;
        y_col = 3;
    else
        error('Trajectory data must have at least 3 columns.');
    end
end

function plot_obstacles(obs)
    th = linspace(0, 2*pi, 120);
    for i = 1:size(obs,1)
        ox = obs(i,1);
        oy = obs(i,2);
        rr = obs(i,3);

        xx = ox + rr*cos(th);
        yy = oy + rr*sin(th);

        fill(xx, yy, [0.35 0.35 0.35], ...
            'FaceAlpha', 0.30, ...
            'EdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
end

function s = sanitize_name(s)
    s = char(s);
    s = lower(s);
    s = strrep(s, ' ', '_');
    s = strrep(s, '-', '');
    s = strrep(s, '/', '_');
    s = strrep(s, '\', '_');
end
