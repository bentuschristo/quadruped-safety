%% compare_backup_methods_update_v2.m
%   col 1: step
%   col 2: x_ref
%   col 3: y_ref
%   col 4: theta_ref
%   col 5: x_p
%   col 6: y_p
%   col 7: theta
%   col 8: error_x
%   col 9: error_y
%   col 10: error_theta
%   col 11: cmd1
%   col 12: cmd2
%   col 13: cmd3
%
% This script therefore loads data.m from each method/scenario folder.

clear; clc; close all;

%% ========================= USER SETTINGS =========================

% Put this script inside:
%   /results/Backup CBFs/update/
% and leave root_dir = pwd.
root_dir = pwd;

% Or use an absolute path, for example:
% root_dir = '/home/bgerodac/quadruped/a1_waypoint/sample/results/Backup CBFs/update';

save_figures = false;  % set true after plots look correct

fig_dir = fullfile(root_dir, 'figures_compare_methods');
if save_figures && ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

methods = { ...
    'Backup CBF - Baseline', ...
    'Backup CBF - LSE', ...
    'Backup CBF - LSE - CF', ...
    'Blended Inputs', ...
    'Optimal Interpolation', ...
    'Optimal Interpolation - QP', ...
};

method_labels = { ...
    'bCBF-QP', ...
    'bCBF-QP-LSE', ...
    'Backup CBF - LSE - CF', ...
    'Blended Inputs', ...
    'Optimal Interpolation', ...
    'Optimal Interpolation - QP', ...
};

scenarios = { ...
    'Backup in Y', ...
    'Zero-Control' ...
};

scenario_titles = { ...
    'Backup in Y', ...
    'Zero-Control Backup' ...
};

Linewidth = 2.0;
Fontsize = 12;

colors = { ...
    [0.00 0.00 0.00], ...   % bCBF black
    [0.85 0.33 0.10], ...   % bCBF-LSE orange
    [0.90 0.50 0.70], ...   % bCBF-LSE-CF
    [0.00 0.45 0.74], ...   % blended blue
    [0.47 0.67 0.19] ...    % OI green
    [0.47 0.67 0.60] ...    % OI - QP
};

linestyles = {'-', ':', '--','--', '-.', '--'};

% Obstacles. Keep consistent with your MuJoCo XML / C filter.
obs = [ ...
    -2.0,  1.0, 0.35;
    -4.0, -1.0, 0.35;
    -2.0, -1.0, 0.35;
    -4.0,  1.0, 0.35
];

% Uncomment if your run includes the two middle obstacles:
% obs = [obs;
%     -2.8,  0.0, 0.30;
%     -3.3,  0.0, 0.30];

draw_obstacles = true;

%% ========================= LOAD DATA =========================

D = struct();

for s = 1:numel(scenarios)
    for m = 1:numel(methods)
        folder = fullfile(root_dir, methods{m}, scenarios{s});

        D(s,m).method = methods{m};
        D(s,m).method_label = method_labels{m};
        D(s,m).scenario = scenarios{s};
        D(s,m).folder = folder;
        D(s,m).all_data = [];
        D(s,m).ctrl = [];

        [all_data, ok_data] = load_highlevel_data_m(folder);
        if ok_data
            D(s,m).all_data = all_data;
        else
            warning('Could not load high-level data.m in: %s', folder);
        end

        txt_file = fullfile(folder, 'data.txt');
        if exist(txt_file, 'file') == 2
            D(s,m).ctrl = readmatrix(txt_file);
        else
            warning('Missing controller log: %s', txt_file);
        end
    end
end

%% ========================= FIGURE 1 and FIGURE 2 =========================

for s = 1:numel(scenarios)

    fig = figure(s);
    clf;
    set(fig, 'Color', 'w', 'Position', [80 80 1250 850]);

    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    %% 1) Trajectory
    nexttile(1);
    hold on; grid on; box on;

    ref_plotted = false;

    for m = 1:numel(methods)
        all_data = D(s,m).all_data;
        if isempty(all_data) || size(all_data,2) < 6
            continue;
        end

        step = all_data(:,1);
        x_ref = all_data(:,2);
        y_ref = all_data(:,3);
        x_p   = all_data(:,5);
        y_p   = all_data(:,6);

        plot(x_p, y_p, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});

        if ~ref_plotted
            plot(x_ref, y_ref, 'k:', 'LineWidth', 1.2, 'DisplayName', 'reference');
            ref_plotted = true;
        end
    end

    if draw_obstacles
        plot_obstacles(obs);
    end

    axis equal;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title(sprintf('%s: trajectory comparison', scenario_titles{s}), ...
        'Interpreter', 'none', 'FontSize', Fontsize, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);

    %% 2) Filtered vx
    nexttile(2);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        plot(ctrl(:,1), ctrl(:,2), ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$v_x^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $v_x$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    legend('Location', 'best', 'FontSize', 9);

    %% 3) Filtered vy
    nexttile(3);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        plot(ctrl(:,1), ctrl(:,3), ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$v_y^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $v_y$', 'Interpreter', 'latex', 'FontSize', Fontsize);

    %% 4) Filtered omega
    nexttile(4);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 4
            continue;
        end
        plot(ctrl(:,1), ctrl(:,4), ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\omega^{safe}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Filtered $\omega$', 'Interpreter', 'latex', 'FontSize', Fontsize);

    %% 5) mu
    nexttile(5);
    hold on; grid on; box on;
    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 11
            continue;
        end
        mu = ctrl(:,11);
        mu(mu < 0) = NaN; % bCBF has mu = -1
        plot(ctrl(:,1), mu, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});
    end
    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\mu$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title('Interpolation / blending parameter', 'Interpreter', 'latex', 'FontSize', Fontsize);
    ylim([-0.05, 1.05]);
    legend('Location', 'best', 'FontSize', 9);

    %% 6) computation time
    nexttile(6);
    hold on; grid on; box on;
    has_time = false;

    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 13
            continue;
        end

        qp_time_ms = 1000 * ctrl(:,13);

        plot(ctrl(:,1), qp_time_ms, ...
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
        legend('Location', 'best', 'FontSize', 9);
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


%% ========================= TIMING-ONLY FIGURES =========================
% Figure 3: Backup in Y timing comparison
% Figure 4: Zero-Control timing comparison
%
% Requires data.txt column 13:
%   qp_time_sec

for s = 1:numel(scenarios)

    fig_timing = figure(2 + s);
    clf;
    set(fig_timing, 'Color', 'w', 'Position', [120 120 1100 500]);

    hold on; grid on; box on;

    has_time = false;

    for m = 1:numel(methods)

        ctrl = D(s,m).ctrl;

        if isempty(ctrl) || size(ctrl,2) < 13
            warning('No timing column for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        step = ctrl(:,1);

        % Convert seconds to milliseconds
        solve_time_ms = 1000.0 * ctrl(:,13);

        % Remove NaN/Inf before computing statistics
        valid_time = solve_time_ms(isfinite(solve_time_ms));

        if isempty(valid_time)
            warning('No valid timing data for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        avg_time = mean(valid_time);
        max_time = max(valid_time);

        legend_label = sprintf('%s, avg = %.4g ms', ...
            method_labels{m}, avg_time);

        % If you also want max in the legend, use this instead:
        % legend_label = sprintf('%s, avg = %.4g ms, max = %.4g ms', ...
        %     method_labels{m}, avg_time, max_time);

        plot(step, solve_time_ms, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', legend_label);

        has_time = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('solve time (ms)', 'FontSize', Fontsize);

    title(sprintf('Computation Time Comparison: %s', scenario_titles{s}), ...
        'FontSize', Fontsize + 2, ...
        'FontWeight', 'bold');

    if has_time
        legend('Location', 'best', 'FontSize', 10);
    else
        text(0.1, 0.5, ...
            'No timing data found. Need column 13 in data.txt.', ...
            'Units', 'normalized', ...
            'FontSize', Fontsize);
    end

    if save_figures
        out_png = fullfile(fig_dir, ...
            sprintf('timing_%s.png', sanitize_name(scenarios{s})));

        out_pdf = fullfile(fig_dir, ...
            sprintf('timing_%s.pdf', sanitize_name(scenarios{s})));

        exportgraphics(fig_timing, out_png, 'Resolution', 300);
        exportgraphics(fig_timing, out_pdf, 'ContentType', 'vector');

        fprintf('Saved timing figure:\n  %s\n  %s\n', out_png, out_pdf);
    end
end

%% ========================= HELPER FUNCTIONS =========================

function [all_data, ok] = load_highlevel_data_m(folder)
% Load high-level all_data by running data.m in the target folder.
% This matches the original main_highlevel_v2.m behavior.
    all_data = [];
    ok = false;

    data_m = fullfile(folder, 'data.m');
    if exist(data_m, 'file') ~= 2
        return;
    end

    old_dir = pwd;

    try
        cd(folder);
        clear all_data; %#ok<CLALL>
        run('data.m');

        if exist('all_data', 'var')
            ok = true;
        end

        cd(old_dir);
    catch ME
        cd(old_dir);
        warning('Failed to run data.m in %s\nReason: %s', folder, ME.message);
        return;
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

        text(ox, oy, sprintf('O%d', i), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
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
