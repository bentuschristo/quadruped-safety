%% compare_backup_methods_update_v6_shaded.m
% One-plot-per-figure version.
%
% Expected high-level data.m columns:
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
% Expected controller log data.txt columns:
%   col 1: step
%   col 2: vx_safe
%   col 3: vy_safe
%   col 4: omega_safe
%   ...
%   col 11: mu or relaxation variable depending on method
%   col 13: solver/core time in seconds
%   col 14: full safety-filter time in seconds (optional, newer logs)
%
% Notes:
% - For Blending/OI methods, column 11 is interpolation parameter mu in [0,1].
% - For the closed-form bCBF-control-dynamics methods, column 11 is NOT the
%   interpolation mu. It is the relaxation variable, so this script does not
%   plot it in the "Interpolation parameter" figure by default.

clear; clc; close all;

%% ========================= USER SETTINGS =========================

root_dir = pwd;

save_figures = false;  % set true after plots look correct

fig_dir = fullfile(root_dir, 'figures_compare_methods_v6_shaded');
if save_figures && ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

methods = { ...
    'bCBF-QP', ...
    'bCBF-LSE-QP', ...
    'bCBF-LSE-CF', ...
    'Blended', ...
    'OI-CF', ...
    'OI-QP' ...
};

method_labels = { ...
    'bCBF-QP', ...
    'bCBF-LSE-QP', ...
    'bCBF-LSE-CF', ...
    'Blended Inputs', ...
    'OI-CF', ...
    'OI-QP' ...
};

% Which methods have true blending/interpolation mu in column 11?
% Closed-form bCBF control-dynamics methods use relaxation mu, not interpolation mu.
is_interpolation_mu = [ ...
    false, ... % Backup CBF - Baseline
    false, ... % Backup CBF - LSE
    false, ... % Backup CBF - LSE
    true,  ... % Blended Inputs
    true,  ... % OI-CF
    true   ... % OI-QP
];

scenarios = { ...
    'Zero-Control' ...
};

scenario_titles = { ...
    'Zero-Control Backup' ...
};

Linewidth = 2.0;
Fontsize = 12;

colors = { ...
    [0.00 0.00 0.00], ...   % bCBF black
    [0.85 0.33 0.10], ...   % bCBF-LSE orange
    [0.90 0.50 0.70], ...   % old CF pink
    [0.00 0.45 0.74], ...   % blended blue
    [0.47 0.67 0.19], ...   % OI green
    [0.47 0.67 0.60] ...    % OI-QP teal
};

linestyles = {'-', ':', '--', '--', '-.', '--'};

%% ----- Figure 1 region / obstacle shading settings -----
% Obstacles for the trajectory plot. Columns:
%   [x_center, y_center, real_radius, CS_barrier_radius, CB_backup_radius]
%
% Interpretation:
%   gray disk     = real obstacle, radius real_radius
%   red annulus   = buffer / C_S obstacle barrier, real_radius -> CS_barrier_radius
%   green annulus = C_S \ C_B near obstacle, CS_barrier_radius -> CB_backup_radius
%
% Because C_B is a subset of C_S, CB_backup_radius should be larger than
% CS_barrier_radius. The blue C_B region is therefore removed around a larger
% circle than the red C_S barrier circle.
obs = [ ...
    -2.0,  1.0, 0.35, 0.45, 0.55;
    -4.0, -1.0, 0.35, 0.45, 0.55;
    -2.0, -1.0, 0.35, 0.45, 0.55;
    -4.0,  1.0, 0.35, 0.45, 0.55
];

% Optional background set shading, similar to the diagram in the slides.
draw_region_shading = true;
draw_obstacles = true;

% Define C_S and C_B rectangles by center and size.
% Size is [width, height]. This is usually easier than editing xmin/xmax.
CS_center = [-3.0, 0.0];
CS_size   = [4.4, 4.6];     % [width, height]
CB_center = [-3.0, 0.0];
CB_size   = [3.4, 3.6];     % [width, height]

% If you prefer explicit boxes, set use_center_size_boxes=false and edit:
%   CS_box = [xmin xmax ymin ymax];
%   CB_box = [xmin xmax ymin ymax];
use_center_size_boxes = false;
CS_box = [-7, 1, -2.3, 2.3];
CB_box = [-7, 1, -1.8, 1.8];

if use_center_size_boxes
    CS_box = center_size_to_box(CS_center, CS_size);
    CB_box = center_size_to_box(CB_center, CB_size);
end

% Shading colors and transparency. Change freely.
CS_color       = [0.50 1.0 0.75];   % light green, safe set C_S
CB_color       = [0.75 0.50 1.00];   % light blue, backup set C_B
obstacle_color = [0.85 0.85 0.85];   % gray real obstacle
buffer_color   = [1.00 0.20 0.20];   % red buffer between obstacle and C_S limit
CS_only_color  = CS_color;           % green C_S \ C_B obstacle annulus
CS_alpha       = 0.35;
CB_alpha       = 0.35;
obstacle_alpha = 0.85;
buffer_alpha   = 0.65;
CS_only_alpha  = 0.55;

show_region_labels = true;
show_obstacle_labels = true;

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

%% ========================= ONE FIGURE PER PLOT =========================

fig_id = 1;

for s = 1:numel(scenarios)

    %% 1) Trajectory comparison
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [80 80 900 500]);
    hold on; grid on; box on;

    ref_plotted = false;

    if draw_region_shading
        plot_region_shading(CS_box, CB_box, CS_color, CB_color, CS_alpha, CB_alpha, show_region_labels);
    end
    if draw_obstacles
        plot_obstacle_shading(obs, obstacle_color, buffer_color, CS_only_color, obstacle_alpha, buffer_alpha, CS_only_alpha, show_obstacle_labels);
    end

    for m = 1:numel(methods)
        all_data = D(s,m).all_data;
        if isempty(all_data) || size(all_data,2) < 6
            continue;
        end

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

    axis equal;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title(sprintf('%s: trajectory comparison', scenario_titles{s}), ...
        'Interpreter', 'none', 'FontSize', Fontsize, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);

    save_if_requested(fig, fig_dir, save_figures, sprintf('trajectory_%s', sanitize_name(scenarios{s})));

    %% 2) Filtered vx
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [100 100 900 420]);
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
    title(sprintf('%s: filtered $v_x$', scenario_titles{s}), ...
        'Interpreter', 'latex', 'FontSize', Fontsize, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);

    save_if_requested(fig, fig_dir, save_figures, sprintf('vx_%s', sanitize_name(scenarios{s})));

    %% 3) Filtered vy
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [120 120 900 420]);
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
    title(sprintf('%s: filtered $v_y$', scenario_titles{s}), ...
        'Interpreter', 'latex', 'FontSize', Fontsize, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);

    save_if_requested(fig, fig_dir, save_figures, sprintf('vy_%s', sanitize_name(scenarios{s})));

    %% 4) Filtered omega
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [140 140 900 420]);
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
    title(sprintf('%s: filtered $\\omega$', scenario_titles{s}), ...
        'Interpreter', 'latex', 'FontSize', Fontsize, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);

    save_if_requested(fig, fig_dir, save_figures, sprintf('omega_%s', sanitize_name(scenarios{s})));

    %% 5) Interpolation / blending mu only
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [160 160 900 420]);
    hold on; grid on; box on;

    has_mu = false;

    for m = 1:numel(methods)
        if ~is_interpolation_mu(m)
            continue;
        end

        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 11
            continue;
        end

        mu = ctrl(:,11);
        mu(mu < 0) = NaN;

        plot(ctrl(:,1), mu, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});

        has_mu = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\mu_{\mathrm{interp}}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title(sprintf('%s: interpolation / blending parameter', scenario_titles{s}), ...
        'Interpreter', 'none', 'FontSize', Fontsize, 'FontWeight', 'bold');
    ylim([-0.05, 1.05]);

    if has_mu
        legend('Location', 'best', 'FontSize', 9);
    else
        text(0.1, 0.5, 'No interpolation/blending \mu data found', ...
            'Units', 'normalized', 'FontSize', Fontsize);
    end

    save_if_requested(fig, fig_dir, save_figures, sprintf('mu_interp_%s', sanitize_name(scenarios{s})));

    %% 6) Relaxation mu for closed-form bCBF only
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [180 180 900 420]);
    hold on; grid on; box on;

    has_relax = false;

    for m = 1:numel(methods)
        if is_interpolation_mu(m)
            continue;
        end

        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 11
            continue;
        end

        mu_relax = ctrl(:,11);

        % Baseline bCBF often logs -1. Skip constant negative placeholders.
        if all(mu_relax < 0 | isnan(mu_relax))
            continue;
        end

        plot(ctrl(:,1), mu_relax, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', method_labels{m});

        has_relax = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('$\mu_{\mathrm{relax}}$', 'Interpreter', 'latex', 'FontSize', Fontsize);
    title(sprintf('%s: closed-form CBF relaxation variable', scenario_titles{s}), ...
        'Interpreter', 'none', 'FontSize', Fontsize, 'FontWeight', 'bold');

    if has_relax
        legend('Location', 'best', 'FontSize', 9);
    else
        text(0.1, 0.5, 'No non-placeholder relaxation \mu data found', ...
            'Units', 'normalized', 'FontSize', Fontsize);
    end

    save_if_requested(fig, fig_dir, save_figures, sprintf('mu_relax_%s', sanitize_name(scenarios{s})));

    %% 7) Computation / solver time
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [200 200 1000 460]);
    hold on; grid on; box on;

    has_time = false;

    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl) || size(ctrl,2) < 13
            warning('No timing column for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        step = ctrl(:,1);
        solve_time_ms = 1000.0 * ctrl(:,13);
        valid_time = solve_time_ms(isfinite(solve_time_ms));

        if isempty(valid_time)
            warning('No valid timing data for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        avg_time = mean(valid_time);
        max_time = max(valid_time);

        legend_label = sprintf('%s, avg %.4g ms, max %.4g ms', ...
            method_labels{m}, avg_time, max_time);

        plot(step, solve_time_ms, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', legend_label);

        has_time = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('solver / core time (ms)', 'FontSize', Fontsize);
    title(sprintf('%s: solver / core time comparison', scenario_titles{s}), ...
        'FontSize', Fontsize + 2, 'FontWeight', 'bold');

    if has_time
        legend('Location', 'best', 'FontSize', 9);
    else
        text(0.1, 0.5, 'No solver timing data found. Need column 13 in data.txt.', ...
            'Units', 'normalized', 'FontSize', Fontsize);
    end

    save_if_requested(fig, fig_dir, save_figures, sprintf('solver_timing_%s', sanitize_name(scenarios{s})));

    %% 8) Full safety-filter time (new logs, column 14)
    fig = figure(fig_id); fig_id = fig_id + 1;
    clf;
    set(fig, 'Color', 'w', 'Position', [220 220 1000 460]);
    hold on; grid on; box on;

    has_full_time = false;

    for m = 1:numel(methods)
        ctrl = D(s,m).ctrl;
        if isempty(ctrl)
            continue;
        end

        step = ctrl(:,1);

        if size(ctrl,2) >= 14
            full_filter_time_ms = 1000.0 * ctrl(:,14);
            time_source = 'full';
        elseif size(ctrl,2) >= 13
            % Backward compatibility for older logs: only solver/core time exists.
            full_filter_time_ms = 1000.0 * ctrl(:,13);
            time_source = 'fallback col13';
        else
            warning('No full-filter timing column for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        valid_time = full_filter_time_ms(isfinite(full_filter_time_ms));

        if isempty(valid_time)
            warning('No valid full-filter timing data for %s / %s', method_labels{m}, scenarios{s});
            continue;
        end

        avg_time = mean(valid_time);
        max_time = max(valid_time);

        if strcmp(time_source, 'full')
            legend_label = sprintf('%s, avg %.4g ms, max %.4g ms', ...
                method_labels{m}, avg_time, max_time);
        else
            legend_label = sprintf('%s, avg %.4g ms, max %.4g ms (%s)', ...
                method_labels{m}, avg_time, max_time, time_source);
        end

        plot(step, full_filter_time_ms, ...
            'LineWidth', Linewidth, ...
            'Color', colors{m}, ...
            'LineStyle', linestyles{m}, ...
            'DisplayName', legend_label);

        has_full_time = true;
    end

    xlabel('step', 'FontSize', Fontsize);
    ylabel('full safety-filter time (ms)', 'FontSize', Fontsize);
    title(sprintf('%s: full safety-filter time comparison', scenario_titles{s}), ...
        'FontSize', Fontsize + 2, 'FontWeight', 'bold');

    if has_full_time
        legend('Location', 'best', 'FontSize', 9);
    else
        text(0.1, 0.5, 'No full-filter timing data found. Need column 14 in data.txt.', ...
            'Units', 'normalized', 'FontSize', Fontsize);
    end

    save_if_requested(fig, fig_dir, save_figures, sprintf('full_filter_timing_%s', sanitize_name(scenarios{s})));
end

%% ========================= HELPER FUNCTIONS =========================

function [all_data, ok] = load_highlevel_data_m(folder)
% Load high-level all_data by running data.m in the target folder.
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

function plot_region_shading(CS_box, CB_box, CS_color, CB_color, CS_alpha, CB_alpha, show_labels)
% Draw C_S as the larger background rectangle and C_B as the inner rectangle.
% Boxes are [xmin xmax ymin ymax]. Patches are hidden from legend.
    draw_box(CS_box, CS_color, CS_alpha, 'off');
    draw_box(CB_box, CB_color, CB_alpha, 'off');

    % Boundary outlines for clarity.
    rectangle('Position', [CS_box(1), CS_box(3), CS_box(2)-CS_box(1), CS_box(4)-CS_box(3)], ...
        'EdgeColor', [0 0 0], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    rectangle('Position', [CB_box(1), CB_box(3), CB_box(2)-CB_box(1), CB_box(4)-CB_box(3)], ...
        'EdgeColor', [0 0 0], 'LineWidth', 1.2, 'HandleVisibility', 'off');

    if show_labels
        text(CS_box(2), CS_box(4), '$C_S$', ...
            'Interpreter', 'latex', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 18, 'FontWeight', 'bold', ...
            'HandleVisibility', 'off');
        text(CB_box(2), CB_box(4), '$C_B$', ...
            'Interpreter', 'latex', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 16, 'FontWeight', 'bold', ...
            'HandleVisibility', 'off');
    end
end

function box = center_size_to_box(center_xy, size_wh)
% Convert [center_x, center_y] and [width, height] into [xmin xmax ymin ymax].
    cx = center_xy(1); cy = center_xy(2);
    w = size_wh(1); h = size_wh(2);
    box = [cx - w/2, cx + w/2, cy - h/2, cy + h/2];
end

function draw_box(box, color, alpha_val, visibility)
    x1 = box(1); x2 = box(2); y1 = box(3); y2 = box(4);
    patch([x1 x2 x2 x1], [y1 y1 y2 y2], color, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none', ...
        'HandleVisibility', visibility);
end

function plot_obstacle_shading(obs, obstacle_color, buffer_color, CS_only_color, obstacle_alpha, buffer_alpha, CS_only_alpha, show_labels)
% Draw obstacle regions for the set diagram.
% obs columns:
%   [x_center, y_center, real_radius, CS_barrier_radius, CB_backup_radius]
%
% If only four columns are supplied, CB_backup_radius is set equal to
% CS_barrier_radius. If only three columns are supplied, all radii are equal.
    th = linspace(0, 2*pi, 240);
    for i = 1:size(obs,1)
        ox = obs(i,1);
        oy = obs(i,2);
        r_real = obs(i,3);

        if size(obs,2) >= 4
            r_CS = obs(i,4);
        else
            r_CS = r_real;
        end

        if size(obs,2) >= 5
            r_CB = obs(i,5);
        else
            r_CB = r_CS;
        end

        % Enforce ordered radii for drawing.
        if r_CS < r_real
            warning('Obstacle %d has CS barrier radius smaller than real radius. Setting r_CS = r_real.', i);
            r_CS = r_real;
        end
        if r_CB < r_CS
            warning('Obstacle %d has CB backup radius smaller than CS barrier radius. Setting r_CB = r_CS.', i);
            r_CB = r_CS;
        end

        % C_S \ C_B annulus near obstacle: backup set excludes a larger radius.
        % Draw this first so it cuts a green hole through the blue C_B shading.
        if r_CB > r_CS
            x_outer = ox + r_CB*cos(th);
            y_outer = oy + r_CB*sin(th);
            x_inner = ox + r_CS*cos(fliplr(th));
            y_inner = oy + r_CS*sin(fliplr(th));
            patch([x_outer x_inner], [y_outer y_inner], CS_only_color, ...
                'FaceAlpha', CS_only_alpha, ...
                'EdgeColor', [0.10 0.25 0.10], ...
                'LineWidth', 0.9, ...
                'HandleVisibility', 'off');
        end

        % Red buffer zone: annulus between real obstacle and C_S barrier radius.
        if r_CS > r_real
            x_outer = ox + r_CS*cos(th);
            y_outer = oy + r_CS*sin(th);
            x_inner = ox + r_real*cos(fliplr(th));
            y_inner = oy + r_real*sin(fliplr(th));
            patch([x_outer x_inner], [y_outer y_inner], buffer_color, ...
                'FaceAlpha', buffer_alpha, ...
                'EdgeColor', [0.20 0.00 0.00], ...
                'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
        end

        % Real obstacle disk.
        x_real = ox + r_real*cos(th);
        y_real = oy + r_real*sin(th);
        fill(x_real, y_real, obstacle_color, ...
            'FaceAlpha', obstacle_alpha, ...
            'EdgeColor', [0.10 0.10 0.10], ...
            'LineWidth', 1.0, ...
            'HandleVisibility', 'off');

        % Boundary outlines.
        if r_CS > r_real
            plot(ox + r_CS*cos(th), oy + r_CS*sin(th), ...
                'Color', [0.20 0.00 0.00], 'LineWidth', 0.9, ...
                'HandleVisibility', 'off');
        end
        if r_CB > r_CS
            plot(ox + r_CB*cos(th), oy + r_CB*sin(th), ...
                'Color', [0.10 0.10 0.10], 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
        end

        if show_labels
            text(ox, oy, sprintf('O%d', i), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'HandleVisibility', 'off');
        end
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

function save_if_requested(fig, fig_dir, save_figures, base_name)
    if ~save_figures
        return;
    end

    out_png = fullfile(fig_dir, [base_name, '.png']);
    out_pdf = fullfile(fig_dir, [base_name, '.pdf']);

    exportgraphics(fig, out_png, 'Resolution', 300);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');

    fprintf('Saved:\n  %s\n  %s\n', out_png, out_pdf);
end
