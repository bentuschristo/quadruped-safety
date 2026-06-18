%% waypoint_forest_scenario_designer.m
% Editable waypoint-based smooth trajectory and circular obstacle preview.
%
% How to use:
%   1. Edit the `waypoints` array below. Rows are followed from TOP to BOTTOM.
%   2. Edit the `obstacles` array below.
%   3. Run this script.
%   4. Inspect the trajectory, heading arrows, curvature, and obstacle clearance.
%
% This script uses a parametric cubic Hermite / PCHIP curve:
%       x = x(s), y = y(s)
% where s is cumulative waypoint distance. This is safer than fitting y(x)
% because the trajectory may reverse direction in x.
%
% Obstacle columns:
%   [x_center, y_center, physical_radius, cbf_radius, backup_radius]
%
% Shading:
%   gray disk  : physical obstacle
%   red ring   : physical radius -> CBF radius
%   green ring : CBF radius -> backup radius

clear; close all; clc;

%% ================= USER SETTINGS =================

% Waypoints are followed from TOP to BOTTOM.
% Keep the first point equal to the robot's initial position.
% First Scenario
waypoints = [
     0.0,  0.0;
     0.0, 1.0;
     0.2, 2.0;
     0.0, 3.0;
     -0.5, 4.0;
     -0.7, 5.0;
     -0.3, 6.0;
     0.0, 7.0;
     0.0, 8.0;
];



% [cx, cy, physical_r, cbf_r, backup_r]
% First Scenario
obstacles = [
    0.5, 1.0, 0.30, 0.40, 0.50
    -0.5, 2.0, 0.30, 0.40, 0.50
    0.0, 4.0, 0.30, 0.40, 0.50
    0.4, 6.0, 0.30, 0.40, 0.50
    -1.0, 5.5, 0.30, 0.40, 0.50
];

% Curve type:
%   "pchip"  : shape-preserving, usually less overshoot
%   "spline" : natural cubic spline, smoother but may overshoot
curve_type = "pchip";

num_samples = 1200;

% Heading display
constant_heading = false;
constant_heading_value = 0.0;   % rad
show_heading_arrows = true;
heading_arrow_stride = 70;
heading_arrow_length = 0.25;

% Plot appearance
physical_color = [0.75, 0.75, 0.75];
cbf_color      = [1.00, 0.20, 0.20];
backup_color   = [0.65, 0.90, 0.65];
physical_alpha = 0.95;
cbf_alpha      = 0.65;
backup_alpha   = 0.45;

show_obstacle_labels = true;
show_backup_radius = true;
show_cbf_radius = true;


% Safe-set shading
show_safe_set_shading = true;   % true/false
show_safe_set_labels  = true;

% Choose either center/size or direct box limits.
use_center_size_boxes = false;

% Center/size mode: [x_center, y_center], [width, height]
CS_center = [-2.5, 0.0];
CS_size   = [6.2, 3.8];

CB_center = [-2.5, 0.0];
CB_size   = [5.2, 3.0];

% Direct box mode: [xmin, xmax, ymin, ymax]
CS_box = [-3, 3, -1, 10];
CB_box = [-2.5, 2.5, -0.9, 9];

CS_color = [0.80, 0.98, 0.82];   % light green
CB_color = [0.78, 0.86, 1.00];   % light blue
CS_alpha = 0.30;
CB_alpha = 0.35;

% Optional workspace box [xmin xmax ymin ymax].
% Set to [] to auto-fit.
plot_limits = [];

% Clearance checking uses the selected radius:
%   "physical", "cbf", or "backup"
clearance_radius = "backup";

%% =============== INPUT VALIDATION ===============

assert(size(waypoints,2) == 2, ...
    'waypoints must be N-by-2: [x, y].');
assert(size(waypoints,1) >= 2, ...
    'At least two waypoints are required.');

assert(isempty(obstacles) || size(obstacles,2) == 5, ...
    'obstacles must be M-by-5: [cx, cy, physical_r, cbf_r, backup_r].');

if ~isempty(obstacles)
    assert(all(obstacles(:,3) > 0), 'Physical radii must be positive.');
    assert(all(obstacles(:,4) >= obstacles(:,3)), ...
        'Each CBF radius must be >= physical radius.');
    assert(all(obstacles(:,5) >= obstacles(:,4)), ...
        'Each backup radius must be >= CBF radius.');
end


% Compute safe-set boxes.
if use_center_size_boxes
    CS_box_plot = [ ...
        CS_center(1)-CS_size(1)/2, ...
        CS_center(1)+CS_size(1)/2, ...
        CS_center(2)-CS_size(2)/2, ...
        CS_center(2)+CS_size(2)/2];

    CB_box_plot = [ ...
        CB_center(1)-CB_size(1)/2, ...
        CB_center(1)+CB_size(1)/2, ...
        CB_center(2)-CB_size(2)/2, ...
        CB_center(2)+CB_size(2)/2];
else
    CS_box_plot = CS_box;
    CB_box_plot = CB_box;
end

assert(CB_box_plot(1) >= CS_box_plot(1) && ...
       CB_box_plot(2) <= CS_box_plot(2) && ...
       CB_box_plot(3) >= CS_box_plot(3) && ...
       CB_box_plot(4) <= CS_box_plot(4), ...
       'C_B box must be contained inside C_S box.');

%% =============== PARAMETRIC CURVE ===============

delta_wp = diff(waypoints,1,1);
segment_lengths = hypot(delta_wp(:,1), delta_wp(:,2));

assert(all(segment_lengths > 1e-9), ...
    'Consecutive waypoints must be distinct.');

s_wp = [0; cumsum(segment_lengths)];
s_query = linspace(s_wp(1), s_wp(end), num_samples);

switch lower(curve_type)
    case "pchip"
        pp_x = pchip(s_wp, waypoints(:,1));
        pp_y = pchip(s_wp, waypoints(:,2));
    case "spline"
        pp_x = spline(s_wp, waypoints(:,1));
        pp_y = spline(s_wp, waypoints(:,2));
    otherwise
        error('curve_type must be "pchip" or "spline".');
end

x_ref = ppval(pp_x, s_query);
y_ref = ppval(pp_y, s_query);

% Numerical derivatives with respect to the path parameter s.
dx_ds = gradient(x_ref, s_query);
dy_ds = gradient(y_ref, s_query);
d2x_ds2 = gradient(dx_ds, s_query);
d2y_ds2 = gradient(dy_ds, s_query);

speed_s = hypot(dx_ds, dy_ds);

if constant_heading
    theta_ref = constant_heading_value * ones(size(s_query));
    dtheta_ds = zeros(size(s_query));
else
    theta_ref = unwrap(atan2(dy_ds, dx_ds));
    dtheta_ds = gradient(theta_ref, s_query);
end

curvature = (dx_ds .* d2y_ds2 - dy_ds .* d2x_ds2) ./ ...
    max(speed_s.^3, 1e-10);

%% =============== CLEARANCE CHECK ===============

if isempty(obstacles)
    selected_radii = [];
else
    switch lower(clearance_radius)
        case "physical"
            selected_radii = obstacles(:,3);
        case "cbf"
            selected_radii = obstacles(:,4);
        case "backup"
            selected_radii = obstacles(:,5);
        otherwise
            error('clearance_radius must be "physical", "cbf", or "backup".');
    end
end

min_clearance = inf;
closest_obstacle = NaN;
closest_path_index = NaN;
clearance_per_obstacle = nan(size(obstacles,1),1);

for j = 1:size(obstacles,1)
    distance_to_center = hypot( ...
        x_ref - obstacles(j,1), ...
        y_ref - obstacles(j,2));
    clearance_j = distance_to_center - selected_radii(j);
    [clearance_per_obstacle(j), idx_j] = min(clearance_j);

    if clearance_per_obstacle(j) < min_clearance
        min_clearance = clearance_per_obstacle(j);
        closest_obstacle = j;
        closest_path_index = idx_j;
    end
end

%% =================== FIGURE ======================

fig = figure( ...
    'Name','Waypoint Forest Scenario Designer', ...
    'Color','w', ...
    'NumberTitle','off');

tabs = uitabgroup(fig);

%% ----- Tab 1: Planar path -----
tab1 = uitab(tabs, 'Title','Trajectory');
ax1 = axes(tab1);
hold(ax1,'on');
axis(ax1,'equal');
grid(ax1,'on');
box(ax1,'on');


% Draw C_S and C_B background regions first.
if show_safe_set_shading
    draw_box_patch(ax1, CS_box_plot, CS_color, CS_alpha, 'k', 1.2);
    draw_box_patch(ax1, CB_box_plot, CB_color, CB_alpha, 'k', 1.2);

    if show_safe_set_labels
        text(ax1, CS_box_plot(2)-0.08, CS_box_plot(4)-0.08, '$\mathcal{C}_S$', ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','top', ...
            'FontSize',16, ...
            'FontWeight','bold');

        text(ax1, CB_box_plot(2)-0.08, CB_box_plot(4)-0.08, '$\mathcal{C}_B$', ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','top', ...
            'FontSize',16, ...
            'FontWeight','bold');
    end
end

% Draw obstacle regions from largest to smallest.
for j = 1:size(obstacles,1)
    cx = obstacles(j,1);
    cy = obstacles(j,2);
    r_phys = obstacles(j,3);
    r_cbf = obstacles(j,4);
    r_backup = obstacles(j,5);

    if show_backup_radius
        draw_filled_circle(ax1, cx, cy, r_backup, ...
            backup_color, backup_alpha, 'k', 1.0);
    end

    if show_cbf_radius
        draw_filled_circle(ax1, cx, cy, r_cbf, ...
            cbf_color, cbf_alpha, 'k', 1.0);
    end

    draw_filled_circle(ax1, cx, cy, r_phys, ...
        physical_color, physical_alpha, 'k', 1.0);

    if show_obstacle_labels
        text(ax1, cx, cy, sprintf('O%d',j), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontWeight','bold');
    end
end

% Curve and waypoints
plot(ax1, x_ref, y_ref, 'LineWidth', 2.2, ...
    'DisplayName','Smooth reference');
plot(ax1, waypoints(:,1), waypoints(:,2), 'o-', ...
    'LineWidth', 1.0, ...
    'MarkerSize', 7, ...
    'MarkerFaceColor','w', ...
    'DisplayName','Waypoints');

% Start and end
plot(ax1, waypoints(1,1), waypoints(1,2), 's', ...
    'MarkerSize',10,'MarkerFaceColor',[0.2 0.8 0.2], ...
    'DisplayName','Start');
plot(ax1, waypoints(end,1), waypoints(end,2), 'd', ...
    'MarkerSize',10,'MarkerFaceColor',[0.9 0.2 0.2], ...
    'DisplayName','End');

% Heading arrows
if show_heading_arrows
    arrow_idx = 1:heading_arrow_stride:num_samples;
    quiver(ax1, ...
        x_ref(arrow_idx), y_ref(arrow_idx), ...
        heading_arrow_length*cos(theta_ref(arrow_idx)), ...
        heading_arrow_length*sin(theta_ref(arrow_idx)), ...
        0, 'LineWidth', 1.0, ...
        'DisplayName','Heading');
end

% Highlight closest point
if isfinite(min_clearance)
    plot(ax1, ...
        x_ref(closest_path_index), y_ref(closest_path_index), ...
        'p', 'MarkerSize',12, ...
        'MarkerFaceColor',[1.0 0.8 0.0], ...
        'DisplayName','Minimum clearance');
end

xlabel(ax1,'X [m]');
ylabel(ax1,'Y [m]');
title(ax1, sprintf('%s path | min %s clearance = %.3f m', ...
    upper(curve_type), clearance_radius, min_clearance));
legend(ax1,'Location','best');

if ~isempty(plot_limits)
    xlim(ax1, plot_limits(1:2));
    ylim(ax1, plot_limits(3:4));
else
    all_x = [x_ref(:); waypoints(:,1)];
    all_y = [y_ref(:); waypoints(:,2)];

    if show_safe_set_shading
        all_x = [all_x; CS_box_plot(1); CS_box_plot(2); ...
                       CB_box_plot(1); CB_box_plot(2)];
        all_y = [all_y; CS_box_plot(3); CS_box_plot(4); ...
                       CB_box_plot(3); CB_box_plot(4)];
    end

    if ~isempty(obstacles)
        all_x = [all_x; obstacles(:,1)-obstacles(:,5); ...
                       obstacles(:,1)+obstacles(:,5)];
        all_y = [all_y; obstacles(:,2)-obstacles(:,5); ...
                       obstacles(:,2)+obstacles(:,5)];
    end

    margin = 0.4;
    xlim(ax1,[min(all_x)-margin,max(all_x)+margin]);
    ylim(ax1,[min(all_y)-margin,max(all_y)+margin]);
end

%% ----- Tab 2: Heading -----
tab2 = uitab(tabs,'Title','Heading');
ax2 = axes(tab2);
plot(ax2,s_query,theta_ref,'LineWidth',1.8);
grid(ax2,'on');
xlabel(ax2,'Path parameter s [m]');
ylabel(ax2,'\theta_{ref} [rad]');
title(ax2, ternary(constant_heading, ...
    'Constant reference heading', ...
    'Tangent reference heading'));

%% ----- Tab 3: Curvature -----
tab3 = uitab(tabs,'Title','Curvature');
ax3 = axes(tab3);
plot(ax3,s_query,curvature,'LineWidth',1.8);
grid(ax3,'on');
xlabel(ax3,'Path parameter s [m]');
ylabel(ax3,'Curvature [1/m]');
title(ax3,'Reference-path curvature');

%% ----- Tab 4: Clearance -----
tab4 = uitab(tabs,'Title','Clearance');
ax4 = axes(tab4);
hold(ax4,'on');
grid(ax4,'on');

for j = 1:size(obstacles,1)
    distance_to_center = hypot( ...
        x_ref - obstacles(j,1), ...
        y_ref - obstacles(j,2));
    plot(ax4,s_query,distance_to_center-selected_radii(j), ...
        'LineWidth',1.2, ...
        'DisplayName',sprintf('Obstacle %d',j));
end

yline(ax4,0,'--','Collision / boundary');
xlabel(ax4,'Path parameter s [m]');
ylabel(ax4,sprintf('Clearance to %s radius [m]',clearance_radius));
title(ax4,'Trajectory clearance');
if ~isempty(obstacles)
    legend(ax4,'Location','best');
end

%% ================= SUMMARY =======================

fprintf('\nWaypoint forest scenario summary\n');
fprintf('--------------------------------\n');
fprintf('Curve type: %s\n',curve_type);
fprintf('Waypoints: %d\n',size(waypoints,1));
fprintf('Obstacles: %d\n',size(obstacles,1));
fprintf('Path parameter length: %.3f m\n',s_wp(end));

if isempty(obstacles)
    fprintf('No obstacles defined.\n');
elseif min_clearance >= 0
    fprintf('PASS: minimum %s clearance = %.4f m at obstacle %d.\n', ...
        clearance_radius,min_clearance,closest_obstacle);
else
    fprintf(2, ...
        'WARNING: trajectory violates %s radius by %.4f m at obstacle %d.\n', ...
        clearance_radius,-min_clearance,closest_obstacle);
end

%% =============== LOCAL FUNCTIONS =================


function draw_box_patch(ax, box_values, face_color, face_alpha, edge_color, line_width)
    xmin = box_values(1);
    xmax = box_values(2);
    ymin = box_values(3);
    ymax = box_values(4);

    patch(ax, ...
        [xmin xmax xmax xmin], ...
        [ymin ymin ymax ymax], ...
        face_color, ...
        'FaceAlpha', face_alpha, ...
        'EdgeColor', edge_color, ...
        'LineWidth', line_width, ...
        'HandleVisibility','off');
end

function draw_filled_circle(ax,cx,cy,r,face_color,face_alpha,edge_color,line_width)
    ang = linspace(0,2*pi,240);
    x = cx + r*cos(ang);
    y = cy + r*sin(ang);

    patch(ax,x,y,face_color, ...
        'FaceAlpha',face_alpha, ...
        'EdgeColor',edge_color, ...
        'LineWidth',line_width, ...
        'HandleVisibility','off');
end

function value = ternary(condition,true_value,false_value)
    if condition
        value = true_value;
    else
        value = false_value;
    end
end
