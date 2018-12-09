function spider_plot(P, P_labels, axes_interval, varargin)
% Create a spider web or radar plot
%
% spider_plot(P, P_labels, axes_interval) creates a spider web plot using
% the points specified in the array P. The column of P contains the data 
% points and the rows of P contain the multiple sets of data points.
% Each point must be accompanied by a label specified in the cell
% P_labels. The number of intervals that separate the axes is specified
% by axes_interval.
% 
% P - [vector | matrix]
% P_labels - [cell of strings]
% axes_interval - [integer]
%
% spider_plot(P, P_labels, axes_interval, line_spec) works the same as
% the function above. Additional line properties can be added in the 
% same format as the default plot function in MATLAB.
%
% line_spec - [character vector]
%
%
% %%%%%%%%%%%%%%%%%%% Example of a Generic Spider Plot %%%%%%%%%%%%%%%%%%%
% % Clear workspace
% close all;
% clearvars;
% clc;
% 
% % Point properties
% num_of_points = 6;
% row_of_points = 4;
% P = rand(row_of_points, num_of_points);
% 
% % Create generic labels
% P_labels = cell(num_of_points, 1);
% 
% for ii = 1:num_of_points
%     P_labels{ii} = sprintf('Label %i', ii);
% end
% 
% % Figure properties
% figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);
% 
% % Axes interval
% axes_interval = 4;
% 
% % Spider plot
% spider_plot(P, P_labels, axes_interval,...
%     'Marker', 'o',...
%     'LineStyle', '-',...
%     'LineWidth', 2,...
%     'MarkerSize', 5);
% 
% % Title properties
% title('Sample Spider Plot',...
%     'Fontweight', 'bold',...
%     'FontSize', 12);
% 
% % Legend properties
% legend('show', 'Location', 'southoutside');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Point Properties %%%
% Number of points
[num_of_rows, num_points] = size(P);

%%% Polar Axes %%%
% Polar increments
polar_increments = 2*pi/num_points;

% Calculate max limit of axes
max_value = max(max(P));
axes_limit = ceil(max_value);

% Polar points
radius = [0; axes_limit];
theta = 0:polar_increments:2*pi;

% Convert polar to cartesian coordinates
[x_axes, y_axes] = pol2cartvect(theta, radius);

% Plot polar axes
grey = [1, 1, 1] * 0.5;
h = line(x_axes, y_axes,...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the line handles
for ii = 1:length(h)
    % Remove polar axes from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Polar Isocurves %%%
% Incremental radius
radius = (0:axes_limit/axes_interval:axes_limit)';

% Convert polar to cartesian coordinates
[x_isocurves, y_isocurves] = pol2cartvect(theta, radius);

% Plot polar isocurves
hold on;
h = plot(x_isocurves', y_isocurves',...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the plot handles
for ii = 1:length(h)
    % Remove polar isocurves from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Data Points %%%
% Iterate through all the rows
for ii = 1:num_of_rows
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cartvect(theta(1:end-1), P(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    % Plot data points
    plot(x_circular, y_circular, varargin{:});
end

%%% Axis Properties %%%
% Figure background
fig = gcf;
fig.Color = 'white';

% Iterate through the isocurve radius
for ii = 1:length(radius)
    % Display axis text for each isocurve
    text(x_isocurves(ii, 1), 0, sprintf('%3.1f', radius(ii)),...
        'Units', 'Data',...
        'Color', 'black',...
        'FontSize', 10,...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle');
end

% Label points
x_label = x_isocurves(end, :);
y_label = y_isocurves(end, :);

% Check if labels equals the number of points
if length(P_labels) == num_points
    % Iterate through each label
    for ii = 1:num_points
        % Angle of point in radians
        theta_point = theta(ii);
        
        % Find out which quadrant the point is in
        if theta_point == 0
            quadrant = 0;
        elseif theta_point == pi/2
            quadrant = 1.5;
        elseif theta_point == pi
            quadrant = 2.5;
        elseif theta_point == 3*pi/2
            quadrant = 3.5;
        elseif theta_point == 2*pi
            quadrant = 0;
        elseif theta_point > 0 && theta_point < pi/2
            quadrant = 1;
        elseif theta_point > pi/2 && theta_point < pi
            quadrant = 2;
        elseif theta_point > pi && theta_point < 3*pi/2
            quadrant = 3;
        elseif theta_point > 3*pi/2 && theta_point < 2*pi
            quadrant = 4;
        end
        
        % Adjust text alignment information depending on quadrant
        switch quadrant
            case 0
                horz_align = 'left';
                vert_align = 'middle';
            case 1
                horz_align = 'left';
                vert_align = 'bottom';
            case 1.5
                horz_align = 'center';
                vert_align = 'bottom';
            case 2
                horz_align = 'right';
                vert_align = 'bottom';
            case 2.5
                horz_align = 'right';
                vert_align = 'middle';
            case 3
                horz_align = 'right';
                vert_align = 'top';
            case 3.5
                horz_align = 'center';
                vert_align = 'top';
            case 4
                horz_align = 'left';
                vert_align = 'top';
        end
        
        % Display text label
        text(x_label(ii), y_label(ii), P_labels{ii},...
            'Units', 'Data',...
            'HorizontalAlignment', horz_align,...
            'VerticalAlignment', vert_align);
    end
else
    % Error message
    error('Error: Please make sure the number of labels is the same as the number of points.');
end

% Axis limits
axis square;
axis([-axes_limit, axes_limit, -axes_limit, axes_limit]);
axis off;
end

%%%%%
function [x,y,z] = pol2cartvect(th,r,z)
if size(th,1) == size(r,1) && size(th,2) == size(r,2)
    x = r.*cos(th); y = r.*sin(th);
else
    x = r*cos(th); y = r*sin(th);
end
end
%%%%%
