function [x_arrowed, y_arrowed] = arrowed_line(x_data, y_data, arrow_window_length, x_scale, y_scale)
% arrowed_line function returns data with arrowheads inserted at
% frequencies determined by the arrow_window_length parameter. x_Scale and
% y_scale parameters are used for manipulating the size and shape of the
% arrowheads. Unlike other solutions, this function actually inserted
% arrowheads into the data itself so the arrows remain in place despite
% window resizing.

% x_data inputted as a column vector.
% y_data can be a one dimensional column vector or contain multiple
% columns for plotting multiple arrowed lines.

if length(x_data) ~= length(y_data)
    error('X and Y data lengths unequal, check input data.');
end

[n_y_rows, n_y_cols] = size(y_data);
[n_x_rows, n_x_cols] = size(x_data);

n_ys = min([n_y_rows n_y_cols]); % assume the length of the data will be greater than the number of vars.
n_xs = min([n_x_rows n_x_cols]); % number of x-axis values.

len_data = length(y_data);

n_arrows = floor(len_data / arrow_window_length); % number of arrowheads to be inserted.

dummy_arrow = [0 0;-1 -1;0 0;-1 1;0 0]; % create dummy arrow pointing right.

if n_ys == 1

    data = [x_data y_data];

    x_range = max(x_data) - min(x_data); y_range = max(y_data) - min(y_data);

    n_arrows = floor(len_data / arrow_window_length);

    i_win_start = 1;

    for i_arrow = 1:n_arrows

        i_win_end = i_win_start + arrow_window_length;

        i_win_local = 1;

        for i_win_global = i_win_start:i_win_end
            if i_win_local == arrow_window_length/2 % if we reach midpoint of the window
                % now need to shift all the data four rows down and insert the
                % arrow head.

                pt = data(i_win_global, :);

                dx = data(i_win_end, 1) - data(i_win_start, 1);
                dy = data(i_win_end, 2) - data(i_win_start, 2); % change over the window

                theta = atan2(dy/y_range, dx/x_range); % get arrow angle

                rotm = [cos(theta) -sin(theta);sin(theta) cos(theta)]; % rotation matrix to rotate arrow

                arrow = ((dummy_arrow * rotm').*([x_range/x_scale y_range/y_scale])) + pt; % transform dummy arrow

                data = [data(1:i_win_global,:);
                    arrow; % insert arrow in data
                    data(i_win_global+1:end,:)];

                i_win_start = i_win_end + length(dummy_arrow); % advance starting index

                break % arrow head zeros inserted for this window
            else
                i_win_local = i_win_local + 1;
            end

        end

    end

    x_arrowed = data(:,1); y_arrowed = data(:,2);

elseif n_ys > 1 && n_xs == 1 % y data is multidimensional but distributed over the same x values.

    % Get the ranges of the figure plot (i.e. across all dimensions).
    y_range = range(y_data,'all'); x_range = range(x_data,'all');

    y_arrowed = zeros((len_data+(n_arrows*length(dummy_arrow))), n_ys);

    for i_var = 1:n_ys

        if n_y_cols < n_y_rows
            y_data_ith_var = y_data(:, i_var);
        elseif n_y_rows < n_y_cols
            y_data_ith_var = y_data(i_var, :)'; % want column vectors to concatonate.
        end

        data = [x_data y_data_ith_var];

        i_win_start = 1;

        for i_arrow = 1:n_arrows

            i_win_end = i_win_start + arrow_window_length;

            i_win_local = 1;

            for i_win_global = i_win_start:i_win_end
                if i_win_local == arrow_window_length/2 % if we reach midpoint of the window
                    % now need to shift all the data four rows down and insert the
                    % arrow head.

                    pt = data(i_win_global, :); % columns are the dimensions

                    dx = data(i_win_end, 1) - data(i_win_start, 1);
                    dy = data(i_win_end, 2) - data(i_win_start, 2); % change over the window

                    theta = atan2(dy/y_range, dx/x_range); % get arrow angle

                    rotm = [cos(theta) -sin(theta);sin(theta) cos(theta)]; % rotation matrix to rotate arrow

                    arrow = ((dummy_arrow * rotm').*([x_range/x_scale y_range/y_scale])) + pt; % transform dummy arrow

                    data = [data(1:i_win_global,:);
                        arrow; % insert arrow in data
                        data(i_win_global+1:end,:)];

                    i_win_start = i_win_end + length(dummy_arrow); % advance starting index

                    break % arrow head zeros inserted for this window
                else
                    i_win_local = i_win_local + 1;
                end

            end

        end

        x_arrowed = data(:,1); % not expecting differing x values
        y_arrowed(:,i_var) = data(:,2);

    end
end