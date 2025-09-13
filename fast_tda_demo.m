% fast_tda demo
% Authors: Charlton Rolle and Jason D. Bakos
% takes a dataset, in the form crafted by Mohsen Gol Zardian
% finds the time delay and window using the method prescribed by Arman Razmarashooli
% fits an ellipse over each window using the conic Fitzgibbon method
% converts each set of conic parameters to "parametric parameters", i.e.
% converts quadratic coefficients to ellipse center, size, and angle
% plots all data and saves to a video file

% read data
data = readtable('Temp_21_output_time.csv');
time = data.Time_s;
output = data.Output;
fprintf('Data loaded: %d rows\n', length(time));

% find the sample rate
dt = time(2)-time(1);

% plot the spectrum
Fs = 1 / dt; % Sampling frequency (Hz)
N = length(output); % Number of samples
Y = fft(output); % Compute FFT
Y_magnitude = abs(Y/N); % Normalize magnitude
if mod(N, 2) == 0
    % Even number of samples
    f = (0:N/2)*(Fs/N); % Frequency vector for single-sided spectrum
    Y_magnitude = Y_magnitude(1:N/2+1); % Take single-sided spectrum
    Y_magnitude(2:end-1) = 2*Y_magnitude(2:end-1); % Double amplitudes (except DC and Nyquist)
else
    % Odd number of samples
    f = (0:(N-1)/2)*(Fs/N); % Frequency vector for single-sided spectrum
    Y_magnitude = Y_magnitude(1:(N+1)/2); % Take single-sided spectrum
    Y_magnitude(2:end) = 2*Y_magnitude(2:end); % Double amplitudes (except DC)
end
figure;
plot(f, Y_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Single-Sided Amplitude Spectrum of Output Signal');
grid on;

% Extract peak frequencies
% assume Arman is defining Fmax as the absolute peak and Fmin as the lowest
% peak having threshold 10% of maximum magnitude
max_mag = max(Y_magnitude);
[peaks, peak_props] = findpeaks(Y_magnitude, ...
    'MinPeakHeight', 0.05 * max_mag, ...
    'MinPeakProminence', 0.1 * max_mag, ...
    'MinPeakDistance', 5);
peak_freqs = f(peak_props);
Fmin = min(peak_freqs);
Fmax = max(peak_freqs);
fprintf ("F_max = %0.2f Hz, F_min = %0.2f Hz\n",Fmax,Fmin);
hold off;

% set up Fast TDA parameters
time_delay = .25/Fmax;
time_delay_in_samples = round(time_delay / dt); % from Arman
window_duration = 1/Fmin; % from Arman
num_points_per_window = round(window_duration / dt); 
num_windows = length(time) - num_points_per_window;
step_size = 50;

% delete the existing video file
if exist('myVideo.mp4', 'file')
    delete('myVideo.mp4');
end

% get ranges for acceleration plot
min_accel = min(output, [], 'all');
max_accel = max(output, [], 'all');

% Create an invisible figure with larger size for better visibility in video
fig = figure('Visible', 'off', 'Renderer', 'opengl', 'Position', [100 100 1200 900]);

% Create subplots (axes) once outside the loop
ax1 = subplot(4,1,1);
plot(ax1, time, output, 'r-', 'LineWidth', 2);
hold(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'time');
ylabel(ax1, 'acceleration');
xlim(ax1, [min(time) max(time)]);
ax2 = subplot(4,1,2);
ax3 = subplot(4,1,3);
ax4 = subplot(4,1,4);

% set up the conic arguments plot
line_handles_conic = gobjects(5,1);
for j = 1:5
    line_handles_conic(j) = plot(ax3, NaN, NaN, 'LineWidth', 1);
    hold(ax3, 'on');
end
% note that we don't plot f, since it's dynamic range is too high
legend(ax3, {"a", "b", "c", "d", "e"});
xlabel(ax3, "time");
xlim(ax3, [min(time) max(time)]);

% set up the parametric parameters plot
line_handles_parametric = gobjects(5,1);
for j = 1:5
    line_handles_parametric(j) = plot(ax4, NaN, NaN, 'LineWidth', 1);
    hold(ax4, 'on');
end
legend(ax4, {"center\_x", "center\_y", "semi-major", "semi-minor", "angle"});
xlabel(ax4, "time");
xlim(ax4, [min(time) max(time)]);

% set up video with lower frame rate to slow down playback
video = VideoWriter('myVideo.mp4', 'MPEG-4');
video.FrameRate = 10;
open(video);

% allocate space for ellipse conic parameters
ellipse_params = zeros(6, num_windows);

% allocate space for parametric parameters
ellipse_params_parametric = zeros(5, num_windows);

% used for the sliding window display
rect_handle = [];

for i = 1:step_size:num_windows
    % build pointcloud
    P = [output(i:i+num_points_per_window-1-time_delay_in_samples), ...
         output(i+time_delay_in_samples:i+num_points_per_window-1)];
    
    % fit the ellipse
    ellipse_params(:,i) = fit_ellipse(P);
    
    % move the sliding window
    if ~isempty(rect_handle)
        delete(rect_handle);
    end
    x_left = i * dt;
    y_bottom = min_accel;
    width = window_duration;
    height = max_accel - min_accel;
    x_rect = [x_left x_left+width x_left+width x_left];
    y_rect = [y_bottom y_bottom y_bottom+height y_bottom+height];
    rect_handle = patch(ax1, x_rect, y_rect, [0.5 0.7 1], 'FaceAlpha', 0.3, ...
                        'EdgeColor', 'b', 'LineWidth', 1);

    % plot the pointcloud and fitted ellipse
    cla(ax2);
    plot_ellipse(P, ellipse_params(:,i), ax2);

    % plot the conic parameters
    t_values = (1:step_size:i) * dt;
    for j = 1:5 % omit 'f'
        set(line_handles_conic(j), 'XData', t_values, 'YData', ellipse_params(j, 1:step_size:i));
    end
    xlim(ax3, [min(time) max(time)]);

    % plot the parametric parameters
    [ellipse_params_parametric(1,i),...
        ellipse_params_parametric(2,i),...
        ellipse_params_parametric(3,i),...
        ellipse_params_parametric(4,i),...
        ellipse_params_parametric(5,i)] = conic_to_parametric(ellipse_params(:,i));

    % check if we fail to find ellipse parametric parameters
    if any(isnan(ellipse_params_parametric(:,i)))
        ellipse_params_parametric(:,i) = zeros(5,1);
        x_left = i * dt;
        y_bottom = min(ellipse_params_parametric,[],"all");
        width = window_duration;
        height = max(ellipse_params_parametric,[],"all") - min(ellipse_params_parametric,[],"all");
        x_rect = [x_left x_left+width x_left+width x_left];
        y_rect = [y_bottom y_bottom y_bottom+height y_bottom+height];
        patch(ax4, x_rect, y_rect, [0.5 0.7 1], ...
                        'EdgeColor', 'r', 'LineWidth', 1);
    end

    for j = 1:4
        set(line_handles_parametric(j), 'XData', t_values, 'YData', ellipse_params_parametric(j, 1:step_size:i));
    end
    xlim(ax4, [min(time) max(time)]);

    drawnow;
    frame = getframe(fig);
    writeVideo(video, frame);
end

close(video);
%close(fig);

function ellipse_params = fit_ellipse(P)
    D = [P(:,1).^2, P(:,1).*P(:,2), P(:,2).^2, P(:,1), P(:,2), ones(size(P,1),1)];
    S = D' * D;
    C = zeros(6,6);
    C(1,3) = 2; C(2,2) = -1; C(3,1) = 2;
    [eigvecs, eigvals] = eig(S, C);
    eigvals = diag(eigvals);

    finite_idx = isfinite(eigvals);
    pos_idx = eigvals > 0;
    idx = find(pos_idx & finite_idx);
    
    if length(idx) ~= 1
        warning('No unique positive finite eigenvalue; using first valid.');
        idx = find(finite_idx, 1);
        if isempty(idx)
            ellipse_params = zeros(6,1);
            return;
        end
    end

    if length(idx) > 1
        1;
    end
    
    v = eigvecs(:, idx(1));
    % Enforce constraint a' C a = 1
    mu = 1 / sqrt(v' * C * v);
    ellipse_params = mu * v;
    
    % Validate ellipse: check 4ac - b^2 ≈ 1 and discriminant
    a = ellipse_params(1); b = ellipse_params(2); c = ellipse_params(3);
    if abs(4*a*c - b^2 - 1) > 1e-5 || (b^2 - 4*a*c) >= 0
        warning('Fit may not be a valid ellipse; forcing to zero.');
        ellipse_params = zeros(6,1);
    end
end

function [] = plot_ellipse(P, ellipse_params, ax)
    axes(ax);
    mins = min(P); maxs = max(P);
    xlim([mins(1) maxs(1)]); ylim([mins(2) maxs(2)]);
    scatter(P(:,1), P(:,2), '.');
    hold on;
        axes(ax);
    mins = min(P); maxs = max(P);
    xlim([mins(1) maxs(1)]); ylim([mins(2) maxs(2)]);

    if ~any(isinf(ellipse_params))

        % Evaluate the ellipse equation on the grid
        a = ellipse_params(1);
        b = ellipse_params(2);
        c = ellipse_params(3);
        d = ellipse_params(4);
        e = ellipse_params(5);
        f = ellipse_params(6);
       
        ellipse = @(x, y) a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f;
        fimplicit(ellipse, [mins(1) maxs(1) mins(2) maxs(2)], 'LineWidth', 2);

        % Plot the fitted ellipse
        %contour(X, Y, Z, [0, 0], 'LineWidth', 2);
        %axis equal;
        %title('Fitted Ellipse and Original Points');
        %xlabel('x');
        %ylabel('y');
    end

    hold off;
end

function [center_x, center_y, semi_major, semi_minor, angle] = conic_to_parametric(params)
    a = params(1); b = params(2); c = params(3);
    d = params(4); e = params(5); f = params(6);
    
    delta = b^2 - 4*a*c;
    if delta >= 0 || abs(a) < 1e-10 || abs(c) < 1e-10
        center_x = NaN; center_y = NaN; semi_major = NaN; semi_minor = NaN; angle = NaN;
        return;
    end
    
    denom = b^2 - 4*a*c;
    center_x = (2*c*d - b*e) / denom;
    center_y = (2*a*e - b*d) / denom;
    center = [center_x; center_y];
    
    if abs(b) < 1e-10
        angle = 0;
    elseif abs(a - c) < 1e-10
        angle = pi/4;
    else
        angle = 0.5 * atan(b / (a - c));
    end
    
    a_prime = a*center_x^2 + b*center_x*center_y + c*center_y^2 + d*center_x + e*center_y + f;
    % if a_prime >= 0 || abs(a_prime) < 1e-10
    %     center = NaN; semi_major = NaN; semi_minor = NaN; angle = NaN;
    %     return;
    % end
    
    lambda1 = (a + c + sqrt((a - c)^2 + b^2)) / 2;
    lambda2 = (a + c - sqrt((a - c)^2 + b^2)) / 2;
    semi_major = sqrt(-a_prime / lambda1);
    semi_minor = sqrt(-a_prime / lambda2);
    
    if semi_minor > semi_major || semi_major < 1e-5 || semi_minor < 1e-5
        temp = semi_major;
        semi_major = semi_minor;
        semi_minor = temp;
        angle = angle + pi/2;
    end
    if semi_major > 1e5 || semi_minor > 1e5
        semi_major = NaN; semi_minor = NaN;
    end
end
