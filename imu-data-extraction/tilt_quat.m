load('ds114.mat')

% using deg90 sample - 114hz
% deg90A - acceleration
% deg90W - angular velocity

% setup
g_0 = [0 0 -1 0];
g_curr = g_0;

% -------- tweakable parameters ---------
sample = npW; % data file in use
t_step = 1/114;
deriv_order = 5;
noise_delineator = 100000; % .65 for 3rd, 1.12 for 4th, 2.2 for 5th
% ---------------------------------------

% keeps 3rd order difference for reference
% magn_track = zeros(4, 1);
magn_arr = zeros(size(sample, 1), 1);

% matrix containing all g vectors throughout computation
g_track = zeros(size(sample, 1), 4);
g_track(1,:) = g_0;
g_track(2,:) = g_0;
g_track(3,:) = g_0;
g_track(4,:) = g_0;
g_track(5,:) = g_0;
g_track(6,:) = g_0;

%diff_old = diff(magn_arr(1:4), 3)
for i = 6:size(sample, 1)
    % normalize pure q vector -> wx, wy, wz
    % conjugate is negation of normalized q vector
    
    % find magnitude of pure q vector -> rot angle
    
    % (q)(g_curr)(q*)
    w = sample(i,:); % grab w vector
    [magn, g_curr] = quaternion_rot(w, g_curr, t_step);

    % compute the snap for comparison to current jerk
    snap = diff(magn_arr(i - deriv_order:i), deriv_order)
    
    if (abs(snap) > noise_delineator)
        g_track(i,:) = g_track(i - 1, :);
        g_curr = g_track(i,:);
    else
        g_track(i,:) = g_curr;
    end    
    
    magn_arr(i) = magn;
end

% % want to limit jerk to filter noise
% diff_magn = diff(magn_track, 3);
fprintf("Final gravity vector: ");
disp(g_track(size(g_track, 1), 2:4)) % last gravity vector

plot(g_track(:, 2:end))

function [magn, g_corrected] = quaternion_rot(w, g_curr, t_step)
    magn = sqrt( w(1)^2 + w(2)^2 + w(3)^2 ); % find magnitude of w
    w_norm = w/magn; % normalize w
    theta = deg2rad(magn * t_step); % angle of rotation
    q = [cos(theta/2) sin(theta/2)*w_norm(1,:)]; % quaternion equivalent of w
    q_conj = [cos(theta/2) -sin(theta/2)*w_norm(1,:)]; % conjugate
    g_corrected = quatmultiply( quatmultiply(q, g_curr), q_conj);
end