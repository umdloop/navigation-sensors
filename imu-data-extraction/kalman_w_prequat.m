load('ds245.mat')
load('ds200.mat')
load('ds133.mat')
load('ds114.mat')
load('ds10.mat')


% Kalman filter implementation for IMU data with isolated acceleration
% for x y and z directions with tilt correction

% matrix vector assignment:
% [ (velocity)  (acceleration)]
% state space
% [ vel(m/s) acceleration(m/s) ]

% sensor data
sensorA = sgimpA;
sensorW = sgimpW;

% noise profile for accleration
np_dataA = npA;

% -------------- BEGIN SETUP --------------
% state, data, and environment parameters
del_t = 1/114;
g = 9.8067;
st_a_var = .5;                  % acceleration state variance
st_v_var = 4;                   % velocity state variance
v_a_covar = 1;                  % v_a_covar
sensor_cov = cov(np_dataA);     % sensor covariance (acceleration)
% -------- TILT CORRECTION --------

g_curr = [0 0 0 1];
g_track = zeros(size(sensorW, 1), 4);
for i = 1:size(sensorA, 1)
    w = sensorW(i,:); % grab w vector
    [magn, g_curr] = quaternion_rot(w, g_curr, del_t); 
    g_track(i,:) = g_curr;
end

corrA = zeros(size(sensorA,1), 3);

for i = 1:size(sensorA, 1)
    corrA(i,1) = sensorA(i, 1) + g_track(i, 2);
    corrA(i,2) = sensorA(i, 2) + g_track(i, 3);
    corrA(i,3) = sensorA(i, 3) - g_track(i, 4);
end


%  -------  BEGIN COVARIANCE SETUP  -------

% state var
state = zeros(6,1);
states = zeros(size(sensorA,1), 6);
states_w_tcorr = zeros(size(sensorA,1), 6);

% prediction matrix
F_k = [eye(3), del_t*g*eye(3);
       zeros(3), eye(3)];
   
% Q_k - Prediction uncertainty
% Q_k =   0.0000001*eye(6);
Q_k =  0.0000001*[eye(3), eye(3);
                   eye(3), eye(3)];

% P_k - state covariance matrix
% P_k = [ st_v_var*eye(2), v_a_covar*eye(2);
%              v_a_covar*eye(2), st_a_var*eye(2)];
P_k = zeros(6);

% sensor covariance matrix
R_k = .1*sensor_cov(1:3, 1:3);
        
% H_k - transformation from state to measurement space
H_k = (1/g)* [0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 1, 0;
              0, 0, 0, 0, 0, 1];
% -------------- END SETUP --------------

for index = 1:size(sensorA,1)
    % ~~~ prediction step ~~~
    state = F_k * state;
    % P_k = F_k * P_k * transpose(F_k);           % without additional uncertainty
    P_k = F_k * P_k * transpose(F_k) + Q_k;   % with additional uncertainty
    
    % ~~~ update step ~~~
    z_k = transpose(corrA(index, 1:3));
    K = P_k*transpose(H_k) * pinv(H_k*P_k*transpose(H_k) + R_k);
    
    state_upd = state + K*(z_k - H_k*state);
    P_upd = P_k - K * H_k * P_k;
    
    states(index,:) = transpose(state_upd);
    state = state_upd;
    P_k = P_upd;
end

plot(sensorA(:, 1:3));
ylim([-1 1]);
title('Raw Acceleration Sensor Data');

% figure
% plot(states(:,4:6))
% ylim([-.1 .1]);
% title('Kalman Filtered Acceleration (w/o tilt correction)');

figure
plot(states(:,4:6))
ylim([-10 10]);
title('Kalman Filtered Acceleration');
ylabel('Acceleration (m/s^2)');

function [magn, g_corrected] = quaternion_rot(w, g_curr, t_step)
    magn = sqrt( w(1)^2 + w(2)^2 + w(3)^2 ); % find magnitude of w
    w_norm = w/magn; % normalize w
    theta = deg2rad(magn * t_step); % angle of rotation
    q = [cos(theta/2) sin(theta/2)*w_norm(1,:)]; % quaternion equivalent of w
    q_conj = [cos(theta/2) -sin(theta/2)*w_norm(1,:)]; % conjugate
    g_corrected = quatmultiply( quatmultiply(q, g_curr), q_conj);
end
