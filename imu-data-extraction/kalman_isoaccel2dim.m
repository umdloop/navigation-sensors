load('ds245.mat')
load('ds200.mat')
load('ds133.mat')
load('ds114.mat')
load('ds10.mat')

% Kalman filter implementation for IMU data with isolated acceleration
% for just x and y direction

% matrix vector assignment:
% [ (velocity)  (acceleration in g)]
% standard deviation for velocity: +- 2 m/s

% sensor data
sensor_datA = sgimpA;

% noise profile for accleration
np_dataA = npA;

% -------------- BEGIN SETUP --------------
% state, data, and environment parameters
del_t = 1/133.33;
g = 9.8067;
st_a_var = .5;                  % acceleration state variance
st_v_var = 4;                   % velocity state variance
v_a_covar = 1;                  % v_a_covar
sensor_cov = cov(np_dataA);     % sensor covariance (acceleration)
%  -------  BEGIN COVARIANCE SETUP  -------

% state var
state = zeros(4,1);
states = zeros(size(sensor_datA,1), 4)

% prediction matrix
F_k = [eye(2), del_t*g*eye(2);
       zeros(2), eye(2)];
   
% Q_k - Prediction uncertainty
%Q_k = .00000001*eye(4);
Q_k =  0.00000001*[eye(2), eye(2);
                   eye(2), eye(2)];

% P_k - state covariance matrix
% P_k = [ st_v_var*eye(2), v_a_covar*eye(2);
%              v_a_covar*eye(2), st_a_var*eye(2)];
P_k = eye(4);

% sensor covariance matrix
R_k = sensor_cov(1:2, 1:2);
        
% H_k - transformation from state to measurement space
H_k = (1/g)* [0, 0, 1, 0;
              0, 0, 0, 1];
% -------------- END SETUP --------------

for index = 1:size(sensor_datA,1)
    % ~~~ prediction step ~~~
    state = F_k * state;
    % P_k = F_k * P_k * transpose(F_k);           % without additional uncertainty
    P_k = F_k * P_k * transpose(F_k) + Q_k;   % with additional uncertainty
    
    % ~~~ update step ~~~
    z_k = transpose(sensor_datA(index, 1:2));
    K = P_k*transpose(H_k) * pinv(H_k*P_k*transpose(H_k) + R_k);
    
    state_upd = state + K*(z_k - H_k*state);
    P_upd = P_k - K * H_k * P_k;
    
    states(index,:) = transpose(state_upd);
    state = state_upd;
    P_k = P_upd;
end

plot(linA(:, 1:2));
ylim([-.1 .1]);
title('Raw Acceleration Sensor Data');

figure
plot(states(:,3:4))
ylim([-.1 .1]);
title('Kalman Filtered Acceleration');
