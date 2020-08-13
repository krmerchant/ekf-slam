function [mu, sigma] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model,
% mu: 2N+3 x 1 vector representing the state mean
% sigma: 2N+3 x 2N+3 covariance matrix
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% TODO: Compute the new mu based on the noise-free (odometry-based) motion model
% Remember to normalize theta after the update (hint: use the function normalize_angle available in tools)
N = (size(mu)(1) - 3)/2; %number of landmarsk

F = [eye(3) zeros(3,2*N)];%3x2N+3

%state transition matrix
A = [u.t*cos(mu(3)  + u.r1);
     u.t*sin(mu(3) +  u.r1); 
      u.r1 + u.r2]; %3x1
%update mu
mu = mu + F'*A; %2N+3x1 + (2N+3)x3*3x1

mu(3) = normalize_angle(mu(3)); %TODO: Go back to this
% TODO: Compute the 3x3 Jacobian Gx of the motion model
Gx_dtheta = [-u.t*sin(mu(3)  + u.r1);
            u.t*cos(mu(3) +  u.r1); 
            0];
G_x = [zeros(3,2) Gx_dtheta];
% TODO: Construct the full Jacobian G
G = eye(2*N+3) + F'*G_x*F; %DxD + Dx3*3x3*3X

% Motion noise
motionNoise = 0.1;
R = [motionNoise, 0, 0; 
     0, motionNoise, 0; 
     0, 0, motionNoise/10];

% TODO: Compute the predicted sigma after incorporating the motion
sigma = G*sigma*G'+ F'*R*F; %DxD*DXD*DxD + Dx3*3x3+3xD;       

end
