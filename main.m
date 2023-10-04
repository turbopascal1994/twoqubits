clear;
N = 3;
val = 2*pi*1e9;
% time grid step
tstep = 5e-14;
% main qubit frequencies
w1 = 5*val;
w2 = 5.2*val;
% anharmonicities
mu1 = 0.25*val;
mu2 = 0.4*val;
% interqubit strength
g = 0.02*val;
% qubit capacities
Cq1 = 1e-12;
Cq2 = 1e-12;
% connection capacities
Cc1 = 4e-16;
Cc2 = 4e-16;
% pulse generation frequencies
wg1 = w1;
wg2 = w2;
% pulse width
tau = 4e-12;
% number of pulses
N1 = 80;
N2 = 0;
% phase (number of grid steps paused on Q2)
phi = 0;
% wait time after pulse
waitq1 = 0;
waitq2 = 0;
% if bip = 0, unipolar
% if bip = 1, bipolar
bip1 = 0;
bip2 = 0;
init = '00';
operation = 'h0';
F = 0;
[Prob00, Prob10, Prob01, Prob20, Prob02, Prob11, F] = ...
     SimulateRegular(N, w1, w2, mu1, mu2, g, Cq1, Cq2, Cc1, Cc2, wg1, wg2,...
     tau, N1, N2, phi, waitq1, waitq2, bip1, bip2, tstep, init, operation);

%SimulateRegularDraw(N, w1, w2, mu1, mu2, g, Cq1, Cq2, Cc1, Cc2, ...
%    wg1, wg2, tau, N1, N2, phi, waitq1, waitq2, bip1, bip2, tstep, init)

% disp(['00: ', num2str(Prob00)]);
% disp(['10: ', num2str(Prob10)]);
% disp(['01: ', num2str(Prob01)]);
% disp(['20: ', num2str(Prob20)]);
% disp(['11: ', num2str(Prob11)]);
% disp(['02: ', num2str(Prob02)]);
disp(['F = ', num2str(F)]);