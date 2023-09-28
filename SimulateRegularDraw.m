function SimulateRegularDraw(N, w1, w2, mu1, mu2, g, Cq1, Cq2, Cc1, Cc2, ...
    wg1, wg2, tau, N1, N2, phi, waitq1, waitq2, bip1, bip2, tstep, init)
% ignore unused vars
%#ok<*NASGU>
% ignore preallocation
%#ok<*AGROW>

% some magic
if bip1 == 1
    N1 = N1/2;
end
if bip2 == 1
    N2 = N2/2;
end

h = 1.054e-34; % planck constant
F0 = 2.06e-15; % magnetic flux quantum

% second order quantizaion
a1 = zeros(N); % a dagger
a2 = zeros(N); % a
for i = 2:1:N
    a1(i, i-1) = sqrt(i-1);
    a2(i-1, i) = sqrt(i-1);
end
I = eye(N);
aa = a1*a2;

% field operator
V0 = F0/tau; % Voltage
Amp1 = Cc1*V0*sqrt(h*w1/(2*Cq1));
Amp2 = Cc2*V0*sqrt(h*w2/(2*Cq2));
V1 = 1i*Amp1*(a2 - a1);
V2 = 1i*Amp2*(a2 - a1);

% now to the hamiltonians
HQ1 = h*w1*aa - h*mu1/2*aa*(aa - I);
HQ2 = h*w2*aa - h*mu2/2*aa*(aa - I);
Hint = h*g*kron((a1 + a2),(a1 + a2));

% No field hamiltonian
H00 = kron(HQ1,I) + kron(I,HQ2) + Hint;

% eigens
[Eigvec, Eigval] = eig(H00); 
[Eigvals, ind] = sort(diag(Eigval)); 
Eigvecs = Eigvec(:,ind); 
WF00 = Eigvecs(:,1); % |00> wavefunction
WF10 = Eigvecs(:,2); % |10> wavefunction
WF01 = Eigvecs(:,3); % |01> wavefunction
switch N
    case 2
        WF11 = Eigvecs(:,4); % |11> wavefunction
    case 3
        WF20 = Eigvecs(:,4); % |20> wavefunction
        WF02 = Eigvecs(:,5); % |02> wavefunction
        WF11 = Eigvecs(:,6); % |11> wavefunction
end

% E00 = Eigvals(1); % |00> energy
% E10 = Eigvals(2); % |10> energy
% E01 = Eigvals(3); % |01> energy
% switch N
%    case 2
%        E11 = Eigvals(4); % |11> energy
%    case 3
%        E20 = Eigvals(4); % |20> energy
%        E02 = Eigvals(6); % |02> energy
%        E11 = Eigvals(5); % |11> energy
% end

% Hamiltonian with field
H10 = H00 + kron(V1,I);
H01 = H00 + kron(I,V2);
H11 = H00 + kron(V1,I) + kron(I,V2);
Hm10 = H00 + kron(-V1,I);
H0m1 = H00 + kron(I,-V2);
Hm11 = H00 + kron(-V1,I) + kron(I,V2); 
H1m1 = H00 + kron(V1,I) + kron(I,-V2);
Hm1m1 = H00 + kron(-V1,I) + kron(I,-V2);

% Operators
U00 = UMatrix(H00, tstep, N^2);
U10 = UMatrix(H10, tstep, N^2);
U01 = UMatrix(H01, tstep, N^2);
U11 = UMatrix(H11, tstep, N^2);
Um10 = UMatrix(Hm10, tstep, N^2);
U0m1 = UMatrix(H0m1, tstep, N^2);
Um11 = UMatrix(Hm11, tstep, N^2);
U1m1 = UMatrix(H1m1, tstep, N^2);
Um1m1 = UMatrix(Hm1m1, tstep, N^2);

PulseString = OperatorGrid(tau, 2*pi/wg1, 2*pi/wg2, N1, N2, phi, ...
    waitq1, waitq2, bip1, bip2, tstep);
%scatter(1:1:length(PulseString), PulseString);
U = eye(N^2);
switch init
    case '00'
        WF = WF00;
    case '01'
        WF = WF01;
    case '10'
        WF = WF10;
    case '11'
        WF = WF11;
    case '20'
        WF = WF20;
    case '02'
        WF = WF02;
end

Probs00 = zeros(1,length(PulseString));
Probs10 = zeros(1,length(PulseString));
Probs01 = zeros(1,length(PulseString));
Probs20 = zeros(1,length(PulseString));
Probs11 = zeros(1,length(PulseString));
Probs02 = zeros(1,length(PulseString));
PulsesQ1 = zeros(1,length(PulseString));
PulsesQ2 = zeros(1,length(PulseString));
for i = 1:1:length(PulseString)
    switch PulseString(i)
        case 0
            Upulse = U00;
            PulsesQ1(i) = 0;
            PulsesQ2(i) = 0;
        case 1
            Upulse = U01;
            PulsesQ1(i) = 0;
            PulsesQ2(i) = 1;
        case 2
            Upulse = U10;
            PulsesQ1(i) = 1;
            PulsesQ2(i) = 0;
        case 3
            Upulse = U11;
            PulsesQ1(i) = 1;
            PulsesQ2(i) = 1;
        case 4 
            Upulse = Um10;
            PulsesQ1(i) = -1;
            PulsesQ2(i) = 0;
        case 5
            Upulse = U0m1;
            PulsesQ1(i) = 0;
            PulsesQ1(i) = -1;
        case 6 
            Upulse = Um11;
            PulsesQ1(i) = -1;
            PulsesQ1(i) = 1;
        case 7
            Upulse = U1m1;
            PulsesQ1(i) = 1;
            PulsesQ1(i) = -1;
        case 8
            Upulse = Um1m1;
            PulsesQ1(i) = -1;
            PulsesQ1(i) = -1;
    end
   WF = Upulse*WF;
   Probs00(i) = abs(ctranspose(WF00)*WF)^2;
   Probs10(i) = abs(ctranspose(WF10)*WF)^2;
   Probs01(i) = abs(ctranspose(WF01)*WF)^2;
   Probs20(i) = abs(ctranspose(WF20)*WF)^2;
   Probs11(i) = abs(ctranspose(WF11)*WF)^2;
   Probs02(i) = abs(ctranspose(WF02)*WF)^2;
end
%figure(1)
%plot(0/1e-9:tstep/1e-9:(length(PulsesQ1)-1)/1e-9*tstep,PulsesQ1)
%title('Pulses on Q1')

%figure(2)
%plot(0:tstep/1e-9:length(length(PulsesQ2)-1)/1e-9*tstep,PulsesQ2)
%title('Pulses on Q2')

%figure(3)
LW = 1.8;
tt = 0:tstep:(length(PulseString) - 1)*tstep;
figure(1)
hold on
plot(tt,Probs00,'LineWidth',LW);
plot(tt,Probs10,'LineWidth',LW);
hold off
legend('00','10')
figure(2)
set(gca, 'YScale', 'log')
hold on
plot(tt,Probs01,'g','LineWidth',LW);
plot(tt,Probs20,'c','LineWidth',LW);
plot(tt,Probs02,'m','LineWidth',LW);
plot(tt,Probs11,'k','LineWidth',LW);
hold off
legend('01','20','02','11','Location','southeast')

disp(['00: ', num2str(Probs00(end))]);
disp(['10: ', num2str(Probs10(end))]);
disp(['01: ', num2str(Probs01(end))]);
disp(['20: ', num2str(Probs20(end))]);
disp(['02: ', num2str(Probs02(end))]);
disp(['11: ', num2str(Probs11(end))]);

function arr = OperatorGrid(w, T1, T2, N1, N2, ...
        phi, waitq1, waitq2, bip1, bip2, step)
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % combines two grids and decides which combined operator  %
    % should be applied on each grid step                     %
    % 1 grid step = 1 grs                                     %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % creating variables
    arr1 = [];
    arr2 = [];
    Nw = floor(w/step);    % pulse width (in grs)
    NT1 = floor(T1/step) - Nw; % distance b/w pulses on Q1 (in grs)
    NT2 = floor(T2/step) - Nw; % distance b/w pulses on Q2 (in grs)
    sp = ones(1,Nw);       % array for +1 pulse
    sm = ones(1,Nw)*(-1);  % array for -1 pulse
    s0_1un = zeros(1,NT1); % array for distance b/w uni pulses on Q1
    s0_2un = zeros(1,NT2); % array for distance b/w uni pulses on Q2
    s0_1bip = zeros(1,floor((NT1 - Nw)/2)); % array for distance b/w bip pulses on Q1
    s0_2bip = zeros(1,floor((NT2 - Nw)/2)); % array for distance b/w bip pulses on Q2
    % phase between pulses (in grs)
    if phi > 0
        phi_arr = zeros(1,phi);
    end
    % placing pulses and the distance b/w them
    % note: the last pulse is placed outside the loop, because we don't 
    % need to wait additional time T1 or T2 after it
    for j = 1:1:(N1-1)
        if bip1 == 0
            arr1 = [arr1 sp]; 
            arr1 = [arr1 s0_1un];
        elseif bip1 == 1
            arr1 = [arr1 sp];
            arr1 = [arr1 s0_1bip];
            arr1 = [arr1 sm];
            arr1 = [arr1 s0_1bip];
        end
    end
    for j = 1:1:(N2-1)
        if phi > 0
            arr2 = [arr2 phi_arr];
        end
        if bip2 == 0
            arr2 = [arr2 sp]; 
            arr2 = [arr2 s0_2un];
        elseif bip2 == 1
            arr2 = [arr2 sp];
            arr2 = [arr2 s0_2bip];
            arr2 = [arr2 sm];
            arr2 = [arr2 s0_2bip];
        end
    end
    % adding the last pulse
    if N1 > 0
        arr1 = [arr1 sp];
    end
    if N2 > 0
        arr2 = [arr2 sp];
    end
    % wait time after pulses (in grs)
    if waitq1 > 0
        wait1_arr = zeros(1,waitq1);
        arr1 = [arr1 wait1_arr];
    end
    if waitq2 > 0
        wait2_arr = zeros(1,waitq2);
        arr2 = [arr1 wait2_arr];
    end
    % equalizing pulse strings by adding zeros to the lesser one
    if length(arr1) ~= length(arr2)
        difflen = abs(length(arr1) - length(arr2));
        diff = zeros(1,difflen);
        if length(arr1) > length(arr2)
            arr2 = [arr2, diff];
        else
            arr1 = [arr1, diff];
        end
    end
    arr = [];
    % combining into one string
    % legend:
    % 0 = 00, 1 = 01, 2 = 10, 3 = 11
    % 4 = -10, 5 = 0-1, 6 = -11, 7 = 1-1, 8 = -1-1
    for k = 1:1:length(arr1)
        if arr1(k) == 0 && arr2(k) == 0
            arr = [arr, 0];
        elseif arr1(k) == 0 && arr2(k) == 1
            arr = [arr, 1];
        elseif arr1(k) == 1 && arr2(k) == 0
            arr = [arr, 2];
        elseif arr1(k) == 1 && arr2(k) == 1
            arr = [arr, 3];
        elseif arr1(k) == -1 && arr2(k) == 0
            arr = [arr, 4];
        elseif arr1(k) == 0 && arr2(k) == -1
            arr = [arr, 5];
        elseif arr1(k) == -1 && arr2(k) == 1
            arr = [arr, 6];
        elseif arr1(k) == 1 && arr2(k) == -1
            arr = [arr, 7];
        elseif arr1(k) == -1 && arr2(k) == -1
            arr = [arr, 8];
        end
    end
end
end