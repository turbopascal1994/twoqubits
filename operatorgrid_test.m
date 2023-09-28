clear;
step = 1;
w = 1;
T1 = 8;
T2 = 8;
N1 = 4;
N2 = 4;
phi = 0;
waitq1 = 0;
waitq2 = 0;
% if bip = 0, unipolar
% if bip = 1, bipolar
bip1 = 0;
bip2 = 1;

[res,a1,a2] = ...
    OperatorGrid(w, T1, T2, N1, N2, phi, waitq1, waitq2, bip1, bip2, step);

function [arr, arr1, arr2] = OperatorGrid(w, T1, T2, N1, N2, ...
        phi, waitq1, waitq2, bip1, bip2, step)
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % combines two grids and decides which combined operator  %
    % should be applied on each grid step                     %
    % 1 grid step = 1 grs                                     %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % creating variables
    arr1 = [];
    arr2 = [];
    Nw = floor(w/step); % pulse width (in grs)
    NT1 = floor(T1/step) - Nw; % distance b/w pulses on Q1 (in grs)
    NT2 = floor(T2/step) - Nw; % distance b/w pulses on Q2 (in grs)
    sp = ones(1,Nw); % +1 pulse
    sm = ones(1,Nw)*(-1); % -1 pulse
    s0_1un = zeros(1,NT1); % distance b/w unipolar pulses on Q1
    s0_2un = zeros(1,NT2); % distance b/w unipolar pulses on Q2
    s0_1bip = zeros(1,floor((NT1 - Nw)/2)); % distance b/w bipolar pulses on Q1
    s0_2bip = zeros(1,floor((NT2 - Nw)/2)); % distance b/w bipolar pulses on Q2
    % phase between pulses
    if phi > 0
        phi_arr = zeros(1,phi);
    end
    % pulse sequence delay
    if waitq1 > 0
        wait1_arr = zeros(1,waitq1);
    end
    if waitq2 > 0
        wait2_arr = zeros(1,waitq2);
    end
    % the last pulse is not placed in the loop, because we don't need 
    % to wait additional time T1 or T2 after it
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
    % wait time after pulses
    if waitq1 > 0
        arr1 = [arr1 wait1_arr];
    end
    if waitq2 > 0
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
    arr = zeros(1,length(arr1));
    % combining into one string
    % 0 = 00, 1 = 01, 2 = 10, 3 = 11
    % 4 = -10, 5 = 0-1, 6 = -11, 7 = 1-1, 8 = -1-1
    for k = 1:1:length(arr1)
        if arr1(k) == 0 && arr2(k) == 0
            arr(k) = 0;
        elseif arr1(k) == 0 && arr2(k) == 1
            arr(k) = 1;
        elseif arr1(k) == 1 && arr2(k) == 0
            arr(k) = 2;
        elseif arr1(k) == 1 && arr2(k) == 1
            arr(k) = 3;
        elseif arr1(k) == -1 && arr2(k) == 0
            arr(k) = 4;
        elseif arr1(k) == 0 && arr2(k) == -1
            arr(k) = 5;
        elseif arr1(k) == -1 && arr2(k) == 1
            arr(k) = 6;
        elseif arr1(k) == 1 && arr2(k) == -1
            arr(k) = 7;
        elseif arr1(k) == -1 && arr2(k) == -1
            arr(k) = 8;
        end
    end
end