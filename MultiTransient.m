clc;
clear all;
close all;

% V and S are obtained using the Newton-Raphson method for the bus given below
V = [1; 0.89847 - 0.10904j; 1.04959 + 0.02952j];
S =  [1.00000 + 0.04722i;
  -4.99999 - 3.99999i;
   4.00000 + 5.37058i];
   
% BUS DATA: from to R X hlc
BusData = [1 3 0 1/40 0; 2 3 0 1/20 0; 1 2 0 1/20 0];

%Find Ybus of data
Ybus = CalcYbus(BusData);

% Internal impedances of Buses with machines
Xd1 = 1.2; Xd3 = 0.2;

% Current vector
J = Ybus * V;

% Find internal EMFs of machines (at buses 1 and 3)
E1 = V(1) + J(1) * Xd1 * 1j;
E3 = V(3) + J(3) * Xd3 * 1j;

% Find steady state mechanical power at each machine
Pm1 = abs(E1 * V(1) / Xd1) * sin(angle(E1) - angle(V(1)));
Pm2 = abs(E3 * V(3) / Xd2) * sin(angle(E3) - angle(V(3)));

Y2 = conj(S(2)) / abs(V(2))^2;

Ybus1 = Ybus;
Ybus1(5, 5) = -1j/Xd3; 
Ybus1(4, 4) = -1j/Xd1; 
Ybus1(1, 1) = Ybus1(1, 1) + -1j/Xd1;
Ybus1(3, 3) = Ybus1(3, 3) + -1j/Xd3;
Ybus1(1, 4) = 1j/Xd1;
Ybus1(4, 1) = 1j/ Xd1;
Ybus1(3, 5) = 1j/ Xd3;
Ybus1(5, 3) = 1j/ Xd3;
Ybus1(2, 2) = Ybus1(2, 2) + Y2;
  
data2 = [angle(E1) 0 Pm1  50 0 0 0.2 5 0.0001];
data3 = [angle(E3) 0 Pm3  50 0 0 0.2 5 0.0001];

%inp = input('Enter d1, wr1, Pmax1, Pmax2, Pmax3, Pm, H, fs, tin, tfault, tclear, tfinal, dt: ');
inp = data2;

d1 = inp(1) * pi/180; 
wr1 = inp(2);
Pmax1 = inp(3); 
Pmax2 = inp(4); 
Pmax3 = inp(5);
Pm = inp(6); 
H = inp(7); 
fs = inp(8);
tin = inp(9); 
tfault = inp(10);
tclear = inp(11);
tfinal = inp(12);
dt = inp(13);

M = 2*pi*fs / H;
n = ceil((tfinal - tin)/dt);



d = zeros(1, n);
wr = zeros(1, n);
t = zeros(1, n);
d(1) = d1; 
wr(1) = wr1; 
t(1) = tin;

for k = 1: n
  
  if (t(k) < tfault)
    Pmax = Pmax1;
  endif
  
  if ((t(k) < tclear) && (t(k) >= tfault))
    Pmax = Pmax2;
  endif
 
 if ((t(k) >= tclear))
    Pmax = Pmax3;
  endif
  
  %if(t(k) <= tclear + dt && t(k) >= tclear - dt)
  %  dclear = d(k);
  %endif

  if(abs(t(k) - tclear) <= 0.001)
    dclear = d(k);
  endif
  
  k1 = (wr(k)) * dt;
  l1 = M * (Pm - Pmax * sin(d(k))) * dt;
   
  k2 = (wr(k) + l1/2) * dt;
  l2 = M * (Pm - Pmax * sin(d(k) + k1/2)) * dt;
   
  k3 = (wr(k) + l2/2) * dt;
  l3 = M * (Pm - Pmax * sin(d(k) + k2/2)) * dt;
   
  k4 = (wr(k) + l3) * dt;
  l4 = M * (Pm - Pmax * sin(d(k) + k3)) * dt;

  wr(k+1) = wr(k) + (l1 + 2*(l2 + l3) + l4)/6;
  d(k+1) = d(k) + (k1 + 2*(k2 + k3) + k4)/6;
  t(k + 1) = t(k) + dt;


endfor
   
figure(1);
plot(t, d);
figure(2);
plot(t, wr);
   
   
