clc;
clear all;
close all;

% Works only for the 3 bus system done in class
% V and S are found using the NR method
%V = [1; 0.89847 - 0.10904j; 1.04959 + 0.02952j];
%S =  [1.00000 + 0.04722i;
 % -4.99999 - 3.99999i;
 %  4.00000 + 5.37058i];

BusData = input('Enter bus data: from to R X hlc: ');
PQdat = input('Enter PQ bus data: P Q for each PQ bus: ');
PVdat = input('Enter PV bus data: P V for each PV bus: ');
%machDat = input('Enter machine data: node | internal Z | H: ');
SSdat = input('Enter fs, tin, tfault, tclear, tfinal, dt: ');
%faultDat = input('Grounded nodes during fault: ');
%cbDat = input('Lines removed by CB: ');

if(min(size(BusData)) == 0 || min(size(PQdat)) == 0 || min(size(PVdat)) == 0)
{ 
  % default values
  % BUS DATA: from to R X hlc
  BusData = [1 3 0 1/40 0; 2 3 0 1/20 0; 1 2 0 1/20 0];

  % PQ bus data: P Q for each PQ bus
  PQdat = [-5 -4];

  % PV bus data: P V for each PV bus
  PVdat = [4 1.05]; 
}
endif

%Find Ybus of data
[S V Ybus] = NewtonRaphson(BusData, PQdat, PVdat);

% Internal impedances of Buses with machines
Xd1 = 1.2; Xd3 = 0.2;

% Current vector
J = Ybus * V;

E1 = V(1) + J(1) * Xd1 * 1j;
E3 = V(3) + J(3) * Xd3 * 1j;

Pm1 = abs(E1) * abs(V(1)) * sin(angle(E1) - angle(V(1))) / abs(Xd1) ;
Pm3 = abs(E3) * abs(V(3)) * sin(angle(E3) - angle(V(3))) / abs(Xd3) ;

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

% Data required for SS analysis 
%SSdat = input('Enter fs, tin, tfault, tclear, tfinal, dt');
SSdat = [50 0 0 0.2 5 0.0001];

fs = SSdat(1);
tin = SSdat(2); 
tfault = SSdat(3);
tclear = SSdat(4);
tfinal = SSdat(5);
dt = SSdat(6);

% Data Format: initial delta (in radians), initial \omega_r, H  
data1 = [angle(E1) 0 7.5];
data3 = [angle(E3) 0 5  ];

din1 = data1(1);  
wrin1 = data1(2);
H1 = data1(3); 
M1 = 2*pi*fs / H1;

din3 = data3(1); 
wrin3 = data3(2);
H3 = data3(3); 
M3 = 2*pi*fs / H3;


n = ceil((tfinal - tin)/dt);
t(1) = tin;
t = zeros(1, n);

d1 = zeros(1, n);
wr1 = zeros(1, n);
d1(1) = din1; 
wr1(1) = wrin1; 

d3 = zeros(1, n);
wr3 = zeros(1, n);
d3(1) = din3; 
wr3(1) = wrin3; 



Vk = zeros(size(Ybus1, 1), 1);
Vk(1:3) = V;
Vk(4) = E1;
Vk(5) = E3;

% Ybus after fault clearance (line connecting nodes 1 and 3 is removed)

Ypf = Ybus1;
Ypf(1, 1) = Ypf(1, 1) + Ypf(1, 3);
Ypf(3, 3) = Ypf(3, 3) + Ypf(1, 3);
Ypf(1, 3) = 0;
Ypf(3, 1) = 0;

%ind is a row of node numbers

ind = 1: size(Ybus1, 1);

for k = 1: n

  if ((t(k) < tclear) && (t(k) >= tfault))	
	  % node 1 is grounded
    Vk(1) = 0;
    Vk(2) = -Ybus1(2, ind ~= 2) * Vk(ind ~= 2) / Ybus1(2, 2); 
    Vk(3) = -Ybus1(3, ind ~= 3) * Vk(ind ~= 3) / Ybus1(3, 3);	
  endif
 
  if ((t(k) >= tclear)) 
    Vk(1) = -Ypf(1, ind ~= 1) * Vk(ind ~= 1) / Ypf(1, 1);
    Vk(2) = -Ypf(2, ind ~= 2) * Vk(ind ~= 2) / Ypf(2, 2); 
    Vk(3) = -Ypf(3, ind ~= 3) * Vk(ind ~= 3) / Ypf(3, 3);	
  endif
  
  Pmax1 = abs(Vk(4) * Vk(1) / Xd1);
  Pmax3 = abs(Vk(5) * Vk(3) / Xd3);
  
  % For machine 1
  k11 = (wr1(k)) * dt;
  l11 = M1 * (Pm1 - Pmax1 * sin(d1(k))) * dt;
   
  k21 = (wr1(k) + l11/2) * dt;
  l21 = M1 * (Pm1 - Pmax1 * sin(d1(k) + k11/2)) * dt;
   
  k31 = (wr1(k) + l21/2) * dt;
  l31 = M1 * (Pm1 - Pmax1 * sin(d1(k) + k21/2)) * dt;
   
  k41 = (wr1(k) + l31) * dt;
  l41 = M1 * (Pm1 - Pmax1 * sin(d1(k) + k31)) * dt;

  wr1(k+1) = wr1(k) + (l11 + 2*(l21 + l31) + l41)/6;
  d1(k+1) = d1(k) + (k11 + 2*(k21 + k31) + k41)/6;
  
  % For machine 3
  k13 = (wr3(k)) * dt;
  l13 = M3 * (Pm3 - Pmax3 * sin(d3(k))) * dt;
   
  k23 = (wr3(k) + l13/2) * dt;
  l23 = M3 * (Pm3 - Pmax3 * sin(d3(k) + k13/2)) * dt;
   
  k33 = (wr3(k) + l23/2) * dt;
  l33 = M3 * (Pm3 - Pmax3 * sin(d3(k) + k23/2)) * dt;
   
  k43 = (wr3(k) + l33) * dt;
  l43 = M3 * (Pm3 - Pmax3 * sin(d3(k) + k33)) * dt;

  wr3(k+1) = wr3(k) + (l13 + 2*(l23 + l33) + l43)/6;
  d3(k+1) = d3(k) + (k13 + 2*(k23 + k33) + k43)/6;
  
  Vk(4) = abs(E1) * exp(1j * d1(k));
  Vk(5) = abs(E3) * exp(1j * d3(k));
  
  t(k + 1) = t(k) + dt;


endfor
   
figure(1);
plot(t, d1);
figure(2);
plot(t, wr1);

figure(3);
plot(t, d3);
figure(4);
plot(t, wr3);
   
   
