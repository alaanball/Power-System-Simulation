
dat1 = [1 2 0 0.1 0.002];
dat2 = [1.5 0.8];
dat3 = [];

dt = 0.001; 
X = 0.2*1j;

[S V] = NewtonRaphson(dat1, dat2, dat3);
E = V(1) + (V(2) - V(1)) * X / (dat1(4) * 1j);
%E = 1.2773;
%V(2) = 0.9083 - 0.1471 * 1j;
%d0 = 0.2323;
d0 = angle(E);
Pm = abs(E) * abs(V(2)) * sin(d0) / (0.2 + 0.1);
Pmax1 = Pm / sin(d0);
Pmax2 = Pmax1 / 6;
Pmax3 = 0.75 * Pmax1;
dmax = pi - asin(Pm / Pmax3);

data2 = [d0, 0, Pmax1, Pmax2, Pmax3, Pm, 5, 50, 0, 0, 0.5, 5, dt];
dclr = RK4(data2);

A1 = Pm * (dclr - d0) + Pmax2 * (cos(dclr) - cos(d0));
A2 = -Pmax3 * (cos(dmax) - cos(dclr)) - Pm * (dmax - dclr);