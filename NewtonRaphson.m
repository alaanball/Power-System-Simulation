function [S V Ybus] = NewtonRaphson(dat1, dat2, dat3)
  % Order of buses
  % slack; PQ; PV
  % dat 1: network data - from to R X hlc
  % dat 2: P and Q for each PQ bus
  % dat 3: P and V for each PV bus
  
  data = dat1;
  Ybus = CalcYbus(data);
  nbus = size(Ybus, 1);
  maxiter = 10;
  n = 0;

  PQdat = dat2; 
  PVdat = dat3; 

  if(min(size(PVdat)) != 0)
    Psp = [PQdat(:, 1); PVdat(:, 1)];
    Qsp = [PQdat(:, 2)];
  else
    Psp = [PQdat(:, 1)] ;
    Qsp = [PQdat(:, 2)];
  endif

  nPQ = size(PQdat, 1);
  nPV = size(PVdat, 1);
  SlackVoltage = 1;
  
  if(nPV == 0)
    V = [SlackVoltage; ones(nbus - 1, 1)];
  else
    V = [SlackVoltage; ones(nPQ, 1); PVdat(:, 2)];
  endif
    
  %J1 = zeros(nPQ + nPV, nPQ + nPV);
  %J2 = zeros(nPQ + nPV, nPQ);
  %J3 = zeros(nPQ, nPQ + nPV);
  %J4 = zeros(nPQ, nPQ);

   J1 = zeros(nbus - 1, nbus - 1);
   J2 = zeros(nbus - 1, nbus - 1);
   J3 = zeros(nbus - 1, nbus - 1);
   J4 = zeros(nbus - 1, nbus - 1);
    
  while( n <= maxiter)
    Scon = conj(V) .*(Ybus * V);
    inp = [angle(V(2: end));  abs(V(2: 2 + nPQ - 1))];
    Pcal = real(Scon);
    Qcal = -imag(Scon);
    delP = Psp - Pcal(2: end);
    delQ = Qsp - (Qcal(2: 2 + nPQ - 1 ));

    delPQ = [delP; delQ];
    
    if(abs(delPQ) < 0.001 )
      break;
    endif
    
    % All partial derivatives are calculated without considering the type
    % of buses but while constructing J only the essential rows and columns
    % are used
    for i = 2: nbus
      for k  = 2: nbus
        if(i == k)
          PQdii = -1j * Scon(i) + 1j * V(i) * conj(V(i)) * Ybus(i, i);
          PQvii = Scon(i) + V(i) * conj(V(i)) * Ybus(i, i);
          J1(i - 1, i - 1) = real(PQdii);
          J3(i - 1, i - 1) = -imag(PQdii);
          J2(i - 1, i - 1) = real(PQvii);
          J4(i - 1, i - 1) = -imag(PQvii);
        else
          PQdik = 1j * V(k) * conj(V(i)) * Ybus(i, k);
          PQvik = V(k) * conj(V(i)) * Ybus(i, k);
          J1(i - 1, k - 1) = real(PQdik);
          J3(i - 1, k - 1) = -imag(PQdik);
          J2(i - 1, k - 1) = real(PQvik);
          J4(i - 1, k - 1) = -imag(PQvik);
        endif    
      endfor
    endfor
    
    J = [J1(1: nbus - 1, 1: nbus - 1) J2(1: nbus - 1, 1: nPQ); ...
         J3(1: nPQ, 1: nbus - 1)      J4(1: nPQ, 1: nPQ)     ];

    correc = inv(J) * delPQ;
    correc(nPQ + nPV + 1 : end) = correc(nPQ + nPV + 1: end) .* abs(V(2: 2 + nPQ - 1));
    inp = inp + correc;
    
    V(2: 2 + nPQ - 1) = inp(nPQ + nPV + 1: end) .* exp(1j * inp(1: nPQ));
    
    if(nPV ~= 0)
      V(2 + nPQ: end) = abs(V(2 + nPQ: end)) .* exp(1j * inp(nPQ + 1: nPQ + nPV));
    endif
    
    n++;
  endwhile
  
 V = [SlackVoltage ; V(2: end)];
 S = V .* conj(Ybus * V) ;
 
 disp(n);
 
 endfunction
 
