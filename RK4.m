function dclear = RK4(data2)
  
  %data2 = [.4605, 0, 1.8, 0.65, 1.3, 0.8, 5, 60, 0, 0, 0.2, 0.5, .01];

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
   
   
 end