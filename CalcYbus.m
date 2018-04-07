%data = input('Enter the network data (from, to, R, X, hlch): '); 

function Ybus = CalcYbus(data)
  
ftbus = data(:, 1:2);  

z = data(:, 3) + 1j * data(:, 4);   
y = 1 ./ z;                         
hl = data(:, 5) * 1j;

g = (min(min(ftbus)) == 0); %1 if ground node is present, 0 otherwise
nodes = max(max(ftbus)) + g; 
    
Ybus = zeros(nodes, nodes);  

for k = 1: length(ftbus(:,1))     
  
  ynew = y(k);   
  el = ftbus(k, :) + g; 
  Ybus(el, el) = Ybus(el, el) + ynew * [1 -1; -1 1] + hl(k) * eye(2);
  
end

disp(Ybus(1 + g: end, 1 + g: end));  

end