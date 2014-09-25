%Initialize all constants
g_k = 36;
g_na = 120;
g_l = .3;
e_k = -12;
e_na = 115;
e_l = 10.6;
v_rest = -70;
c_m = 1;
I = 0;

%Introduce all variables
% v_m = membrane voltage
% m = probability of active sodium channels
% h = probability of inactive sodium channels
% n = probability of active potassium channels

v_m = 0; %start membrane at 0 for the purposes of solving the differential equations

step = .01;
t = 0:step:100;

%Initial conditions
a_m = .1*((25-v_m)/(exp((25-v_m)/10)-1));
b_m = 4*exp(-v_m/18);
a_h = .07*exp(-v_m/20);
b_h = 1/(exp((30-v_m)/10)+1);
a_n = .01*((10-v_m)/(exp((10-v_m)/10)-1));
b_n = .125*exp(-v_m/80);

m = a_m/(a_m + b_m);
h = a_h/(a_h + b_h);
n = a_n/(a_n + b_n);

G_na = (m^3)*g_na*h;
G_k = (n^4)*g_k;
G_l = g_l;

%main loop over 100 ms: calculate probability coeffs, ion channel currents,
%and run euler's method to obtain next point in time
for i = 1:length(t)-1

a_m(i) = .1*((25-v_m(i))/(exp((25-v_m(i))/10)-1));
b_m(i) = 4*exp(-v_m(i)/18);
a_h(i) = .07*exp(-v_m(i)/20);
b_h(i) = 1/(exp((30-v_m(i))/10)+1);
a_n(i) = .01*((10-v_m(i))/(exp((10-v_m(i))/10)-1));
b_n(i) = .125*exp(-v_m(i)/80);

i_na = (m(i)^3)*g_na*h(i)*(v_m(i) - e_na);
G_na(i+1) = (m(i)^3)*g_na*h(i);
i_k = (n(i)^4)*g_k*(v_m(i) - e_k);
G_k(i+1) = (n(i)^4)*g_k;
i_l = g_l*(v_m(i) - e_l);
G_l = g_l;
i_ion = I - i_na - i_k - i_l;

v_m(i+1) = v_m(i) + (step*(i_ion/c_m));
m(i+1) = m(i) + (step*((a_m(i)*(1-m(i)))-(b_m(i)*m(i))));
h(i+1) = h(i) + (step*((a_h(i)*(1-h(i)))-(b_h(i)*h(i))));
n(i+1) = n(i) + (step*((a_n(i)*(1-n(i)))-(b_n(i)*n(i))));

end

v_m = v_m + v_rest; %shift resting potential down to -70mV

plot(t,v_m);
figure
plot(t,G_k,'r',t,G_na,'b');



