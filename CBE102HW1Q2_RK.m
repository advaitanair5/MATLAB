h = 1e-3; %increments of time
t=0:h:10; %time elapsed
L1 = 1; %piston bottom to cylinder bottom initial distance 
x0 = L1; %initial position
v0 = 0; %initial velocity
R = 8.314; %gas constant (J/(mol*K))
T = 310; %Temperature of gas (k)
n = 0.323; %moles of gas (moles)
m = 1.5; %Mass of piston 
Patm = 101325; 
x=zeros(1,length(t));
x(1)=x0;
v=zeros(1,length(t));
v(1)=v0;
M = 9; 
mu = 1450; %change to 0 for reversible!)
r = 0.05; %Cylinder radius (m)
L2 = 0.03; %piston thickness
g = 9.81; %gravity

Tc = 33.19; %K
Pc = 1.313*10^6; %pa
phi = 0.42748;
Omega = 0.08664;
omega = -0.216;
sig = 1;
eps = 0;

Tr = T/Tc;
alpha = (Tr)^(-1/2);

% change the function depending on the system

a = phi*alpha*R^2*Tc^2/Pc; %a(T)
b = Omega*(R*Tc/Pc);

% change the function depending on the system
P =@(x) ((R*T)./(((pi*r^2.*x)/n) -b)) - (a/((sig-eps)*b))*((1./(pi*r^2.*x) +eps*b) - (1./((pi*r^2.*x)+sig*b)));
F_tx=@(x,v) v;
F_tv=@(x,v) -(mu*2*pi*r*L2/(M+m)).*v-g-((Patm-P(x)).*pi*r^2/(M+m));

% calculation loop
for i=1:(length(t)-1)
    
    kx1 = F_tx(x(i),v(i))*h;
    kv1 = F_tv(x(i),v(i))*h;
    kx2 = F_tx(x(i)+kx1/2,v(i)+kv1/2)*h;
    kv2 = F_tv(x(i)+kx1/2,v(i)+kv1/2)*h;
    kx3 = F_tx((x(i)+kx2/2),v(i)+kv2/2)*h;
    kv3 = F_tv((x(i)+kx2/2),(v(i)+kv2/2))*h;
    kx4 = F_tx((x(i)+kx3),v(i)+kv3)*h;
    kv4 = F_tv((x(i)+kx3),(v(i)+kv3))*h;
    
    x(i+1) = x(i) + (kx1+2*kx2+2*kx3+kx4)/6; 
    v(i+1) = v(i) + (kv1+2*kv2+2*kv3+kv4)/6;

end

%displacement vs time
figure;
plot(t, x);
xlabel('Time Elapsed [s]');
ylabel('Displacement [m]');
title('Displacement (RK) vs. Time Elapsed');
grid on; 

%velocity vs time
figure;
plot(t, v);
xlabel('Time Elapsed [s]');
ylabel('Velocity [m/s]');
title('Velocity (RK) vs. Time Elapsed');
grid on; 

%acceleration vs time
accel = F_tv(x,v);
figure;
plot(t, accel);
xlabel('Time Elapsed [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration (RK) vs. Time Elapsed');
grid on; 

volume = (pi*r^2.*x);
figure;
plot(t, volume);
xlabel('Time Elapsed [s]');
ylabel('Volume [m^3]');
title('Gas Volume (RK) vs. Time Elapsed');
grid on; 

Pres = ((R*T)./(((pi*r^2.*x)/n) -b)) - (a/((sig-eps)*b))*((1./(pi*r^2.*x) +eps*b) - (1./((pi*r^2.*x)+sig*b)));
figure;
plot(t, Pres);
xlabel('Time Elapsed [s]');
ylabel('Pressure [Pa]');
title('Gas Pressure (RK) vs. Time Elapsed');
grid on; 

Wg = (M+m)*g*(L1-x);
figure;
plot(t, Wg);
xlabel('Time Elapsed [s]');
ylabel('Work [J]');
title('Gravity Work (RK) vs. Time Elapsed');
grid on; 

Watm = Patm*pi*r^2*(L1-x);
figure;
plot(t, Watm);
xlabel('Time Elapsed [s]');
ylabel('Work [J]');
title('Atmospheric Work (RK) vs. Time Elapsed');
grid on; 

W_fr = -cumtrapz(t, mu*2*pi*r*L2*v.^2);
figure;
plot(t, W_fr);
xlabel('Time Elapsed [s]');
ylabel('Work [J]');
title('Friction Work (RK) vs. Time Elapsed');
grid on; 

integrand = (R * T) ./ ((pi * r^2 .* x) / n - b) - ...
            (a ./ (((pi * r^2 .* x) / n + eps * b)) .* ((pi * r^2 .* x / n + sig * b)));
W_gas = pi*r^2*cumtrapz(x, integrand);
figure;
plot(t, W_gas);
xlabel('Time Elapsed [s]');
ylabel('Work [J]');
title('Gas Work (RK) vs. Time Elapsed');
grid on; 

tot_work = W_fr+W_gas+Watm+Wg;
KE = 0.5*(M+m).*(v.^2);
figure;
plot(t,tot_work, 'b');
hold on;
plot(t, KE, 'r');
hold off;
xlabel('Time Elapsed [s]');
ylabel('[J]');
legend('Total Work', 'Kinetic Energy');
title('Total Work & KE (RK) vs. Time Elapsed');

[M1_RK, x_real_RK, W_gas_RK] = Quasi_work(phi, Omega, sig, eps, alpha);

figure;
plot(x_real_RK, W_gas_RK);
xlabel('Displacement (m)');
ylabel('Work (J)');
title('Quasi Static Work (RK) vs Displacement');
grid on;
