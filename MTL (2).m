%Part 01
e_r = 9.8; 
h = 1.6e-3; 
Z_inf = 50;
V_iso = 10e-3;
I_iso = 0.4e-3;
f = 10e9;
L = 100e-3;
z = 0:L/1000:L;
c = 3*1e8;

A = (Z_inf/60)*sqrt((e_r+1)/2)+((e_r-1)/(e_r+1))*(0.23*+(0.11/e_r));
B = (377*pi)/(2*Z_inf*sqrt(e_r));

W_h_1 = (8*exp(A))/(exp(2*A)-2);
W_h_2 = (2/pi)*(B-1-log(2*B-1)+((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r)));


if W_h_1<2
    W=h*W_h_1
elseif W_h_2>2
    W=h*W_h_2
end

Z0 = V_iso/I_iso;
Zin = Z0;

Gamma_In = (Zin-Z_inf)/(Zin+Z_inf);
e_eff= (e_r+1)/2 + (e_r-1)/2 *1/sqrt(1+12*h/W);
Vp = c/sqrt(e_eff);
K = 2*pi*f/Vp;
Gamma_L = Gamma_In*exp(2*i*K*L);
ZL = Z_inf*(1+Gamma_In)/(1-Gamma_In);
Vz = V_iso .* (1+ (Gamma_In*exp(2*i*K.*z))) .*exp(2*i*K.*z)./(1+Gamma_In);
Vz = abs(Vz);

figure(1)
plot(z, Vz, 'LineWidth',2); 
title('V_{iso}(z) Along Line'); grid on;
xlabel('z (m)'); ylabel('V_{iso}(z) (V)');



%Part_02
V1_0 = V_iso;
I1_0 = I_iso;
V2_0 = 0;
I2_0 = 0;
s = 10e-3;

A = (Z_inf/60)*sqrt((e_r+1)/2)+((e_r-1)/(e_r+1))*(0.23*+(0.11/e_r));
B = (377*pi)/(2*Z_inf*sqrt(e_r));

W_h_1 = (8*exp(A))/(exp(2*A)-2);
W_h_2 = (2/pi)*(B-1-log(2*B-1)+((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r)));

if W_h_1<2
    W=h*W_h_1
elseif W_h_2>2
    W=h*W_h_2
end

e_eff= (e_r+1)/2 + (e_r-1)/2 *1/sqrt(1+12*h/W);
Vp = c/sqrt(e_eff);
K = 2*pi*f/Vp;
Theta = K*L;
E = exp(-2*1i*K*L);

A= [1 1;1 -1];
B = [V1_0;V2_0];
C= [I1_0;I2_0];

Vm0 = 1/sqrt(2)*A*B;
Im0= 1/sqrt(2)*A*C;

% From the Given website, we have the following:
Z_even = 45.5421;
Z_odd = 44.4239;
k_even = 6.80943;
k_odd = 6.58146;

%Using the above value, we calculate:
Ke = K*sqrt(k_even);
Ko = K*sqrt(k_odd);

% For Even Line
Vo_even = Vm0(1);
Io_even = Im0(1);
Vo_even_plus = 0.5.*(Vo_even + Z_even.*Io_even);
Vo_even_minus = 0.5.*(Vo_even - Z_even.*Io_even);
V_even = Vo_even_plus*exp(-1i*Ke.*z)+Vo_even_minus*exp(1i*Ke.*z);

% For Odd Line
Vo_odd = Vm0(2);
Io_odd = Vm0(2);
Vo_odd_plus = 0.5.*(Vo_odd + Z_odd.*Io_odd);
Vo_odd_minus = 0.5.*(Vo_odd - Z_odd.*Io_odd);
V_odd = Vo_odd_plus*exp(-1i*Ko.*z)+Vo_odd_minus*exp(1i*Ko.*z); 

V = (1/sqrt(2)).*A*[V_even; V_odd];
V1 = abs(V(1,:));
V2 = abs(V(2,:));

figure(2);
plot(z,V1, 'LineWidth',2);
title('V_{1}(z) and V_{2}(z) for s = 10mm');
grid on; xlabel('z (m)'); ylabel('V(z)');
hold on;
plot(z, V2, 'LineWidth',2);
legend ('V_1(z)','V_2(z)');
