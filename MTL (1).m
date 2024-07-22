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

W_h_1 = (8*exp(A))/(exp(2*A)-2)
W_h_2 = (2/pi)*(B-1-log(2*B-1)+((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r)))


if W_h_1<2
    W=h*W_h_1;
elseif W_h_2>2
    W=h*W_h_2;
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

plot(z, Vz, 'LineWidth',2); 
title('V_{iso}(z) Along Line'); grid on;
xlabel('z (m)'); ylabel('V_{iso}(z) (V)');







