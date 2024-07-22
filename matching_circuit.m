%% Problem nr.1
Z0 = 50;
ZL1 = 100;
ZL2 = 200;
Zc1 = sqrt(Z0*ZL1); % ZC1 = 70.7107
Zc2 = sqrt(Z0*ZL2); % ZC2 = 100
c = 3e8;
F0 = 5e9;
f = 0:10e6:10e9; %Range from [0, 2*f0]

lambda0 = c./F0;
L = lambda0/4;
lamda = c./f;
K = (2*pi)./(lamda);

GammaB1 = (ZL1-Zc1)/(ZL1+Zc1);
GammaA1 = GammaB1*exp(-2i*K*L);
ZA1 = Zc1*(1+GammaA1)./(1-GammaA1);
GammaA_minus1 = (ZA1 -Z0)./(ZA1 +Z0);

GammaB2 = (ZL2-Zc2)/(ZL2+Zc2);
GammaA2 = GammaB2*exp(-2i*K*L);
ZA2 = Zc2*(1+GammaA2)./(1-GammaA2);
GammaA_minus2 = (ZA2 -Z0)./(ZA2 +Z0);

plot(f/F0, abs(GammaA_minus1), 'linewidth', 2); hold on;
plot(f/F0, abs(GammaA_minus2), 'linewidth', 2); hold off;
yline(0.1, 'LineWidth', 2); % 20log(S_{11}(f)) < -20dB
legend('ZL1 = 100 ohms and ZC1 = 70.7107',...
'ZL2 = 200 ohms and ZC2 = 100', '20log(S_{11}(f)) < -20dB');
grid on; xlabel('Frequency f/f0'); ylabel('S_{11}(f)');
title('S_{11}(f) as function of f/f0')


%% Problem nr.2
Z0 = 50;
ZR1 = 50;
ZL = 100;
ZR2 = 100;
ZC1 = 60;
ZC2 = 84;
c = 3e8;
F0 = 5e9;
f = 0:10e6:10e9; %Range from [0, 2*f0]

lambda0 = c./F0;
lamda = c./f;
L1 = lambda0/4;
K1 = (2*pi)./(lamda);
L2 = lambda0/4;
K2 = (2*pi)./(lamda);

%For Section L1
A1 = cos(K1*L1);
B1 = 1i*ZC1*(sin(K1*L1));
C1 = 1i*(1/ZC1)*(sin(K1*L1));
D1 = cos(K1*L1);

%For Section L2
A2 = cos(K2*L2);
B2 = 1i*ZC2*(sin(K2*L2));
C2 = 1i*(1/ZC2)*(sin(K2*L2));
D2 = cos(K2*L2);

Atot = A1.*A2+B1.*C2;
Btot = A1.*B2+B1.*D2;
Ctot = C1.*A2+D1.*C2;
Dtot = C1.*B2+D1.*D2;
ABCD_Total = [Atot Btot; Ctot Dtot];
Total = [A1 B1; C1 D1].*[A2 B2; C2 D2];

S_11_num = (ZR2.*Atot+Btot-ZR1.*ZR2.*Ctot-ZR1.*Dtot);
S_11_den = (ZR2.*Atot+Btot+ZR1.*ZR2.*Ctot+ZR1.*Dtot);
S_11 = S_11_num./S_11_den;

plot(f/F0, abs(S_11), 'linewidth', 2);
yline(0.1, 'LineWidth', 2); % 20log(S_{11}(f)) < -20dB
grid on; xlabel('Frequency f/f0'); ylabel('S_{11}(f)');
title('S_{11}(f) as function of f/f0')
