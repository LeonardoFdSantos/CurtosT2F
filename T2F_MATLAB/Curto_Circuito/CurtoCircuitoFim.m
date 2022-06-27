%% Dados do Alimentador
Ze=13.8^2/120*exp(i*atan(7));%% SCC=120 MVA
Zbase=13.8^2/0.120;
Ztri=2*(0.0048119 +i*0.018511)*Zbase;

%% Dados T2F
Za=Ze+Ztri+rp+i*2*pi*60*lp+Raf;
Zb=Ze+Ztri+rp+i*2*pi*60*lp+Rbf;
Zc=Ze+Ztri+rti+NRE+i*2*pi*60*le+Rcf; %% para curto-circuitos antes do TR do consumidor
Zm=rm+i*2*pi*60*lm;

%% Tensões
Va=13800/sqrt(3)*exp(0);
Vb=13800/sqrt(3)*exp(-i*2*pi/3);
Vc=13800/sqrt(3)*exp(i*2*pi/3);

I1=(Va-Vb-(Vb-Vc)*(Zm-Zb)/(Zb+Zc))/(Za+Zb-2*Zm-(Zm-Zb)^2/(Zb+Zc));
I2=(Vb-Vc-(-Zb+Zm)*I1)/(Zb+Zc);

IA=abs(I1);
IB=abs(I2-I1);
IC=abs(-I2);

