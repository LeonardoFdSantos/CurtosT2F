%% Teste de curto na realidade construida.

clc;
clear all; 
close all;

v=220*sqrt(3);
fp=[1]; %####ALTERAR####
S=[75e3];        %####ALTERAR####

%% Distancia 0.9km

%DIMENCIONAMENTO DAS INDUTANCIAS
f=60; %frequencia da rede       ####ALTERAR####
DI=8.01e-3;   %Di�metro do condutor [m]
RI=1.102;      %Resist�ncia el�trica m�xima CA 60Hz 75�C [Ohm/km]
rmgi=0.00308;  %Raio m�dio geom�trico [m]
Dotima=1.60;
ri=RI*1.609344; % transforma em ohm/milha 
rd=((pi)^2*f*10^(-4))/0.621371; %resistencia de terra ohm/milha 
GMRi=rmgi*3.28084; % transforma para P�S
dij=1.2; %m distancia entre os condutores i e j ####ALTERAR####
% escolhe dist�ncia entre cabos
dij=[ Dotima]; %m distancia entre os condutores i e j ####ALTERAR####
Dij=dij*3.28084; % transforma em P�S
h1=10.372; %distancia condutor mais alto da terra [m] ###ALTERAR####
h=10.372; %distancia condutor mais alto da terra [m] ###ALTERAR####
rti=10;
rtc=10;
%IMPEDANCIAS PR�PRIAS E MUTUAS
zp = ri+rd+(0.12134i)*(log(1/GMRi)+7.93402); %OHM/milha
Zp=zp*0.621371192; % tranfomada pra ohm/km
zm = rd+(0.12134i)*(log(1/Dij)+7.93402); %OHM/milha
Zm=zm*0.621371192; % tranfomada pra ohm/km
% %% Aten��o
% Zm=0; %%######%para teste da eq do CC
%DIMENCIONAMENTO DAS CAPACITANCIAS SHUNT
%Dois condutores
Ri=(DI/2)*3.28084; %raio cond i pés
Rj=(DI/2)*3.28084; %raio cond j pés
%Distancias entre os condutores e suas imagens para dada config de poste
sii=2*h; %m
Sii=sii*3.28084;
Sjj=Sii;
sij=sqrt(((dij)^2)+(sii^2)); %m
Sij=sij*3.28084;
Sji=Sij;
%Calculo Coeficientes de potencias proprios e mutuos
Pii= 11.17689*log(Sii/Ri)*10^6;
Pjj= 11.17689*log(Sjj/Rj)*10^6;
Pij= 11.17689*log(Sij/Dij)*10^6;
Pji= 11.17689*log(Sji/Dij)*10^6;
Mp=[Pii Pij; Pji Pjj];
C = inv(Mp); %capacitancias shunt
C=C*0.621371192; % cconvertendo milha para km 
%capacitancias de balanceamento transversal
Cct=C(1,1)+C(2,1); %C entre condutor e terra Farad/km
Ccc=-C(2,1); %C entre condutores F/km
Ceq=Cct-Ccc; %capacitancia de equalização F/km
%IMPEDANCIA DE COMPENSA��O
Ze = Zp-2*Zm; %OHM/km
%IMPEDANCIAS S�RIE
Zl = Zp-Zm; %OHM/km
dd = [0.9];

for compensada=[0 1 2]
%fazer tudo em fun��o da distancia
    for z = 1:length(dd) % Dist�ncia em KM   
    d=dd(z);
    cequa1=Ceq*d;
    ccc1=Ccc*d;
    cct1=Cct*d;
    rp=real(Zp)*d;
    rm=real(Zm)*d;
    re=(real(Ze)*d)-rti-rtc;
        if re<0
            NRE=1e-12;
        else
            NRE=re;
        end

    % rg=real(Zg)*d;
    lp=(imag(Zp)/(2*pi*f))*d;
    lm=(imag(Zm)/(2*pi*f))*d;
    le=1e-10;%(imag(Ze)/(2*pi*f))*d;
    % lg=(imag(Zg)/(2*pi*f))*d;

        if compensada==0
        cequa=0; 
        NRE=1e-10;
        le=1e-10;
        %Sem Compensa��o 
        elseif compensada==1
        % cequa=Ceq*d; 
        NRE=1e-12;
        %Compensa��o C
        elseif compensada==2
        %Compensa��o RC
        end
    end
end

rm1 = rm;
lp1 = lp;
lm1 = lm;
rp1 = rp;

%% Distancia 0.5km

%DIMENCIONAMENTO DAS INDUTANCIAS
f=60; %frequencia da rede       ####ALTERAR####
DI=8.01e-3;   %Di�metro do condutor [m]
RI=1.102;      %Resist�ncia el�trica m�xima CA 60Hz 75�C [Ohm/km]
rmgi=0.00308;  %Raio m�dio geom�trico [m]
Dotima=1.60;
ri=RI*1.609344; % transforma em ohm/milha 
rd=((pi)^2*f*10^(-4))/0.621371; %resistencia de terra ohm/milha 
GMRi=rmgi*3.28084; % transforma para P�S
dij=1.2; %m distancia entre os condutores i e j ####ALTERAR####
% escolhe dist�ncia entre cabos
dij=[ Dotima]; %m distancia entre os condutores i e j ####ALTERAR####
Dij=dij*3.28084; % transforma em P�S
h1=10.372; %distancia condutor mais alto da terra [m] ###ALTERAR####
h=10.372; %distancia condutor mais alto da terra [m] ###ALTERAR####
rti=10;
rtc=10;
%IMPEDANCIAS PR�PRIAS E MUTUAS
zp = ri+rd+(0.12134i)*(log(1/GMRi)+7.93402); %OHM/milha
Zp=zp*0.621371192; % tranfomada pra ohm/km
zm = rd+(0.12134i)*(log(1/Dij)+7.93402); %OHM/milha
Zm=zm*0.621371192; % tranfomada pra ohm/km
% %% Aten��o
% Zm=0; %%######%para teste da eq do CC
%DIMENCIONAMENTO DAS CAPACITANCIAS SHUNT
%Dois condutores
Ri=(DI/2)*3.28084; %raio cond i pés
Rj=(DI/2)*3.28084; %raio cond j pés
%Distancias entre os condutores e suas imagens para dada config de poste
sii=2*h; %m
Sii=sii*3.28084;
Sjj=Sii;
sij=sqrt(((dij)^2)+(sii^2)); %m
Sij=sij*3.28084;
Sji=Sij;
%Calculo Coeficientes de potencias proprios e mutuos
Pii= 11.17689*log(Sii/Ri)*10^6;
Pjj= 11.17689*log(Sjj/Rj)*10^6;
Pij= 11.17689*log(Sij/Dij)*10^6;
Pji= 11.17689*log(Sji/Dij)*10^6;
Mp=[Pii Pij; Pji Pjj];
C = inv(Mp); %capacitancias shunt
C=C*0.621371192; % cconvertendo milha para km 
%capacitancias de balanceamento transversal
Cct=C(1,1)+C(2,1); %C entre condutor e terra Farad/km
Ccc=-C(2,1); %C entre condutores F/km
Ceq=Cct-Ccc; %capacitancia de equalização F/km
%IMPEDANCIA DE COMPENSA��O
Ze = Zp-2*Zm; %OHM/km
%IMPEDANCIAS S�RIE
Zl = Zp-Zm; %OHM/km
dd = [0.5];

for compensada=[0 1 2]
%fazer tudo em fun��o da distancia
    for z = 1:length(dd) % Dist�ncia em KM   
    d=dd(z);
    cequa2=Ceq*d;
    ccc2=Ccc*d;
    cct2=Cct*d;
    rp=real(Zp)*d;
    rm=real(Zm)*d;
    re=(real(Ze)*d)-rti-rtc;
        if re<0
            NRE=1e-12;
        else
            NRE=re;
        end

    % rg=real(Zg)*d;
    lp=(imag(Zp)/(2*pi*f))*d;
    lm=(imag(Zm)/(2*pi*f))*d;
    le=1e-10;%(imag(Ze)/(2*pi*f))*d;
    % lg=(imag(Zg)/(2*pi*f))*d;

        if compensada==0
        cequa=0; 
        NRE=1e-10;
        le=1e-10;
        %Sem Compensa��o 
        elseif compensada==1
        % cequa=Ceq*d; 
        NRE=1e-12;
        %Compensa��o C
        elseif compensada==2
        %Compensa��o RC
        end
    end
end

rm2 = rm;
lp2 = lp;
lm2 = lm;
rp2 = rp;

%% Começa Rotina de Curtos Fim de linha
valores_resultados_simulacao_fim_linha_com_comp(1, :) = string({'n' 'tipoCurto' 'Raf1' 'Rbf1' 'Rcf1' 'Raf2' 'Rbf2' 'RcF2' 'IA_Sistema' 'IB_Sistema' 'IC_Sistema'});
valores_resultados_simulacao_meio_linha_com_comp_1(1, :) = string({'n' 'm1' 'tipoCurto' 'Raf1' 'Rbf1' 'Rcf1' 'Raf2' 'Rbf2' 'RcF2' 'IA_Sistema' 'IB_Sistema' 'IC_Sistema'});
valores_resultados_simulacao_meio_linha_com_comp_2(1, :) = string({'n' 'm2' 'tipoCurto' 'Raf1' 'Rbf1' 'Rcf1' 'Raf2' 'Rbf2' 'RcF2' 'IA_Sistema' 'IB_Sistema' 'IC_Sistema'});

TempoSimulacao = 5;
FatorSimulacao = 1e-3;

rtc1 = 8;
rtc2 = 8;

b = 0;
c = 1;
e = 1;
Parametros_testes = [0.001 .1 .2 .3 .4 .5 .6 .7 .8 .9 .999];
for b = [1:1:10]
    if b == 1
        Raf1 = 1e-5;
        Rbf1 = 1e-5;
        Rcf1 = 1e-5;
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        tipoCurto = b;
    elseif b == 2
        Raf1 = 1e-5;
        Rbf1 = inf;
        Rcf1 = 1e-5;
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        tipoCurto = b;
    elseif b == 3
        Raf1 = inf;
        Rbf1 = 1e-5;
        Rcf1 = 1e-5;
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        tipoCurto = b;
    elseif b == 4
        Raf1 = 1e-5;
        Rbf1 = 1e-5;
        Rcf1 = inf;
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        tipoCurto = b;
    elseif b == 5
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        tipoCurto = b;
    elseif b == 6
        Raf2 = 1e-5;
        Rbf2 = 1e-5;
        Rcf2 = 1e-5;
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        tipoCurto = b;
    elseif b == 7
        Raf2 = 1e-5;
        Rbf2 = inf;
        Rcf2 = 1e-5;
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        tipoCurto = b;
    elseif b == 8
        Raf2 = inf;
        Rbf2 = 1e-5;
        Rcf2 = 1e-5;
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        tipoCurto = b;
    elseif b == 9
        Raf2 = 1e-5;
        Rbf2 = 1e-5;
        Rcf2 = inf;
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        tipoCurto = b;
    elseif b == 10
        Raf2 = inf;
        Rbf2 = inf;
        Rcf2 = inf;
        Raf1 = inf;
        Rbf1 = inf;
        Rcf1 = inf;
        tipoCurto = b;
    end
    sim('.\GhendyanoT2F_SemFuse.slx')
    CorrenteSistemaT2F_lista = abs(CorrenteSistemaT2F)/sqrt(2);
    valores_resultados_simulacao_fim_linha_com_comp(c+1, :) = [c tipoCurto Raf1 Rbf1 Rcf1 Raf2 Rbf2 Rcf2 CorrenteSistemaT2F_lista];
    c = c + 1;
end
fprintf('Terminado fim de Linha! \n');
%% Meio Linha 1
c = 1;
for m1 = Parametros_testes
    for b = [1:1:10]
        if b == 1
            Raf1 = 1e-5;
            Rbf1 = 1e-5;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 2
            Raf1 = 1e-5;
            Rbf1 = inf;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 3
            Raf1 = inf;
            Rbf1 = 1e-5;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 4
            Raf1 = 1e-5;
            Rbf1 = 1e-5;
            Rcf1 = inf;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 5
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 6
            Raf2 = 1e-5;
            Rbf2 = 1e-5;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 7
            Raf2 = 1e-5;
            Rbf2 = inf;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 8
            Raf2 = inf;
            Rbf2 = 1e-5;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 9
            Raf2 = 1e-5;
            Rbf2 = 1e-5;
            Rcf2 = inf;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 10
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        end
        sim('.\GhendyanoT2FMeioLinha1_SemFuse.slx')
        CorrenteSistemaT2F_lista = abs(CorrenteSistemaT2F)/sqrt(2);
        valores_resultados_simulacao_meio_linha_com_comp_1(c+1, :) = [c m1 tipoCurto Raf1 Rbf1 Rcf1 Raf2 Rbf2 Rcf2 CorrenteSistemaT2F_lista];
        c = c + 1;
    end
end

%% Meio de Linha 2
c = 1;
for m2 = Parametros_testes
    for b = [1:1:10]
        if b == 1
            Raf1 = 1e-5;
            Rbf1 = 1e-5;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 2
            Raf1 = 1e-5;
            Rbf1 = inf;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 3
            Raf1 = inf;
            Rbf1 = 1e-5;
            Rcf1 = 1e-5;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 4
            Raf1 = 1e-5;
            Rbf1 = 1e-5;
            Rcf1 = inf;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 5
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            tipoCurto = b;
        elseif b == 6
            Raf2 = 1e-5;
            Rbf2 = 1e-5;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 7
            Raf2 = 1e-5;
            Rbf2 = inf;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 8
            Raf2 = inf;
            Rbf2 = 1e-5;
            Rcf2 = 1e-5;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 9
            Raf2 = 1e-5;
            Rbf2 = 1e-5;
            Rcf2 = inf;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        elseif b == 10
            Raf2 = inf;
            Rbf2 = inf;
            Rcf2 = inf;
            Raf1 = inf;
            Rbf1 = inf;
            Rcf1 = inf;
            tipoCurto = b;
        end
        sim('.\GhendyanoT2FMeioLinha2_SemFuse.slx')
        CorrenteSistemaT2F_lista = abs(CorrenteSistemaT2F)/sqrt(2);
        valores_resultados_simulacao_meio_linha_com_comp_2(c+1, :) = [c m2 tipoCurto Raf1 Rbf1 Rcf1 Raf2 Rbf2 Rcf2 CorrenteSistemaT2F_lista];
        c = c + 1;
    end
end

writematrix(valores_resultados_simulacao_fim_linha_com_comp, 'ResultadoFimLinha_SemFuse.csv');
writematrix(valores_resultados_simulacao_meio_linha_com_comp_1, 'ResultadoMeioLinha1_SemFuse.csv');
writematrix(valores_resultados_simulacao_meio_linha_com_comp_2, 'ResultadoMeioLinha2_SemFuse.csv');

fprintf('Terminado meio de Linha! \n');