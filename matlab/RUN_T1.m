%% Solu��o de um problema de otimiza��o via LINPROG
% Victor F. Carvalho     25/10/2023
%%
clear all;  
close all;  
clc
tic

global DGER DLIN DEMANDA SW h delta Sbase

% DADOS DE GERA��O
     % BARRA  CUSTO($/MWh)  MAX(MW) MIN(MW)  RAMPA(MW)    C02(m3/MWh) TIPO DEFICIT($/MWh)
DGER=[  1        10             30      05       10            90      0      200
        2        30             40      15       05            10      1      600
        3       100             40      00       03            70      2     2000];

%###########################################
% fiquem a vontade para usar outros valores:    
h=0.25; % fator de convers�o emiss�o - $/MW 
delta=0.9; % fator de pondera��o emiss�o     
%##########################################

% DADOS DA REDE
      %  DE   PARA   SUSCEPT�NCIA(OHMS) CONDUT�NCIA(OHMS) LIMITES (MW)  
 DLIN=[   1     2          33                25                20;
          1     3          50                20                25;
          1     3          50                20                25;
          2     3          50                20                30];
           
      
 % DADOS DE DEMANDA
           %HORA BARRA1  BARRA2  BARRA3
           %  h   (MW)    (MW)    (MW)   
 DEMANDA= [   1     0      1.5*40      1.5*30;
              2     0      1.5*43      1.5*25;
              3     0      1.5*25      1.5*25;
              
              0     1       2       3]; % �ltima linha refere-se � posi��o

DEMANDA(1:end-1,2:end)= 1.0*  DEMANDA(1:end-1,2:end)  ;    
          
NL = length(DLIN(:,1)); % N�mero de Circuitos      
NG = length(DGER(:,1)); % N�mero de Geradores 




%% Inser��o de perdas por meio de uma constante Pperdaskm = gkm(Tetak - Tetam) * CP

          % Perdas k - m   % Perdas m - k
CP =[[DLIN(:,1:2);[DLIN(:,2) DLIN(:,1)]] 0.005*[DLIN(:,5); DLIN(:,5)]]; 

%% OP��ES
itermax = 100000;
%--------------------------------------------------------------------------
B_=[];
Beq=[];
lb=[];
ub=[];
f=[];
fcost=[];
femission=[];
A=[];
Aeq=[];
MRampa=[];
Rampa=[];

for TimeD = 1 : length(DEMANDA(:,1))-1

%% FLOW DC
DataDC_T1; % Obten��o da B'

%--------------------------------------------------------------------------       

%% LIMITES DAS VARI�VEIS

%       PG      |     D�ficit     |     Teta     | 
%% Desilgualdade de Gera��o

LGM = DGER(:,3)'; % Limite superior de gera��o
LGm = DGER(:,4)'; % Limite inferior de gera��o

%--------------------------------------------------------------------------
%% Determina��o dos D�ficits


DefS = PD'*Sbase;  % Limite Superior de D�ficit
DefI = zeros(1,length(DEMANDA(1,2:end))); % Limite Inferior de D�ficit


%--------------------------------------------------------------------------
%% Limite para os �ngulos

VAngLS = pi* ones(1,NG-1); % Vetor para os limites superiores dos �ngulos diferentes de SW

VAngLI = -pi* ones(1,NG-1); % Vetor para os limites inferiores dos �ngulos diferentes de SW


%--------------------------------------------------------------------------
%% RESTRI�OES DE DESIGUALDADES LINEARES  (A*x <= b)

%% Determina��o dos Fluxos

MatFluxoM_KM=zeros(NL,NG);
MatFluxom_KM=zeros(NL,NG);

MatFluxoM_MK=zeros(NL,NG);
MatFluxom_MK=zeros(NL,NG);

LimS_KM=zeros(1,NL);
LimI_KM=zeros(1,NL);

LimS_MK=zeros(1,NL);
LimI_MK=zeros(1,NL);

infoMaior_KM=zeros(NL,2);
infoMenor_KM=infoMaior_KM;

infoMaior_MK=zeros(NL,2);
infoMenor_MK=infoMaior_MK;


for i = 1 : NL
   % Fluxo K - M     fik <= F <= Fik
    infoMaior_KM(i,1:2)=[DLIN(i,1) DLIN(i,2)];

    MatFluxoM_KM(i,DLIN(i,1)) = BbusLin(DLIN(i,1),DLIN(i,2));
    MatFluxoM_KM(i,DLIN(i,2)) = -BbusLin(DLIN(i,1),DLIN(i,2));
    
    infoMenor_KM(i,1:2)=[DLIN(i,1) DLIN(i,2)];
    MatFluxom_KM(i,DLIN(i,1)) = -BbusLin(DLIN(i,1),DLIN(i,2));
    MatFluxom_KM(i,DLIN(i,2)) = BbusLin(DLIN(i,1),DLIN(i,2)); 
    
    LimS_KM(1,i) = DLIN(i,5)/Sbase;
    LimI_KM(1,i) = DLIN(i,5)/Sbase;
    
    
    % Fluxo M - K    
    infoMaior_MK(i,1:2)=[DLIN(i,2) DLIN(i,1)];

    MatFluxoM_MK(i,DLIN(i,2)) = BbusLin(DLIN(i,2),DLIN(i,1));
    MatFluxoM_MK(i,DLIN(i,1)) = - BbusLin(DLIN(i,2),DLIN(i,1));
    
    infoMenor_MK(i,1:2)=[DLIN(i,2) DLIN(i,1)];
    MatFluxom_MK(i,DLIN(i,2)) = -BbusLin(DLIN(i,2),DLIN(i,1));
    MatFluxom_MK(i,DLIN(i,1)) = BbusLin(DLIN(i,2),DLIN(i,1)); 
    
    LimS_MK(1,i) = DLIN(i,5)/Sbase;
    LimI_MK(1,i) = DLIN(i,5)/Sbase;
        
end

fluxos=[[infoMaior_KM;infoMaior_MK];[infoMenor_KM;infoMenor_MK]];

Mflu_Original=[[MatFluxoM_KM;MatFluxoM_MK];[MatFluxom_KM;MatFluxom_MK]];

MatFluxo=Mflu_Original;

MatFluxo(:,SW)=[];
PGflu=zeros(length(MatFluxo(:,1)),NG);
PDeflu=zeros(length(MatFluxo(:,1)),NG);

%         |   PG  |    D�ficit  |   Teta   | 

Desig_Fluxo=[ PGflu  , PDeflu , MatFluxo];
LimF=[LimS_KM, LimS_MK, LimI_KM, LimI_MK];

clear PGflu PDeflu infoMaior_KM infoMaior_MK infoMenor_KM infoMenor_MK MatFluxom_KM MatFluxoM_KM
clear LimS_KM LimS_MK LimI_KM LimI_MK MatFluxom_MK MatFluxoM_MK


%--------------------------------------------------------------------------
%% Determina��o das Limita��es de Rampa

%    |  PG  |    D�ficit   |    Teta    |        -P  <= PG(h+1)-PG(h) <= P
M1 = [zeros(NG) zeros(NG)   zeros(NG,NG-1)];

M2 = [-eye(NG)  zeros(NG)   zeros(NG,NG-1)];

M3 = [eye(NG)   zeros(NG)   zeros(NG,NG-1)];

MATRIZ1=[];
for k = 1 : length(DEMANDA(:,1))-1
    if k==TimeD
        MATRIZ1=[MATRIZ1,M2];
    elseif k == TimeD+1
        MATRIZ1=[MATRIZ1,M3];
    else
        MATRIZ1=[MATRIZ1,M1];
    end
end
if TimeD <= length(DEMANDA(:,1))-2
    Rampa=[Rampa,DGER(:,5)'];
    MRampa=[MRampa;MATRIZ1];
end

clear M1 M2 M3
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% Restri��es de Igualdade Lineares  (Aeq*x = beq)


%% Restri��o de Igualdade PGi + PDefi - sum(Pik) + Pperdasi = PDi

%-------- Fluxo em cada linha
SomaFluxo=zeros(NG,NG-1);
for i = 1 : NG
    for k = 1 : length(fluxos(:,1))/2
        if fluxos(k,1) == DGER(i,1)
            SomaFluxo(i,:)=SomaFluxo(i,:) - MatFluxo(k,:);
        end 
    end 
end

%-------- Perdas em cada linha                  

SomaPerdas=zeros(NG,NG-1);
for i = 1 : NG
    for k = 1 : length(fluxos(:,1))/2
        if fluxos(k,1) == DGER(i,1)
   
            % PikTotal      =  Pik           +     CP/2 *  gik *( Bik *(Tetai-Tetak)/Bik)

            SomaPerdas(i,:)= SomaPerdas(i,:) - (CP(k,3)/2)*Gbus(fluxos(k,1),fluxos(k,2))*(MatFluxo(k,:)./BbusLin(fluxos(k,1),fluxos(k,2)));
        end 
    end 
end

%--------------------------------------------------------
% Inseridas as Perdas
%           |   PG  |       D�ficit         |   Teta   | 
MatIgualdade=[eye(NG)/Sbase eye(NG)/Sbase SomaFluxo+SomaPerdas];
Big1=PD';

%% Somat�rio das Pot�ncias com os D�ficits sum(PG)+sum(Pdef)+sum(Pperdas)=PD


M=[ones(1,NG) ones(1,NG) sum(SomaPerdas)];
Big2=sum(PD)*Sbase;%sum(DEMANDA(1,2:end));

clear SomaFluxo MatPG Matfluxos


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constru��o da Matriz de Igualdades e Desigualdades;

[Lines,Colums]=size([M; MatIgualdade]);
[Lines2,Colums2]=size(Desig_Fluxo);
A1=[];
Aeq1=[];
for k=1:length(DEMANDA(:,1))-1
    if k == TimeD
        Aeq1=[Aeq1,[MatIgualdade;M ]];
        A1=[A1, Desig_Fluxo];
    else 
        Aeq1=[Aeq1,zeros(Lines,Colums)];
        A1=[A1, zeros(Lines2,Colums2)];
    end


end
A=[A;A1];
Aeq=[Aeq;Aeq1];

%% FOB
%                          PG                    |       PDEf        |    Teta
    f = [f,delta*DGER(:,6)'+(1-delta)*h*DGER(:,2)'  delta*DGER(:,8)' zeros(1,length(DGER(:,1))-1)]; % Fun��o Objetivo total
    
    
    fcost = [fcost,DGER(:,6)',DGER(:,8)',zeros(1,length(DGER(:,1))-1)];
    femission = [femission,DGER(:,2)',DGER(:,8)',zeros(1,length(DGER(:,1))-1)]; 
    
    Beq=[Beq,Big1,Big2];
    B_=[B_,LimF];
    lb=[lb,LGm,DefI,VAngLI];
    ub=[ub,LGM,DefS,VAngLS];
end
A=[A;MRampa;-MRampa];
B_=[B_,Rampa,Rampa];

%% RESOLVENDO
options = optimset('Display', 'iter', 'Algorithm', 'interior-point'); % Defina as op��es
[x,fval,exitflag,output,lambda] = linprog(f, A, B_, Aeq, Beq, lb,ub);

clear LGm DefI VAngLI LGM DefS VAngLS Big2 Big1 LimF i C M MatIgualdade Lines Lines2 Mflu_Original SomaPerdas
clear Desig_Fluxo A1 Aeq1 Colums Colums2 Des_DEFS Des_GERS Desig_Ang fluxos MatFluxo Linha k MRampa Rampa MATRIZ1


%% An�lise de Resultados

Custo_Hora=[];   % C�lculo do Custo por Hora
Emissao_Hora=[]; % C�lculo da Emiss�o de CO2 por Hora
Carga=[];        % Carga por Hora
Ger=[];          % Gera��o de cada unidade por Hora
Deficit=[];      % D�ficit de cada barra por Hora
C_Ger=[];        % C�lculo do Custo de Cada Gerador por Hora
E_Ger=[];        % C�lculo da Emiss�o de Cada Gerador por Hora

for i = 1 : length(DEMANDA(:,1))-1

    %                   Custo calculado na hora i
Custo_Hora=[Custo_Hora fcost(1+(2*NG+NG-1)*(i-1):2*NG+NG-1+(2*NG+NG-1)*(i-1))*x(1+(2*NG+NG-1)*(i-1):2*NG+NG-1+(2*NG+NG-1)*(i-1))];

    %                   Emiss�o calculada na hora i
Emissao_Hora=[Emissao_Hora femission(1+(2*NG+NG-1)*(i-1):2*NG+NG-1+(2*NG+NG-1)*(i-1))*x(1+(2*NG+NG-1)*(i-1):2*NG+NG-1+(2*NG+NG-1)*(i-1))];

    %                   Demanda na hora i
Carga=[Carga ;DEMANDA(i,2:end)];

    %                   Gera��o total na hora i
Ger=[Ger; x(1+(2*NG+NG-1)*(i-1):NG+(2*NG+NG-1)*(i-1))'];

    %                   Custo de Gera��o de cada gerador na hora i
C_Ger=[C_Ger; DGER(:,2)'.*x(1+(2*NG+NG-1)*(i-1):NG+(2*NG+NG-1)*(i-1))'];

    %                   Emiss�o de cada gerador na hora i
E_Ger=[E_Ger; DGER(:,6)'.*x(1+(2*NG+NG-1)*(i-1):NG+(2*NG+NG-1)*(i-1))'];


    %                   D�ficit de cada barra na hora i
Deficit=[Deficit; DGER(:,8)'.*x(1+NG+(2*NG+NG-1)*(i-1):NG+NG+(2*NG+NG-1)*(i-1))'];
end



toc
clc

%--------------------------------------
figure(1)
bar(Custo_Hora)
title('Custo Total de Gera��o x Hora')
ylabel('Custo ($)')
xlabel('Hora')

axis([0 4 0 1.1*max(Custo_Hora)])


%--------------------------------------
figure(2)
bar(Emissao_Hora)
title('Emiss�o Total de CO_2 x Hora')
ylabel('Emiss�o (m�)')
xlabel('Hora')

axis([0 4 0 1.1*max(Emissao_Hora)])


%--------------------------------------
figure(3)
bar(Carga)
title('Demanda do Sistema x Hora')
ylabel('Demanda (MW)')
xlabel('Hora')
legend('Barra 1','Barra 2','Barra 3')
axis([0 4 0 1.1*max(max(Carga))])

%--------------------------------------
figure(4)
bar(Ger)
title('Gera��o de Cada Gerador x Hora')
ylabel('Gera��o (MW)')
xlabel('Hora')
legend('Gerador 1','Gerador 2','Gerador 3')
axis([0 4 0 1.1*max(max(Ger))])
%--------------------------------------
figure(5)
bar(Deficit)
title('Custo de Deficit do Sistema x Hora')
ylabel('Custo ($)')
xlabel('Hora')
legend('Barra 1','Barra 2','Barra 3')

axis([0 4 0 1.1*max(max(Deficit))])

%--------------------------------------
figure(6)
bar(C_Ger)
title('Custo dos Geradores x Hora')
ylabel('Custo ($)')
xlabel('Hora')
legend('Gerador 1','Gerador 2','Gerador 3')

axis([0 4 0 1.1*max(max(C_Ger))])

%--------------------------------------
figure(7)
bar(E_Ger)
title('Emiss�o dos Geradores x Hora')
ylabel('Emiss�o m�)')
xlabel('Hora')
legend('Gerador 1','Gerador 2','Gerador 3')

axis([0 4 0 1.1*max(max(E_Ger))])
