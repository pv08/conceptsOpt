%% IMPLEMENTAÇÃO FLUXO DC COM E SEM PERDAS
% clc 
% clear all
global DGER DLIN DEMANDA Sbase SW


% Potência máxima do sistema 
Sbase = 100;

% C número de barras

% DADOS DE BARRA
% 1ª coluna: número da barra
% 2ª coluna: tipo de barra: 0 (barra SW), 1 (barra PV) ou 2 (barra PQ)
% 3ª coluna: capacidade máxima de geração de potência ativa (em pu)
% 4ª coluna: capacidade máxima de geração de potência reativa (em pu)
% 5ª coluna: potência ativa da carga (em pu)
% 6ª coluna: potência reativa da carga (em pu)
% 7ª coluna: potência reativa de banco de capacitores(+) ou indutores(-) (em pu)
% 8ª coluna: tensão da barra (em pu)- valor 1 caso não seja definido
% 9ª coluna: angulo da tensão (em graus)- valor 0 caso não seja definido

C=length(DGER(:,1));
B=DGER(:,7); 
PG=DGER(:,3)/Sbase;   %% Inserir Geração
PD=DEMANDA(TimeD,2:end)'/Sbase;

SW=find(B==0);
   
% DADOS DAS LINHAS
% 1ª e 2ª coluna: número da barra
% 3ª coluna: resistência da linha (em pu)
% 4ª coluna: reatância da linha (em pu)
% 5ª coluna: suceptância shunt da linha (em pu)
% 6ª coluna: Defasagem do transformador ( 1:e^jPhs ) em graus
   

%      D P    R       X        Sh     Phs
Linha(:,1)= DLIN(:,1);  
Linha(:,2)= DLIN(:,2); 

Linha(:,4)=imag(1./((DLIN(:,4)-1j*DLIN(:,3))/Sbase));
Linha(:,3)=real(1./((DLIN(:,4)-1j*DLIN(:,3))/Sbase));

Linha(:,5)=zeros(length(Linha(:,1)),1);   
Linha(:,6)=zeros(length(Linha(:,1)),1);  


% bs=  potência aparente base (MVA)
bs= Sbase;

% [Gbus,Bbus,BbusLin,Ang,Perdas,Fl,PinjTotal]=FlowDC_T1(C,PG,PD,bs,Linha);
[Gbus,Bbus,BbusLin]=BbusLin_T1(C,Linha);
