%% IMPLEMENTA��O FLUXO DC COM E SEM PERDAS
% clc 
% clear all
global DGER DLIN DEMANDA Sbase SW


% Pot�ncia m�xima do sistema 
Sbase = 100;

% C n�mero de barras

% DADOS DE BARRA
% 1� coluna: n�mero da barra
% 2� coluna: tipo de barra: 0 (barra SW), 1 (barra PV) ou 2 (barra PQ)
% 3� coluna: capacidade m�xima de gera��o de pot�ncia ativa (em pu)
% 4� coluna: capacidade m�xima de gera��o de pot�ncia reativa (em pu)
% 5� coluna: pot�ncia ativa da carga (em pu)
% 6� coluna: pot�ncia reativa da carga (em pu)
% 7� coluna: pot�ncia reativa de banco de capacitores(+) ou indutores(-) (em pu)
% 8� coluna: tens�o da barra (em pu)- valor 1 caso n�o seja definido
% 9� coluna: angulo da tens�o (em graus)- valor 0 caso n�o seja definido

C=length(DGER(:,1));
B=DGER(:,7); 
PG=DGER(:,3)/Sbase;   %% Inserir Gera��o
PD=DEMANDA(TimeD,2:end)'/Sbase;

SW=find(B==0);
   
% DADOS DAS LINHAS
% 1� e 2� coluna: n�mero da barra
% 3� coluna: resist�ncia da linha (em pu)
% 4� coluna: reat�ncia da linha (em pu)
% 5� coluna: sucept�ncia shunt da linha (em pu)
% 6� coluna: Defasagem do transformador ( 1:e^jPhs ) em graus
   

%      D P    R       X        Sh     Phs
Linha(:,1)= DLIN(:,1);  
Linha(:,2)= DLIN(:,2); 

Linha(:,4)=imag(1./((DLIN(:,4)-1j*DLIN(:,3))/Sbase));
Linha(:,3)=real(1./((DLIN(:,4)-1j*DLIN(:,3))/Sbase));

Linha(:,5)=zeros(length(Linha(:,1)),1);   
Linha(:,6)=zeros(length(Linha(:,1)),1);  


% bs=  pot�ncia aparente base (MVA)
bs= Sbase;

% [Gbus,Bbus,BbusLin,Ang,Perdas,Fl,PinjTotal]=FlowDC_T1(C,PG,PD,bs,Linha);
[Gbus,Bbus,BbusLin]=BbusLin_T1(C,Linha);
