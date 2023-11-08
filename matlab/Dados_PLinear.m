%% CONCEITOS CL�SSICOS DE OTIMIZA��O
%% Dados Programa��o Linear - TRABALHO 1 
global DGER DLIN DEMANDA

% DADOS DE GERA��O
     % BARRA  CUSTO($/MWh)  MAX(MW) MIN(MW)  RAMPA(MW)    C02(m3/MWh)
DGER=[  1        10             30      05       10            90;
        2        30             40      15       05            10;
        3       100             40      00       03            70];

%###########################################
% fiquem a vontade para usar outros valores:    
h=0.25; % fator de convers�o emiss�o - $/MW 
delta=1; % fator de pondera��o emiss�o     
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
 DEMANDA= [   1     0      40      30;
              2     0      43      25;
              3     0      25      25];