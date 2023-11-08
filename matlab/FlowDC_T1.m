%% FLUXO DC
function [Gbus,Bbus,BbusLin,A,Perdas,Fl,PinjTotal]=FlowDC_T1(C,PG,PD,bs,Linha)
%% *Construção de YBUS
global DLIN SW

L=length(DLIN(:,1));
BARRADE    = DLIN(:,1);           % barra de
BARRAPARA  = DLIN(:,2);           % barra para

DE = BARRADE;
PARA = BARRAPARA;

nlin = size(DLIN(:,1));     % Número de Circuitos do Sistema
R = Linha(:,3);
X = Linha(:,4);
YKM = 1./(R + 1j*X);        % ykm de cada circuito
% -----------------------------------------------------------


Ybus = zeros(C);
BbusLin = zeros(C);
for s = 1:nlin
    
    k = DE(s);
    m = PARA(s);
    % Ybus Tradicional
    Ybus(k,k) = Ybus(k,k) + YKM(s,1);
    Ybus(k,m) = Ybus(k,m) - YKM(s,1);
    Ybus(m,k) = Ybus(m,k) - YKM(s,1);
    Ybus(m,m) = Ybus(m,m) + YKM(s,1);
    
    % Bbus' do fluxo DC
    BbusLin(k,k) = BbusLin(k,k) + 1./(1j*X(s,1));
    BbusLin(k,m) = BbusLin(k,m) - 1./(1j*X(s,1));
    BbusLin(m,k) = BbusLin(m,k) - 1./(1j*X(s,1));
    BbusLin(m,m) = BbusLin(m,m) + 1./(1j*X(s,1));
    if k == SW
        BbusLin(k,k) = 1j*10^15;
    end
    
end

 % Ybus = G + j*B
 Gbus = real(Ybus);
 Bbus = imag(Ybus);
 
 BbusLin = imag(BbusLin);
%-------------------------------------------------------------------------- 
%% Injeções de Potência

Pinj=PG-PD;   %Potência injetada na barra 'i' como a diferença entre a geração e a demanda

Pinj(SW,1)=sum(-Pinj)+Pinj(SW,1);  % Manter P(SW) como a soma exata do negativo das demais potências injetadas

A=-BbusLin^-1*Pinj;  % Ângulos calculados pelo modelo linearizado (radianos)
A(SW,1)=0; % Ângulo da barra SW mantido na referência


%--------------------------------------------------------------------------
%% Inclusão das Perdas no Sistema

Perdas=zeros(C,1);
    for i=1:C
        incremento=0;
        for k=1:C
            incremento=incremento-1/2*(A(i,1)-A(k,1))^2*Gbus(i,k);  % Cálculo das perdas acumuladas ( barra 'i' para barra 'k')
        end
        Perdas(i,1)=incremento;   
    end
    
PinjTotal=Pinj-Perdas;   % Potências ativas injetadas, incluindo perdas

PinjTotal(SW,1)=sum(-PinjTotal)+PinjTotal(SW,1);  % Manter P(SW) como a soma exata do negativo das demais potências
 
AN=-BbusLin^-1*PinjTotal;    % Ângulo calculado incluindo as perdas
AN(SW,1)=0; % Ângulo da barra SW mantido na referência


PinjTotal=Pinj;    % Potências ativas injetadas, perdas incluídas como injeção da barra 'SW'
PinjTotal(SW,1)=sum(-Pinj)+Pinj(SW,1)+sum(Perdas); % Total de injeção da barra SW

disp('Ângulos Resultantes: Com Perdas (°)')
disp(AN*180/pi)
disp('Potências Ativas Injetadas: Com Perdas (MW)')
disp(PinjTotal*bs)

%% Fluxo de Potência com Perdas
Flkm=zeros(L,3);
Flmk=zeros(L,3);
for s=1:L  
    
    k = DE(s);
    m = PARA(s); 
    if k~=m
        % Fluxo K-M
        Flkm(s,1)=k;
        Flkm(s,2)=m;
        Flkm(s,3)= imag(1./(1j*X(s,1)))*(AN(DLIN(s,1))-AN(DLIN(s,2)));
        % Fluxo M-K        
        Flmk(s,1)=m;
        Flmk(s,2)=k;
        Flmk(s,3)= imag(1./(1j*X(s,1)))*(AN(DLIN(s,2))-AN(DLIN(s,1)));
    end
end
Fl=[Flkm;Flmk];
end
