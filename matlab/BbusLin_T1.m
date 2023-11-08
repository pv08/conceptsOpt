%% FLUXO DC
function [Gbus,Bbus,BbusLin]=BbusLin_T1(C,Linha)
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

end
