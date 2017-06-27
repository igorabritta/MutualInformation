function [I,In]=mi(A,B,method)
% ALGORITMO TESTADO NO MATLAB R2015b
% =========================================================================
% FUNÇÃO REALIZAR O CALCULO DE INFORMAÇÃO MÚTUA ENTRE DUAS DISTRIBUIÇÕES
% =========================================================================
% ENTRADA:      A = Distribuição 1
%               B = Distribuição 2
%               method = Método de Estimação de Densidades
%                       ->'HIST' = Estimação utilizando HISTOGRAMA 
%                       ->'ASH' = Estimação utilizando AVERAGE SHIFTED HISTOGRAM 
%                       ->'KDE' = Estimação utilizando KERNEL DENSITY ESTIMATION 
% =========================================================================
%   Class support for inputs A, B: 
%      float: double, float, array of columns
% =========================================================================
% SAÍDA:        I = INFORMAÇÃO MÚTUA
%               In = INFORMAÇÃO MÚTUA NORMALIZADA
% =========================================================================
% FUNÇÕES NECESSÁRIAS: fastKDE.m;ashN.m;ash2D.m;
% =========================================================================

%==========================================================================
% Escolha dos métodos de construção das Probabildades utilizadas no cálculo
% da Informação mútua:
%==========================================================================
%% Method 1): Histograma
if strcmp(method,'HIST')
    [na,Aedges,~] = histcounts(A,'Normalization', 'probability','BinMethod','fd');
    xA = (Aedges(1:end-1) + Aedges(2:end))/2;
    na = interp1(xA(na>0),na(na>0),xA,'nearest','extrap');
    [nb,Bedges,~] = histcounts(B,Aedges,'Normalization', 'probability');
    xB = (Bedges(1:end-1) + Bedges(2:end))/2;
    nb = interp1(xB(nb>0),nb(nb>0),xB,'nearest','extrap');
    [n2,~,~,~] = histcounts2(A,B,Aedges,Bedges,'Normalization', 'probability');
    % n2 = smoothn(n2,3);
end
%% Method 2): Kernel Density Etimation
if strcmp(method,'KDE')
    [~,na] = fastKDE(A,200,1);
    [~,nb] = fastKDE(B,200,1);
    [~,n2] = fastKDE([A B],200,1);
    n2 = n2';
    na = na/sum(na);
    nb = nb/sum(nb);
    n2 = n2/sum(sum(n2));
end
%% Method 3): Average Shifted Histogram
if strcmp(method,'ASH')
    [xa,na] = ashN(A,5);
    [xb,nb] = ashN(B,5);    
    [newgridX,newgridY,fash] = ash2D([A B],5,'linear');
    for i=1:length(xa)
        n2(i,:) = interp2(newgridX,newgridY,fash,repmat(xa(i),1,length(xb)),xb,'linear',min(min(fash(fash>0))));
    end
    na = na/sum(na);
    nb = nb/sum(nb);
    n2 = n2/sum(sum(n2));
end

%==========================================================================
% Cálculo da Informação mútua:
%==========================================================================
%% Entropia:
ha=-sum(na(na>0).*log2(na(na>0)));
hb=-sum(nb(nb>0).*log2(nb(nb>0)));
%% Probabilidade feita com as marginais:
nanb=(na'*nb);
nanb(nanb==0)=1e-12;
n2(n2==0)=1e-12;
%% Informação Mútua:
I = sum(sum(n2.*(log2(n2./nanb))));
%% Normalização(Literatura): 
In = (I/sqrt(ha*hb));
% In =(2*I)/(ha+hb);
% In = (I/max([ha hb]))