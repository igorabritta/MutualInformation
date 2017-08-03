function [X,pdf] = fastKDE(data,nPoint,f)
% l=1;
doNorm = 0;

tic
[nd,n,data]= format_data(data);

Div.L=10;
Div.S=2;
Div.D=5;

if nd ==1
    [hist_pdf,Xedges,~] = histcounts(data,'Normalization', 'pdf','BinMethod','fd');
    xh = (Xedges(1:end-1) + Xedges(2:end))/2;
    optN = length(hist_pdf);
    proby=interp1(xh(find(hist_pdf~=0)),hist_pdf(find(hist_pdf~=0)),data,'linear','extrap');           %cálculo das probabilidades dos eventos de kernel
    lambda=exp((length(data)^-1)*sum(log(proby(proby~=0)))); %calculo do lambda otimo pela teoria do paper
    h =((4/(nd+2))^(1/(nd+4)))*std(data)*n^(-1/(nd+4)); % cálculo do h("rule-of-thumb") ótimo pelo AMISE
    X = linspace(min(data),max(data),nPoint); %range X do kernel
else
    [hist_pdf,Xedges,Yedges,~] = histcounts2(data(1,:),data(2,:),'Normalization', 'pdf','BinMethod','fd');
    xh = (Xedges(1:end-1) + Xedges(2:end))/2;
    yh = (Yedges(1:end-1) + Yedges(2:end))/2;
    optN = length(hist_pdf);
    proby=interp2(xh,yh,hist_pdf',data(1,:),data(2,:),'linear',min(min(hist_pdf(hist_pdf>0))));           %cálculo das probabilidades dos eventos de kernel
    lambda=exp((length(data)^-1)*sum(log(proby(proby~=0)))); %calculo do lambda otimo pela teoria do paper
    for i=1:nd
        h(i) =((4/(nd+2))^(1/(nd+4)))*std(data(i,:))*n^(-1/(nd+4)); % cálculo do h("rule-of-thumb") ótimo pelo AMISE
        X(i,:) = linspace(min(data(i,:)),max(data(i,:)),nPoint); %range X do kernel
    end
end

h=h*f;
[hi] = hihj(h,lambda,proby,n);         % cálculo da banda variável
Hi =hi.^2;


%==========================================================================
% Fazendo os Cálculos do KERNEL ND de banda Fixa e Variável
%==========================================================================
if nd ==1
    %% Dimensões = 1
    if length(data)>=69998
        pdf=[];
        for j=1:Div.L
            pdf=[pdf ((1/n)*sum((repmat((Hi.^(-1/2)),nPoint/Div.L,1).*Kn(repmat((Hi.^(-1/2)),nPoint/Div.L,1).*((repmat(X((1+(nPoint/Div.L)*(j-1)):(nPoint/Div.L)*(j)),length(data),1)')-repmat(data,nPoint/Div.L,1)))),2))'];
        end
    else
        pdf=[];
        for j=1:Div.S
            pdf=[pdf ((1/n)*sum((repmat((Hi.^(-1/2)),nPoint/Div.S,1).*Kn(repmat((Hi.^(-1/2)),nPoint/Div.S,1).*((repmat(X((1+(nPoint/Div.S)*(j-1)):(nPoint/Div.S)*(j)),length(data),1)')-repmat(data,nPoint/Div.S,1)))),2))'];
        end
    end
    
    if doNorm == 1
        [A]= area2d(X,pdf);
        pdf = pdf/A;
    end
else
    %% Dimensões > 1
    if length(data)<=20000
        MH1= repmat((Hi(1,:).^(-1/2)),nPoint,1);
        MH2= repmat((Hi(2,:).^(-1/2)),nPoint,1);
        MD1 = repmat(data(1,:),nPoint,1);
        MD2 = repmat(data(2,:),nPoint,1);
        MX1=repmat(X(1,:),length(data(1,:)),1)';
        MX2=repmat(X(2,:),length(data(2,:)),1)';
        
        pdf=[(1/n)*(MH1.*Kn(MH1.*(MX1-MD1)))*(MH2.*Kn(MH2.*(MX2-MD2)))']';
        if doNorm == 1
            [V]= volume_pts(X(1,:),X(2,:),pdf,nPoint);
            pdf = pdf/V;
        end
    else
        pdf=[];
        for jj=1:Div.D
            pdfy=[];
            for kk=1:Div.D
                MH1= repmat((Hi(1,:).^(-1/2)),nPoint/Div.D,1);
                MH2= repmat((Hi(2,:).^(-1/2)),nPoint/Div.D,1);
                MD1 = repmat(data(1,:),nPoint/Div.D,1);
                MD2 = repmat(data(2,:),nPoint/Div.D,1);
                MX1=repmat(X(1,((1+(nPoint/Div.D)*(jj-1)):(nPoint/Div.D)*(jj))),length(data(1,:)),1)';
                MX2=repmat(X(2,((1+(nPoint/Div.D)*(kk-1)):(nPoint/Div.D)*(kk))),length(data(2,:)),1)';
                
                pdf1=[(1/n)*(MH1.*Kn(MH1.*(MX1-MD1)))*(MH2.*Kn(MH2.*(MX2-MD2)))']';
                pdfy=[pdfy; pdf1];
            end
            pdf=[pdf pdfy];
        end
        if doNorm == 1
            [V]= volume_pts(X(1,:),X(2,:),pdf,nPoint);
            pdf = pdf/V;
        end
        
    end
    
    
    
end
toc
end


function [nd,n,data]= format_data(data)

[nd,n] = size(data);

if nd>n
    data = data';
    [nd,n] = size(data);
end

end


function [hi] = hihj(h,lambda,fpi,n)
for i=1:n
    hi(:,i)=abs((h).*(sqrt(lambda./fpi(i))));
end
end

function [K]=Kn(u)
K=((2*pi)^(-1/2))*exp((-(u.^2)/2)); %Gaussian
end

function [x,y] = data_normalized(x,bin)
[yh,xh]= hist(x,bin);
ah= area2d(xh,yh);
x = xh;
y=yh/ah;
end

function [ area ] = area2d(x,y)
tbin=min(diff(x));
area=sum(abs(y))*tbin;
end

function [V]= volume_pts(x,y,z,pts);

xgrid = min(x):range(x)/pts:max(x);
ygrid = min(y):range(y)/pts:max(y);

h = waitbar(0,'Remaking GRID points');
for i=1:length(xgrid)
    zgrid(i,:)= interp2(x,y,z,xgrid(i),ygrid,'linear');
    waitbar(i/length(xgrid))
end
close(h)
dx=min(diff(xgrid));
dy=min(diff(ygrid));
V=sum(sum(dx*dy*zgrid));
end
