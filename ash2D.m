function [newgridX,newgridY,fash] = ash2D(data,m,type)

[~,~,data]= format_data(data);

[na,~,~] = histcounts2(data(1,:),data(2,:),'Normalization', 'probability','BinMethod','fd');
bin1 = size(na,1);
bin2 = size(na,2);

h1=(max(data(1,:))-min(data(1,:)))/bin1;
h2=(max(data(2,:))-min(data(2,:)))/bin2;

t01=linspace(0,(h1-(h1/m)),m)';
t02=linspace(0,(h2-(h2/m)),m)';
for k=1:m
    
    gridx=(min(data(1,:))+t01(k)):h1:(max(data(1,:))+t01(k));
    gridy=(min(data(2,:))+t02(k)):h2:(max(data(2,:))+t02(k));
    
    for i=1:bin1
        for j=1:bin2
            a=find(data(1,:)>=gridx(i) & data(1,:)<=(gridx(i+1)));
            b=find(data(2,:)>=gridy(j) & data(2,:)<=(gridy(j+1)));
            
            
            f(j,i,k)=length(intersect(a,b));
            y(k,j)=(gridy(j)+gridy(j+1))/2;
        end
        x(k,i)=(gridx(i)+gridx(i+1))/2;
        
    end
end

newgridX=x(:)';
newgridY=y(:)';

for k=1:m
    for i=1:length(newgridX)
        fsh(:,i,k)=interp2(x(k,:),y(k,:),f(:,:,k),newgridX(i),newgridY,type,0);
    end
end

fash=mean(fsh,3);

[V]= volume_pts(newgridX,newgridY,fash,length(fash));
fash=fash/V;

end

function [nd,n,data]= format_data(data)

[nd,n] = size(data);

if nd>n
    data = data';
    [nd,n] = size(data);
end

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
