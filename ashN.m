function [x,y] = ashN(x,M)
[xh,yh]= ash(x,M,'linear');
ah= area2d(xh,yh);
x = xh;
y=yh/ah;
end

function [ area ] = area2d( x,y )
%tbin=abs((x(2)-x(1)));
tbin=min(diff(x));
%area=trapz(x,y);
area=sum(abs(y))*tbin;

end

