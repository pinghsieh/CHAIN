function [xout, yout, xout2, yout2] = fill_interior(xmax,yinput,xlimit,ylimit,Delta)

xout = [];
yout = [];
xout2 = [];
yout2 = [];

for id=0:length(yinput)-1
   vec = double(id)*Delta*ones(1,round((yinput(id+1))/Delta) + 1);
   xout = cat(2,xout,vec);
   yout = cat(2,yout,0:int16(Delta):int16(yinput(id+1)));
   xout2 = cat(2,xout2,double(id)*Delta*ones(1,round((ylimit-yinput(id+1))/Delta)));
   yout2 = cat(2,yout2,(yinput(id+1)+Delta):Delta:ylimit);
end
for j=xmax+Delta:Delta:xlimit
   xout2 = cat(2,xout2,double(j)*ones(1,round(ylimit/Delta)+1));
   yout2 = cat(2,yout2,0:Delta:ylimit);
end
end