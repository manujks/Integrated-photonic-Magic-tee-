function [nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminp21(dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx x xip dimension 
%Ly y xip dimension 
%l coupling length
%g gap coupling distance
%w1 waveguide width up
%w2 waveguide width down
%dx,dy discretization
%nBody refractive index waveguide
%nBkgrnd refractive index enviroment
%m upper margin distance
% example Lx=100; Ly=50; l=30; g=0.5; w = 0.3; m=4; dx= 0.01; dy = 0.01; 
%size in pixels of the filter hsize=20
l=dz*k;
%max([gl+w1+w3,gr+2*w2])
yf = 2*sptw+gl+w11+w31;
yf=yf-rem(yf,dx);

wul=w3;       wur=w2;
wdl=w1;       wdr=w4;
%yzL=yzL+max(abs(yzL))+5*max([wdr,wur]);
%yu=(yzL(end)+wur-yzL(1)-wul)/(xL(end)-xL(1))*(xL(1)-xL(1))+yzL(1)+wul+gl;

%yf = yu+5*max([wdr,wur]);
%yf=yf-rem(yf,dx);


x=dz:dz:l+dz;
y=dx/2:dx:yf-dx/2;
[xx,yy]=meshgrid(x,y);


%dd=du-wdl+recta(xx,0,0,1,wdl-wdr);
%du=spline(xL,yzL,xx);

%dd=repmat(ypos(1:k),length(y),1)+sptw;
dd=repmat(ypos(1:k+1),length(y),1)+sptw;
du=repmat(dd((1:size(dd,1)),1),1,length(x))+wdl-recta(xx,0,0,1,(wdl-wdr)/2);
ud=du+recta(xx,0,gl,1,gr);
uu=ud+wul+recta(xx,0,0,1,wur-wul);

botu =1-fd(-yy+ud,dy);
topu =fd(-yy+uu,dy);
top=botu+topu;

topd=fd(-yy+du,dy);
botd=1-fd(-yy+dd,dy);
down=botd+topd;

tot=top.*down;
N1 = (1-tot)*(nCore-nAbove)+nAbove;




%Lx x xip dimension 
%Ly y xip dimension 
%l coupling length
%g gap coupling distance
%w1 waveguide width up
%w2 waveguide width down
%dx,dy discretization
%nBody refractive index waveguide
%nBkgrnd refractive index enviroment
%m upper margin distance
% example Lx=100; Ly=50; l=30; g=0.5; w = 0.3; m=4; dx= 0.01; dy = 0.01; 
%size in pixels of the filter hsize=20
l=dz*k;
%max([gl+w11+w31,gr+2*w21]);
yf1 = 2*sptw+gl+w11+w31;
yf1=yf1-rem(yf1,dx);

wul1=w31;       wur1=w21;
wdl1=w11;       wdr1=w41;
%yzL=yzL+max(abs(yzL))+5*max([wdr,wur]);
%yu=(yzL(end)+wur-yzL(1)-wul)/(xL(end)-xL(1))*(xL(1)-xL(1))+yzL(1)+wul+gl;

%yf = yu+5*max([wdr,wur]);
%yf=yf-rem(yf,dx);


x=dz:dz:l+dz;
y=dx/2:dx:yf-dx/2;
[xx,yy]=meshgrid(x,y);


%dd=du-wdl+recta(xx,0,0,1,wdl-wdr);
%du=spline(xL,yzL,xx);

%dd1=repmat(ypos1(1:k),length(y),1)+sptw;
dd1=repmat(ypos1(1:k+1),length(y),1)+sptw-(wdl1-wdr1)/2;
du1=repmat(dd1((1:size(dd1,1)),1),1,length(x))+wdl1-recta(xx,0,0,1,(wdl1-wdr1)/2);
ud1=du1+recta(xx,0,gl,1,gr);
uu1=ud1+wul1+recta(xx,0,0,1,wur1-wul1);

botu1 =1-fd(-yy+ud1,dy);
topu1 =fd(-yy+uu1,dy);
top1=botu1+topu1;

topd1=fd(-yy+du1,dy);
botd1=1-fd(-yy+dd1,dy);
down1=botd1+topd1;

tot1=top1.*down1;
N2 = (1-tot1)*(nCore-nAbove)+nAbove;






[nn,pp]=size(xx(:,1:k+1));

nlyrs=zeros(4,nn,pp);
nlyrs(1,:,:)=nBelow;
nlyrs(2,:,:)=N1;
nlyrs(3,:,:)=N2;
nlyrs(4,:,:)=nAbove;

dlyrsy=[yCladBot,hWg1,hWg2,yCladTop]';
dlyrsy=repmat(dlyrsy,[1,pp]);
dlyrsx=ones(nn,pp)*dx;

if strcmp(ifplot, 'yes')
 figure(15);
 imagesc(x,y,N1);
 figure;
 imagesc(x,y,N2);
 
end

end

function y=recta(x,x0,y0,x1,y1)
y=(y1-y0)/(x1-x0)*(x-x0)+y0;
end

function y = fd(x,dy)
a=2;
y=1./(exp((a*x/(1.5*dy))*2*log(3))+1);  %Fermi-dirac distribution to smooth index jump/reduce discretization errors
end
