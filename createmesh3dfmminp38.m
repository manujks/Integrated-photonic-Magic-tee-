function [nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminp38(ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
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

yf = 2*sptw+max([gl+w1+w3,gr+2*w2]);
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
dd=repmat(ypos(1:k+1),length(y),1)+sptw;
du=repmat(dd((1:size(dd,1)),1),1,length(x))+wdl-recta(xx,0,0,1,(wdl-wdr)/2);
ud=du+recta(xx,0,gl,1,gr);
uu=ud+wul+recta(xx,0,0,1,wur-wul);


for ii=1:ix
    dd=spline(xalpha(ii,:),dd,x);
    du=spline(xalpha(ii,:),du,x);
    ud=spline(xalpha(ii,:),ud,x);
    uu=spline(xalpha(ii,:),uu,x);
end




%dd=repmat(ypos(1:k),length(y),1)+sptw;


%dd=repmat(ypos(1:k),length(y),1)+sptw;
%du=repmat(dd((1:size(dd,1)),1),1,length(x))+wdl-recta(xx,0,0,1,(wdl-wdr)/2);
%ud=du+recta(xx,0,gl,1,gr);
%uu=ud+wul+recta(xx,0,0,1,wur-wul);

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

yf1 = 2*sptw+max([gl+w11+w31,gr+2*w21]);
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

dd1=repmat(ypos1(1:k+1),length(y),1)+sptw-(wdl1-wdr1)/2;
du1=repmat(dd1((1:size(dd1,1)),1),1,length(x))+wdl1-recta(xx,0,0,1,(wdl1-wdr1)/2);
ud1=du1+recta(xx,0,gl,1,gr);
uu1=ud1+wul1+recta(xx,0,0,1,wur1-wul1);
%dd1=repmat(ypos(1:k),length(y),1)+sptw*0.8;
%du1=dd1+wdl1-rectan(x,e1,e2,0,0,1,wdl-wdr);
%ud1=du1+rectagn(x,e1,e2,0,gl,1,gr);
%uu1=ud1+wul1+rectan2(x,e1,e2,0,0,1,wur-wul,w3,wmin);

for ii=1:ix
    dd1=spline(xalpha(ii,:),dd1,x);
    du1=spline(xalpha(ii,:),du1,x);
    ud1=spline(xalpha(ii,:),ud1,x);
    uu1=spline(xalpha(ii,:),uu1,x);
end


%dd1=repmat(ypos1(1:k),length(y),1)+sptw;


%dd1=repmat(ypos1(1:k),length(y),1)+sptw;
%du1=repmat(dd1((1:size(dd1,1)),1),1,length(x))+wdl1-recta(xx,0,0,1,(wdl1-wdr1)/2);
%ud1=du1+recta(xx,0,gl,1,gr);
%uu1=ud1+wul1+recta(xx,0,0,1,wur1-wul1);

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
 figure;
 imagesc(x,x,N1);
 figure;
 imagesc(x,x,N2);
 axis image;
 set(gca,'YDir','normal');
 title(['size(n) = ', num2str(size(N2))])
 colormap(jet);
 colorbar;
end
l=22;
x=l*x;
x=x-l-dz-25;
bv =[[x'; fliplr(x)'], [dd' ; fliplr(du)']]*1e-6;
tv=[[x'; fliplr(x)'], [ud' ; fliplr(uu)']]*1e-6;
bv=bv(:,1:2);
tv=tv(:,1:2);
figure;
plot(bv(:,1),bv(:,2))
%convert units to meters
bendy=0;bendx=0;lam=1.55;wout=0;
w1=w1*1e-6;
w2=w2*1e-6;
w3=w3*1e-6;
w4=w4*1e-6;
gl=gl*1e-6;
gr=gr*1e-6;
wout=wout*1e-6;
lam=lam*1000;
bendy=bendy*1e-6;
bendx=bendx*1e-6;
bendxin=3e-6;


bv1 =[[x'; fliplr(x)'], [dd1' ; fliplr(du1)']]*1e-6;
tv1=[[x'; fliplr(x)'], [ud1' ; fliplr(uu1)']]*1e-6;
bv1=bv1(:,1:2);
tv1=tv1(:,1:2);figure;
plot(bv1(:,1),bv1(:,2))
%convert units to meters
bendy=0;bendx=0;lam=1.55;wout=0;
w11=w11*1e-6;
w21=w21*1e-6;
w31=w31*1e-6;
w41=w41*1e-6;
gl=gl*1e-6;
gr=gr*1e-6;
wout=wout*1e-6;
lam=lam*1000;
bendy=bendy*1e-6;
bendx=bendx*1e-6;
bendxin=3e-6;
cof=0;





end

function y=recta(x,x0,y0,x1,y1)
y=(y1-y0)/(x1-x0)*(x-x0)+y0;
end

function y = fd(x,dy)
a=2;
y=1./(exp((a*x/(1.5*dy))*2*log(3))+1);  %Fermi-dirac distribution to smooth index jump/reduce discretization errors
end