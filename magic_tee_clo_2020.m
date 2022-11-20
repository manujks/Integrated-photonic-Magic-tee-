%%
close all
clear 
 
dx=0.01;
dy=0.01;
dz=0.0025;
%dz=0.01;
w1=0.38;%WD1
w2=0;%WUR
w3=0;%wul
%w4=0.84;%WDRfirst
w4=0.38;%WDR
%w11=0.38;
w11=0.1;
w21=0;
w31=0;
w41=0.38;
%w41=0.05
gl=0;
gr=0;
sptw=2*w1;
yCladBot=0.68;
hWg1=0.15;
hWg2=0.07;
yCladTop=0.64;%0.64
nBelow = 1.444;  
nCore  =  3.48;  

nAbove =  1.444;   
ifplot='yes';
k=(1/dz);


%yposf=spline(xL,yzL,xx(1,1:k));
x=dz:dz:1+dz;
kmin=1e-8;
%ypos=[0 1 -1 0];
%k=length(ypos);
k0 = 2*pi/1.55;
%ypos=zeros(1,k);
ypos=recta(x,0,0,1,-(-w1+w4)/2);%%%%%%%%%%%errrorr
%[K,B,nslice] = fmm3dTEmodificat(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot, k0, ...
 %              struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
  %                    'BC',[0 0 0 0],'enginever','m2wcylR2', ...
   %                   'operver','m2dpmloperR2'));  
 ypos1=recta(x,0,0,1,-(-w11+w41)/2);  
Fold=0;
nz=1;
%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm21(dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm21(dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminp21(dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx
%Lx x
%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminp21(L,ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%
%Lx x xip dimension 
C=cat(1,nlyrs,nlyrs);figure;
[K,B,nslice] = fmm3d(nlyrs, dlyrsx, dlyrsy, dx, dy, k0, ...
               struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                      'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','eigmode','b','epsavg','no')); 
z = [dz:dz:1+dz];

zmult = [1:1:200];
figure;
S = fmm3dprop(K,B,z,zmult);
NMODES = size(S,1)/2;


S11 = abs(squeeze(S(NMODES+1,1,:).^2)); S22 = abs(squeeze(S(NMODES+2,2,:).^2)); S21 = abs(squeeze(S(NMODES+1,2,:).^2)); S31 = abs(squeeze(S(NMODES+1,3,:).^2));   % Forwards
S12 = abs(squeeze(S(NMODES+2,1,:).^2)); S23 = abs(squeeze(S(NMODES+3,2,:).^2)); S32 = abs(squeeze(S(NMODES+3,2,:).^2)); 
%S11b = abs(squeeze(S(1,NMODES+1,:).^2)); S21b = abs(squeeze(S(2,NMODES+1,:).^2));  % Backwards
S31 = abs(squeeze(S(NMODES+1,3,:).^2));S13 = abs(squeeze(S(NMODES+3,1,:).^2));
%figure; subplot(1,2,1); loglog(zmult, S11,'LineWidth',2.5); hold on; loglog(zmult, S21,'LineWidth',2.5);  loglog(zmult, S31,'LineWidth',2.5);
%xlabel('Length (\mum)'); ylabel('Power transmission coeff'); legend('S11','S21','S31') ; grid on;
%set(gca,'fontsize',16)

figure;loglog(zmult, S11,'LineWidth',3); hold on; loglog(zmult, S22,'LineWidth',2.5);  loglog(zmult, S32,'LineWidth',2.5);
xlabel('Length (\mum)'); ylabel('Power transmission coeff'); legend('S11','S22','S32') ; grid on;
set(gca,'fontsize',16)

%figure;loglog(zmult, S31,'LineWidth',3); hold on; loglog(zmult, S13,'LineWidth',2.5);  loglog(zmult, S32,'LineWidth',2.5);
%label('Length (\mum)'); ylabel('Power transmission coeff'); legend('S31','S13','S32') ; grid on;
%set(gca,'fontsize',16)



figure;
plot(dz:dz:1+dz,B/(2*pi/1.55),'LineWidth',2.5); xlabel('length [\mum]');ylabel('N_{eff}');legend('TE1','TE2','TM1');
set(gca,'fontsize',16)


%figure(5); plot(dz:dz:1-dz,abs(squeeze(K(1,2,:))),'LineWidth',2.5);hold on%;plot(z(2:end-1),squeeze(K(1,3,:)),'LineWidth',2);
%plot(xL,zL,'ok','LineWidth',2.5);
%plot(dz:dz:1-dz,abs(squeeze(K(2,1,:))),'ok','LineWidth',2);
%xlabel('Length (\mum)'); ylabel('Coupling coeff');
%legend('\kappa_{12}','\kappa_{21}')
%set(gca,'fontsize',16)

figure; plot(dz:dz:1,abs(squeeze(K(2,3,:))),'LineWidth',2.5);hold on%;plot(z(2:end-1),squeeze(K(1,3,:)),'LineWidth',2);
%plot(xL,zL,'ok','LineWidth',2.5);
plot(dz:dz:1,abs(squeeze(K(3,2,:))),'ok','LineWidth',2);
xlabel('Length (\mum)'); ylabel('Coupling coeff');
legend('\kappa_{23}','\kappa_{32}')
set(gca,'fontsize',16)
figure;
ste=200;K23=K(2,3,:);

K23=K23/max(squeeze(K23))*100;

%K23=K(2,3,:)./(B(2)-B(3));

KK=K23+ste;
figure;
plot(dz:dz:1,abs(squeeze(KK)),'ok','LineWidth',2);

alpha=KK/sum(KK);
xalpha=cumsum(alpha,3);
[ix,~]=size(xalpha);
xalpha(:,:,401)=0;
%xalpha(:,:,401)=0;
%xalpha(:,:,1)=0.01;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm35(xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k);
figure;
L=25;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm_138_35(L,ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx x xip dime
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminp38(ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%L
L=22
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmminput138_revised(L,ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx x xip dimension 
%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm35(xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx x xip dimension 
%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm_138_1(L,ste,xalpha,ix,dx,dy,dz,w1,w2,w3,w4,w11,w21,w31,w41,gl,gr,sptw,yCladBot,hWg1,hWg2,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,ypos1,k)
%Lx x xip dimension 







%%
Fold=0;
nz=1;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm2(dx,dy,dz,w1,w2,w3,w4,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  


[F,kl] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
               struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                      'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                      'operver','m2dpmloperR2','epsavg','no'));  
Fold=F; 
e=0.01*dx;%1e-7;
a = 0;

for nz=2
ypos(nz)=a;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
[F,fa] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  

ypos(nz)=a+e;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
[~,fae1] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
ypos(nz)=a-e;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
[~,fae2] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));                           
                         
err=abs(fa);
figure(8);
semilogy(nz,err,'-ok'); hold on; drawnow
figure(11);
plot(nz,ypos(nz),'-pk'); hold on; drawnow
        it=0;
        while err > kmin && it<3
            a = a - 2*e*fa/(fae1-fae2) ;
            ypos(nz)=a;
            [nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
            [F,fa] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
            err=abs(fa);
            if err > kmin
            ypos(nz)=a+e;
            [nlyrse,dlyrsxe, dlyrsye]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
            [~,fae1] = k12(Fold,nz,dx,dy,dz,nlyrse,dlyrsxe, dlyrsye,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
            ypos(nz)=a-e;
            [nlyrse,dlyrsxe, dlyrsye]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
            [~,fae2] = k12(Fold,nz,dx,dy,dz,nlyrse,dlyrsxe, dlyrsye,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
            
            ypos(nz)=a;
            figure(8);
            semilogy(nz,err,'-ok'); hold on; drawnow
            figure(11);
            plot(nz,ypos(nz),'-pk'); hold on; drawnow
            it=it+1;
            end
        end
        figure(8);
        semilogy(nz,err,'-or'); hold on; drawnow
        figure(11);
        plot(nz,ypos(nz),'-or'); hold on; drawnow
         
        
        
Nkt=5;
 xvec=linspace(-ypos(nz)+ypos(2),-ypos(nz)-ypos(2),Nkt);
 K1=zeros(1,Nkt);
 K2=zeros(1,Nkt);
 
 for i=1:Nkt
     ypos(nz)=-xvec(i);
     [nlyrsr,dlyrsxr, dlyrsyr]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
     [kk12,kk21] = ks(Fold,nz,dx,dy,dz,nlyrsr,dlyrsxr, dlyrsyr,  k0, ...
                    struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                           'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                          'operver','m2dpmloperR2','epsavg','no')); 
 K1(1,i)=kk12;                        
 K2(1,i)=kk21;
 end
 %
figure(3)
plot(xvec*1000,K1,'or','LineWidth',2)
hold on
plot(xvec*1000,K2,'ob','LineWidth',2)
plot(xvec*1000,zeros(1,Nkt),'k','LineWidth',2)
plot(-a*1000,K1(round(Nkt/2)),'ok')
%axis([min(xvec)*1000 max(xvec)*1000 min([min(K1),min(K2)]) max([max(K1),max(K2)])])
%axis([0 400 0 1 -0.6 0.6])
zlabel('Coupling')
ylabel('z (\mum)')
xlabel('tilt \rho (nm)')
legend('\kappa_{12}','\kappa_{21}')
drawnow
%}
ypos(nz)=a;
Fold=F; 
        %}
end
   %figure(10); plot((dz+dz/2:dz:nz*dz-dz/2),kl,'-ok'); drawnow;
 
%[nlyrs,dlyrsx, dlyrsy]=createmeshopt3d(dx,dy,dz,w1,w2,w3,gl,gr,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot);

%[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,1/dz); 

%[nlyrs,dlyrsx, dlyrsy]=createmeshopt3dbuumid(dx,dy,dz,xL,yL,wul,wdl,wur,wdr,dwU,dwD,lU,lD,gl,gr,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,'yes');
%{
nz=1;
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
[F,fa] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
Fold=F;                         
for nz=2:1/dz
[nlyrs,dlyrsx, dlyrsy]=createmesh3dfmm(dx,dy,dz,w1,w2,w3,gl,gr,sptw,yCladBot,hWg,yCladTop,nBelow,nCore,nAbove,ifplot,ypos,nz);  
[F,fa] = k12(Fold,nz,dx,dy,dz,nlyrs,dlyrsx, dlyrsy,  k0, ...
                     struct('mu_guess',3*k0,'NMODES_CALC',6,'tol',1e-16, ...
                            'BC',[0 0 0 0],'enginever','m2wcylR2', ...
                             'operver','m2dpmloperR2','epsavg','no'));  
Fold=F; 
            figure(8);
            semilogy(nz,abs(fa),'-pg'); hold on; drawnow
            figure(11);
            plot(nz,ypos(nz),'-pg'); hold on; drawnow
end
%}

%}


