%2D FDTD TE mode with a POINT source and a PML abc

%     2-D FDTD TE code with PML absorbing boundary conditions
%***********************************************************************
clear all
R=42;      %%%ROW
C=42;      %%%%COLUMN
h=42;
ic=21;jc=21;kc=21;     %%%CENTER OF YEE GRID%
dx=.00000005;dt=dx/6e8;     %%%%TIME STEP
nmax=200;
pi=3.14;epsz=8.8e-12;
c=3.0e+14;     %%%micrometer/sec%%%%%%
lambda=2*1.0e-01;      %%%%%%%WAVELENGHT IN MICROMETER%
w0=c/lambda;       %%%%%%FREQUENCY OF SOURCE%
E0=3.0e-01;       %%%%AMPLITUDE OF SOURCE%
etta=376.5;       %%%%IMPEDANCE OF FREE SPACE%
for j=1:C-1
    for i=1:R-1
        for k=1:h-1
             ez(i,j,k)=0;
             ex(i,j,k)=0;
             ey(i,j,k)=0;
             hz(i,j,k)=0;
             hx(i,j,k)=0;
             hy(i,j,k)=0;
             px(i,j,k)=0;
             py(i,j,k)=0;
             pz(i,j,k)=0;
             qx(i,j,k)=0;
             qy(i,j,k)=0;
             qz(i,j,k)=0;
             
        end
    end
end
%%%%%%%%%%%  pml paramaters    %%%%%%%%
for i=1:R
    ni1(i)=0;
    ni2(i)=1;
    ni3(i)=1;
    mi1(i)=0;
    mi2(i)=1;
    mi3(i)=1;
end
for j=1:C
    nj1(j)=0;
    nj2(j)=1;
    nj3(j)=1;
    mj1(j)=0;
    mj2(j)=1;
    mj3(j)=1;
end

for k=1:h
    nk1(k)=0;
    nk2(k)=1;
    nk3(k)=1;
    mk1(k)=0;
    mk2(k)=1;
    mk3(k)=1;
end
npml=1;       %%%%%%%%%%  NUMBER OF PML CELLS     %%%%%%%%%%%%%
for i=1:npml;
    IN=npml-i;
    xd=npml;
    IM=IN/xd;
    xn=.3*IM^3;     %%%%%%%%   xn=(function of electric conductivity in pml)*( dt/2*epsz)   %%%%%%%%
     mi1(i)=xn;
    mi1(R-1-i)=xn;
    ni2(i)=1/(1+xn);
    ni2(R-1-i)=1/(1+xn);
    ni3(i)=(1-xn)/(1+xn);
    ni3(R-i-1)=(1-xn)/(1+xn);
    IM=(IN)/xd;
    xn=.3*IM^3;
    ni1(i)=xn;
    ni1(R-1-i)=xn;
    mi2(i)=1/(1+xn);
    mi2(R-1-i)=1/(1+xn);
    mi3(i)=(1-xn)/(1+xn);
    mi3(R-1-i)=(1-xn)/(1+xn);
   
end
for j=1:npml;
    IN=npml-j;
    xd=npml;
    IM=IN/xd;
    xn=.3*IM^3;     %%%%%%%%   xn=(function of electric conductivity in pml)*( dt/2*epsz)   %%%%%%%%
     mj1(j)=xn;
    mj1(C-1-j)=xn;
    nj2(j)=1/(1+xn);
    nj2(C-1-j)=1/(1+xn);
    nj3(j)=(1-xn)/(1+xn);
    nj3(C-j-1)=(1-xn)/(1+xn);
    IM=(IN)/xd;
    xn=.3*IM^3;
    nj1(j)=xn;
    nj1(C-1-j)=xn;
    mj2(j)=1/(1+xn);
    mj2(C-1-j)=1/(1+xn);
    mj3(j)=(1-xn)/(1+xn);
    mj3(C-1-j)=(1-xn)/(1+xn);
end
for k=1:npml;
    IN=npml-k;
    xd=npml;
    IM=IN/xd;
    xn=.3*IM^3;     %%%%%%%%   xn=(function of electric conductivity in pml)*( dt/2*epsz)   %%%%%%%%
     mk1(k)=xn;
    mk1(h-1-k)=xn;
    nk2(k)=1/(1+xn);
    nk2(h-1-k)=1/(1+xn);
    nk3(k)=(1-xn)/(1+xn);
    nk3(h-k-1)=(1-xn)/(1+xn);
    IM=(IN-.5)/xd;
    xn=.3*IM^3;
    nk1(k)=xn;
    nk1(h-1-k)=xn;
    mk2(k)=1/(1+xn);
    mk2(h-1-k)=1/(1+xn);
    mk3(k)=(1-xn)/(1+xn);
    mk3(h-1-k)=(1-xn)/(1+xn);
end
rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/2,gcf,rect);


ez(ic,jc,kc)=E0*sin(2*pi*w0*.5*dt);



   nmax=100;    %%%%%%%%%%%%     NUMBER OF STEP   %%%%%%%%%%%%%
    for n=1:nmax
          for j=1:C-2
    for k=1:R-2
    for i=1:npml-1
            delta(i,j,k)=ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k);
            qx(i,j,k)=qx(i,j,k)+delta(i,j,k);
            hx(i,j,k)=mj3(j)*mk3(k)*hx(i,j,k)+mj2(j)*mk2(k)*.5*(delta(i,j,k)+mi1(i)*qx(i,j,k));
    end
    end
    end
    for j=1:C-2
    for k=1:R-2
    for i=npml:h-1-npml
             delta(i,j,k)=ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k);
           
            hx(i,j,k)=mj3(j)*mk3(k)*hx(i,j,k)+mj2(j)*mk2(k)*.5*delta(i,j,k);
    end
    end
    end
    for j=1:C-2
    for k=1:R-2
    for i=h-npml:h-2
            delta(i,j,k)=ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k);
            qx(i,j,k)=qx(i,j,k)+delta(i,j,k);
            hx(i,j,k)=mj3(j)*mk3(k)*hx(i,j,k)+mj2(j)*mk2(k)*.5*(delta(i,j,k)+mi1(i)*qx(i,j,k));
    end
    end
    end
       
       
   
    for k=1:C-2
    for i=1:R-2
    for j=1:npml-1
            delta(i,j,k)=ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k);
            qy(i,j,k)=qy(i,j,k)+delta(i,j,k);
            hy(i,j,k)=mi3(i)*mk3(k)*hy(i,j,k)+mi2(i)*mk2(k)*.5*(delta(i,j,k)+mj1(j)*qy(i,j,k));
    end
    end
    end
    for k=1:C-2
    for i=1:R-2
    for j=npml:h-1-npml
            delta(i,j,k)=ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k);
            hy(i,j,k)=mi3(i)*mk3(k)*hy(i,j,k)+mi2(i)*mk2(k)*.5*delta(i,j,k);
    end
    end
    end
           for k=1:C-2
    for i=1:R-2
        for j=h-npml:h-2
             delta(i,j,k)=ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k);
            qy(i,j,k)=qy(i,j,k)+delta(i,j,k);
            hy(i,j,k)=mi3(i)*mk3(k)*hy(i,j,k)+mi2(i)*mk2(k)*.5*(delta(i,j,k)+mj1(j)*qy(i,j,k));
            
        end
    end
           end
            for j=1:C-2
    for i=1:R-2
    for k=1:npml-1
            delta(i,j,k)=ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k);
            qz(i,j,k)=qz(i,j,k)+delta(i,j,k);
            hz(i,j,k)=mi3(i)*mj3(j)*hz(i,j,k)+mi2(i)*mj2(j)*.5*(delta(i,j,k)+mk1(k)*qz(i,j,k));
    end
    end
    end
         
    for  j=1:C-2
    for i=1:R-2
    for k=npml:h-1-npml
            delta(i,j,k)=ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k);
            hz(i,j,k)=mi3(i)*mj3(j)*hz(i,j,k)+mi2(i)*mj2(j)*.5*delta(i,j,k);
    end
    end
    end
    for j=1:C-2
    for i=1:R-2
    for k=h-npml:h-2
             delta(i,j,k)=ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k);
            qz(i,j,k)=qz(i,j,k)+delta(i,j,k);
            hz(i,j,k)=mi3(i)*mj3(j)*hz(i,j,k)+mi2(i)*mj2(j)*.5*(delta(i,j,k)+mk1(k)*qz(i,j,k));
    end
    end
    end
  
            for i=1:R-1
            for k=1:h-1
                ez(i,1,k)=0;ex(i,1,k)=0;ez(i,C-1,k)=0;ex(i,C-1,k)=0;
            end
            end
        for i=1:R-1
            for j=1:C-1
                ex(i,j,1)=0;ey(i,j,1)=0;ex(i,j,h-1);ey(i,j,h-1)=0;
            end
        end
        for j=1:C-1
            for k=1:h-1
                ey(1,j,k)=0;ez(1,j,k)=0;ey(R-1,j,k)=0;ez(R-1,j,k)=0;
            end
        end
           
    for j=2:C-1
    for k=2:R-1
    for i=1:npml-1
             delta(i,j,k)=hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1);
            px(i,j,k)=px(i,j,k)+delta(i,j,k);
            ex(i,j,k)=nj3(j)*nk3(k)*ex(i,j,k)+nj2(j)*nk2(k)*.5*(delta(i,j,k)+ni1(i)*px(i,j,k));
    end
    end
    end
      for j=2:C-1
      for k=2:R-1
      for i=npml:h-1-npml
            delta(i,j,k)=hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1);
            ex(i,j,k)=nj3(j)*nk3(k)*ex(i,j,k)+nj2(j)*nk2(k)*.5*delta(i,j,k);
      end
      end
      end
    for j=2:C-1
    for k=2:R-1
    for i=h-npml:h-2
            delta(i,j,k)=hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1);
            px(i,j,k)=px(i,j,k)+delta(i,j,k);
            ex(i,j,k)=nj3(j)*nk3(k)*ex(i,j,k)+nj2(j)*nk2(k)*.5*(delta(i,j,k)+ni1(i)*pz(i,j,k));
    end
    end
    end
    for k=2:C-1
    for i=2:R-1
    for j=1:npml-1:h-2
             delta(i,j,k)=hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k);
            py(i,j,k)=py(i,j,k)+delta(i,j,k);
            ey(i,j,k)=ni3(i)*nk3(k)*ey(i,j,k)+ni2(i)*nk2(k)*.5*(delta(i,j,k)+nj1(j)*py(i,j,k));
    end
    end
    end
    for k=2:C-1
    for i=2:R-1
    for j=npml:h-1-npml
      
            delta(i,j,k)=hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k);
            ey(i,j,k)=ni3(i)*nk3(k)*ey(i,j,k)+ni2(i)*nk2(k)*.5*delta(i,j,k);
    end
    end
    end
    for k=2:C-1
    for i=2:R-1
    for j=h-npml:h-2
            delta(i,j,k)=hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k);
            py(i,j,k)=py(i,j,k)+delta(i,j,k);
            ey(i,j,k)=ni3(j)*nk3(k)*ey(i,j,k)+ni2(i)*nk2(k)*.5*(delta(i,j,k)+nj1(j)*pz(i,j,k));
    end
    end
    end
       
        
    for j=2:C-1
    for i=2:R-1
    for k=1:npml-1
            delta(i,j,k)=hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k);
            pz(i,j,k)=pz(i,j,k)+delta(i,j,k);
            ez(i,j,k)=ni3(i)*nj3(j)*ez(i,j,k)+ni2(i)*nj2(j)*.5*(delta(i,j,k)+nk1(k)*pz(i,j,k));
    end
    end
    end
    for j=2:C-1
    for i=2:R-1
    for k=npml:h-1-npml
             if i==ic&&j==jc&&k==kc
            ez(i,j,k)=E0*sin(2*w0*pi*(n+.5)*dt);
             else
            delta(i,j,k)=hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k);
            ez(i,j,k)=ni3(i)*nj3(j)*ez(i,j,k)+ni2(i)*nj2(j)*.5*delta(i,j,k);
             end
    end
    end
    end
    for j=2:C-1
    for i=2:R-1
    for k=h-npml:h-2
            delta(i,j,k)=hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k);
            pz(i,j,k)=pz(i,j,k)+delta(i,j,k);
            ez(i,j,k)=ni3(i)*nj3(j)*ez(i,j,k)+ni2(i)*nj2(j)*.5*(delta(i,j,k)+nk1(k)*pz(i,j,k));
    end
    end
    end
    
  
          if mod(n,2)==0;
      timestep=int2str(n);
tview(:,:)=ez(:,:,kc);
sview(:,:)=ez(:,jc,:);

subplot('position',[0.15 0.45 0.7 0.45]),pcolor(tview');
shading flat;
caxis([-.0001 .0001]);
colorbar;
axis image;
title(['Ez(i,j,k=kc), time step = ',timestep]);
xlabel('i coordinate');
ylabel('j coordinate');

nn=n/2;
M(:,nn)=getframe(gcf,rect);
 
       
       end
end
movie(gcf,M,20,15,rect);


  
    
            
          


   

