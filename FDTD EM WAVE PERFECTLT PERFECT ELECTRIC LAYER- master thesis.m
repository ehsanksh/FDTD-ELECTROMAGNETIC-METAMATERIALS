%2D FDTD TE mode with a POINT source and a PEC abc

%     2-D FDTD TE code with PEC boundary conditions
%*********************************************************************** 
clear all
R=42;      %%%ROW
C=42;      %%%%COLUMN
h=42;
ic=21;jc=21;kc=21;     %%%CENTER OF YEE GRID%
dx=.000005;dt=dx/6e8;     %%%%TIME STEP
pi=3.14;epsz=8.8e-12;
c=3.0e+14;     %%%micrometer/sec%%%%%%
lambda=2*1.0e+01;      %%%%%%%WAVELENGHT IN MICROMETER%
w0=c/lambda;       %%%%%%FREQUENCY OF SOURCE%
E0=3.0;%%%%AMPLITUDE OF SOURCE
nmax=350;
etta=376.5;       %%%%IMPEDANCE OF FREE SPACE%
for j=12:C-1
    for i=1:R-1
        for k=1:h-1
             ez(i,j,k)=0;
             ex(i,j,k)=0;
             ey(i,j,k)=0;
             hz(i,j,k)=0;
             hx(i,j,k)=0;
             hy(i,j,k)=0;
             
        end
    end
end
rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/2,gcf,rect);

%%%%%%%%%%%%%% yee algoritm %%%%%%%%%%%%%%%%%%%
ez(ic,jc,kc)=E0*sin(2*pi*w0*.5*dt);

for n=1:nmax
    for i=1:R-1
    for j=1:C-2
        for k=1:h-2
            hx(i,j,k)=hx(i,j,k)+(.5/etta)*(ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k));
        end
    end
    end
     for i=1:R-2
    for j=1:C-1
        for k=1:h-2
            hy(i,j,k)=hy(i,j,k)+(.5/etta)*(ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k));
        end
    end
     end
        for j=1:C-2
    for i=1:R-2
        for k=1:h-1
            hz(i,j,k)=hz(i,j,k)+(.5/etta)*(ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k));
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
     for k=2:h-1
         for j=2:C-1
             for i=1:R-1
                 
                 ex(i,j,k)=ex(i,j,k)+.5*etta*(hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1));
                  
             end
         end
     end
      for k=2:h-1
         for j=1:C-1
             for i=2:R-1
                  
                 ey(i,j,k)=ey(i,j,k)+.5*etta*(hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k));
                  
             end
         end
      end
       for k=1:h-1
         for j=2:C-1
             for i=2:R-1
                   if i==21&&j==21&&k==21
            ez(i,j,k)=E0*sin(2*w0*pi*(n+.5)*dt);
            else
                 ez(i,j,k)=ez(i,j,k)+.5*etta*(hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k));
             end
         end
         end
       end
       if mod(n,2)==0;
      timestep=int2str(n);
tview(:,:)=ez(:,:,kc);
sview(:,:)=ez(:,jc,:);

subplot('position',[0.15 0.15 0.7 0.65]),pcolor(tview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['Ez(i,j,k=5), time step = ',timestep]);
xlabel('i coordinate');
ylabel('j coordinate');


nn=n/2;
M(:,nn)=getframe(gcf,rect);
 
       
       end
end
movie(gcf,M,0,10,rect);


 
       

    
    


            
           
            