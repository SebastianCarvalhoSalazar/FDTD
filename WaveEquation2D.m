clear;close all;clc;

% Domain/Space

Lx = 20;
Ly = 20;
dx = 0.1;
dy = dx;
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);

% Time
T = 16; 

% Field Variable

wn = zeros(nx,ny);
wnm1 = wn; % w at time n-1
wnp1 = wn; % w at time n+1

% Parameters (CLF = c*dt/dx)

CLF = 0.5; 
c = 1;
dt = (CLF*dx)/c;

% Time Stepping Loop

t = 0;
f = 15;

fh = figure(1);
set(fh, 'Color', 'white'); 
set(fh,'Position',[8 45 1350 630]);

aviobj = VideoWriter('2D_interference.avi');
aviobj.Quality = 100;
aviobj.FrameRate = 10;
open(aviobj);

while(t < T)
    
    % Reflecting Boundary Conditions
    
%     wn(:,[1 end]) = 0;
%     wn([1 end],:) = 0;

    % Absorbing boundary conditions
    
    wnp1(1,:) = wn(2,:)+((CLF-1)/(CLF+1))*(wnp1(2,:)-wn(1,:));
    wnp1(end,:) = wn(end-1,:)+((CLF-1)/(CLF+1))*(wnp1(end-1,:)-wn(end,:));
    wnp1(:,1) = wn(:,2)+((CLF-1)/(CLF+1))*(wnp1(:,2)-wn(:,1));
    wnp1(:,end) = wn(:,end-1)+((CLF-1)/(CLF+1))*(wnp1(:,end-1)-wn(:,end));
    
    % Solution
    
    t = t+dt;
    wnm1 = wn; wn = wnp1;
    
    % Source
    A = 260;
    
    wn(100,60*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,59*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,58*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,57*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,56*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,55*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,54*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,53*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,52*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,51*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,50*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,49*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,48*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,47*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,46*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,45*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,44*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,43*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,42*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,41*2) = dt^2*A*sin(2*pi*f*t/T);
    wn(100,40*2) = dt^2*A*sin(2*pi*f*t/T);

    
    for i = 2:nx-1 
        for j = 2:ny-1
            wnp1(i,j) = 2*wn(i,j) - wnm1(i,j) ...
                + CLF^2 * (wn(i+1,j) + wn(i,j+1) - 4*wn(i,j) + wn(i-1,j) + wn(i,j-1));
        end
    end
   
    
    % Plot 
   
    subplot(1,2,1);
    pcolor(x, y, wn');
    colormap jet;
    title1 = '          2D VIEW          SOURCES: 21       ';
    title(title1 ,'fontsize',14);    
    colorbar('location','eastoutside','fontsize',14); 
    colorbar('YLim',[-1 1],'fontsize',14);
    shading interp
    
%     hold on;
%     scatter(10,10.1,20,'kd','filled');
%     hold off;

    subplot(1,2,2);
    surf(x, y, wn','EdgeColor','none','LineStyle','none','FaceLighting','phong'); 
    title2 = [' TIME = ',sprintf('%.2f',t), ' SECONDS', '             3D VIEW'];
    title(title2 ,'fontsize',14);   shading interp;
    colorbar('location','eastoutside','fontsize',14);
    colorbar('YLim',[-1 1],'fontsize',14);
    axis([0 Lx 0 Ly -1 1]); 
    grid minor;
    shg; pause(0.01);

    F = getframe(fh);
    writeVideo(aviobj,F); 
    
end

close(fh);
aviobj.close();