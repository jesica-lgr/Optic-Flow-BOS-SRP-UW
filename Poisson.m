function [Po]=Poisson(Dm,lx,ly,m,t,pp,a)

%------------------------------------------------------------------------
%---Poisson Equation Solver using Successive Overrelaxation Method-------
%------------------------------------------------------------------------


%-------------------------------------------------------------------------
%--Define dimension of the 2-d box and -----------------------------------
%-------------------------------------------------------------------------

loop=1;
        
% Nx=a*floor(ly/m)-a+1;
% Ny=a*floor(lx/m)-a+1;
Nx=size(Dm,1);
Ny=size(Dm,2);

% M=40000; % maximum iteration value
M=10000000;
%--------------------------------------------------------------------------
%--Initialize V-matrix and Rho-matrix--------------------------------------
%---V(i,j)=potential inside the 2-D box------------------------------------
%---rho(i,j)=charge density inside the 2-D box-----------------------------
%---i is along X-axis and j is along Y-axis--------------------------------
%--------------------------------------------------------------------------
Po(1:Nx,1:Ny)=0.0;
des=Dm;
% rho(25,25)=25; % charge at the center of the box

%--------------------------------------------------------------------------
% Condiciones de Frontera
% t=1.0002926;
% t=1.29;

Po(1,:)=t; 
Po(Nx,:)=t; 
Po(:,1)=t; 
Po(:,Ny)=t;
%-------------------------------------------------------------------------

w=cos(pi/Nx)+cos(pi/Ny); % Converging Term
Ncount=0;
loop=1;

%%Escoge las fronteras de la condición de Neumann
choice = questdlg('¿En qué fronteras deseas tener Neumann?', ...
	'Fronteras-Neumann', ...
	'TODAS','Este-Oeste');
% Handle response
switch choice
%     case 'Norte-Sur';
%         disp([choice 'Las fronteras Norte y Sur'])
%         
%         while loop==1;
%             Rmin=0; 
%             for i=2:Nx0
%                 for j=2:Ny-1
%                     Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+Po(i,j-1)+Po(i,j+1)+((1/(Nx+1))*(1/(Ny+1)))*des(i,j))-Po(i,j));
%                     Rmin=Rmin+abs(Residue);
%                     Po(i,j)=Po(i,j)+Residue;
%                 end
%             end
% %         
% %         u(1,14:42)=u(2,14:42);%NeumannNorth
% %         u(Nx,14:42)=u(Nx-1,14:42);%NeumannSouth
%             Po(1,:)=Po(2,:);%NeumannNorth
%             Po(Nx,:)=Po(Nx-1,:);%NeumannSouth
%         
%             Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
% 
% %           ConvC=0.00001;
%             ConvC=0.0001;
%             
%             if(Rmin>=ConvC)
%                 Ncount=Ncount+1;
%                 if mod(Ncount,5000) == 0
%                     Jesy = 18;
%                 end
%                 
%                 if(Ncount>M)
%                     loop=0;
%                     disp(['solution doesnt converge in ',num2str(M),' iter'])
%                 end
%             else
%                 loop=0;
%                 disp(['solution converges in ',num2str(Ncount) ,' iteration'])
%             end
%             
%         end
        
        
        
        
        
    case 'Este-Oeste';
        disp([choice 'Las fronteras Este y Oeste'])
        
        while loop==1;
            Rmin=0; 
            for i=2:Nx-1
                for j=2:Ny-1
                    Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+Po(i,j-1)+Po(i,j+1)+((1/(Nx+1))^2)*des(i,j))-Po(i,j));
                    Rmin=Rmin+abs(Residue);
                    Po(i,j)=Po(i,j)+Residue;
                end
            end
        
            Po(:,1)=Po(:,2);%NeumannWest
            Po(:,Ny)=Po(:,Ny-1);%NeumannEast

            Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
            ConvC=0.0001;

            if(Rmin>=ConvC)
                Ncount=Ncount+1;
                if(Ncount>M)
                    loop=0;
                    disp(['solution doesnt converge in ',num2str(M),' iter'])
                end
            else
                loop=0;
                disp(['solution converges in ',num2str(Ncount) ,' iteration'])
            end
        end
        
   
    case 'TODAS';
        disp('ALL')
        
        while loop==1;
            Rmin=0; 
            for i=2:Nx-1
                for j=2:Ny-1
                    Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+(Po(i,j-1)+Po(i,j+1))+((1/(Nx+1))^2)*des(i,j))-Po(i,j));
                    Rmin=Rmin+abs(Residue);
                    Po(i,j)=Po(i,j)+Residue;
                end
            end
        
            Po(:,1)=Po(:,2);%NeumannWest
            Po(:,Ny)=Po(:,Ny-1);%NeumannEast
            Po(1,:)=Po(2,:);%NeumannNorth
            Po(Nx,:)=Po(Nx-1,:);%NeumannSouth

            Rmin=Rmin/(Nx*Ny); % Average Residue per grid point

            %ConvC=0.00001;
            %ConvC=0.0001;
            %ConvC=0.001;
            %ConvC=0.1;
            ConvC=0.4;
            if(Rmin>=ConvC)
                Ncount=Ncount+1;
                if(Ncount>M)
                    loop=0;
                    disp(['solution doesnt converge in ',num2str(M),' iter ',num2str(Rmin), 'Rmin  '])
                end
            else
                loop=0;
                disp(['solution converges in ',num2str(Ncount) ,' iteration'])
            end
        end
        
end %switch


%------------------------------------------------------------------------
% Plot the result
%------------------------------------------------------------------------

X=1:Nx;
Y=1:Ny;
% Rmin
%  contour(X,Y,V,0:1:50,'linewidth',2)
%  h=gca; 
% get(h,'FontSize') 
% set(h,'FontSize',12)
% colorbar('location','eastoutside','fontsize',12);
% axis([1 Nx 1 Ny])
% xlabel('X','fontSize',12);
% ylabel('Y','fontSize',12);
% title('Electric Potential Distribution, V(X,Y)','fontsize',12);
% % figure(pp+3);
% % if pp == 1;
% %     surf(X,Y,u')
% %     title('Densidad')
% %     
% % else
% % 
% %     surf(X,Y,u')
% %     title('Índice de Refracción')
end