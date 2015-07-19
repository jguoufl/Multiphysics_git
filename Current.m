classdef Current < handle & CommonFunctions
    %Current Summary of this class goes here
    %   Detailed explanation goes here
    %Current is used to calculate the electrical current and heat
    %generation in the device. Since the original grid is FEM grid, it will
    %do a uniform interpolation for the charge density, potential and
    %temperature profiles. It returns current(Id), heat generation
    %(Joule_heat), and the coordiates of the uniform interpolation grid.
    properties
        Device
    end
    
    methods
        
        function Currentobj = Current(device)
            Currentobj = Currentobj@CommonFunctions();
            Currentobj.Device = device;
        end
        
        function [Id, Joule_heat, Xlin, Ylin] = Current_Heat(Currentobj,FEMgrid,Ne_bias,Ph_bias,U,T)
            
            Dn = Currentobj.Device.mu * Currentobj.Device.kBT;
            Dnv = Currentobj.Device.muv * Currentobj.Device.kBT;
            %mobility degradation factor
            factor = Currentobj.Device.md_factor;
            
            Ne = full(Ne_bias);
            Ph = full(Ph_bias);
            
            %convert electron and hole density and potential from the FEM grid
            %to a uniform grid
            vis(:,1) = [FEMgrid.Node.x]';
            vis(:,2) = [FEMgrid.Node.y]';
            vis(:,3) = log(Ne);
            [xlin ylin Ne2D] = uniform_interpolation(Currentobj,vis,Currentobj.Device.LScale);
            Ne2D = exp(Ne2D);
            
            vis(:,1) = [FEMgrid.Node.x]';
            vis(:,2) = [FEMgrid.Node.y]';
            vis(:,3) = log(Ph);
            [xlin ylin Ph2D] = uniform_interpolation(Currentobj,vis,Currentobj.Device.LScale);
            Ph2D = exp(Ph2D);
            
            vis(:,1) = [FEMgrid.Node.x]';
            vis(:,2) = [FEMgrid.Node.y]';
            vis(:,3) = U;
            [xlin ylin U2D] = uniform_interpolation(Currentobj,vis,Currentobj.Device.LScale);
            
            vis(:,1) = [FEMgrid.Node.x]';
            vis(:,2) = [FEMgrid.Node.y]';
            vis(:,3) = T;
            [xlin ylin T2D] = uniform_interpolation(Currentobj,vis,Currentobj.Device.LScale);
            
            Xlin = xlin;
            Ylin = ylin;
            
            Nx = length(xlin);
            Ny = length(ylin);
            del_x = xlin(2) - xlin(1);
            del_y = ylin(2) - ylin(1);
            
            channel_segment = floor(Currentobj.Device.TChannel / (Currentobj.Device.TOxide + Currentobj.Device.TChannel) * (Ny - 1));

            %calculate the current in the middle cross-section of the channel
            for ii_n = Ny : -1 : (Ny - channel_segment + 1)
                %electron current
                Efield = -((U2D(ii_n,floor(Nx/2))+U2D(ii_n-1,floor(Nx/2)))/2 - (U2D(ii_n,floor(Nx/2)+1)+U2D(ii_n-1,floor(Nx/2)+1))/2)/del_x;
                NeEfield = ((Ne2D(ii_n, floor(Nx/2))+Ne2D(ii_n-1,floor(Nx/2)))/2 + (Ne2D(ii_n,floor(Nx/2)+1)+Ne2D(ii_n-1,floor(Nx/2)+1))/2)/2;
                Jedrift = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_n,floor(Nx/2))/300)^(factor)*Efield*NeEfield*del_y;
                Negradient = -((Ne2D(ii_n,floor(Nx/2))+Ne2D(ii_n-1,floor(Nx/2)))/2 - (Ne2D(ii_n,floor(Nx/2)+1)+Ne2D(ii_n-1,floor(Nx/2)+1))/2)/del_x;
                Jediff = Currentobj.Device.q*Dn*(T2D(ii_n,floor(Nx/2))/300)^(factor)*Negradient*del_y;
                Je(ii_n) = Jedrift + Jediff;
                %hole current
                PhEfield = ((Ph2D(ii_n, floor(Nx/2))+Ph2D(ii_n-1,floor(Nx/2)))/2 + (Ph2D(ii_n,floor(Nx/2)+1)+Ph2D(ii_n-1,floor(Nx/2)+1))/2)/2;
                Jpdrift = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_n,floor(Nx/2))/300)^(factor)*Efield*PhEfield*del_y;
                Phgradient = -((Ph2D(ii_n,floor(Nx/2))+Ph2D(ii_n-1,floor(Nx/2)))/2 - (Ph2D(ii_n,floor(Nx/2)+1)+Ph2D(ii_n-1,floor(Nx/2)+1))/2)/del_x;
                Jpdiff = -Currentobj.Device.q*Dnv*(T2D(ii_n,floor(Nx/2))/300)^(factor)*Phgradient*del_y;
                Jp(ii_n) = Jpdrift + Jpdiff;
            end

            In = sum(Je);
            Ip = sum(Jp);
            Id = -(In + Ip);
            
            %calculate the Joule heat
            Joule_heat = sparse(Ny,Nx);
            for ii_x = 1 : Nx
                for ii_y = Ny : -1 : (Ny - channel_segment) 
                    if ii_y == Ny
                        if ii_x == 1
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
            
                            Ey = (U2D(ii_y,ii_x)-U2D(ii_y-1,ii_x))/del_y;
                            Negry = (Ne2D(ii_y,ii_x)-Ne2D(ii_y-1,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y,ii_x)-Ph2D(ii_y-1,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        elseif ii_x == Nx
                            Ex = (U2D(ii_y,ii_x)-U2D(ii_y,ii_x-1))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x)-Ne2D(ii_y,ii_x-1))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*Ex*(T2D(ii_y,ii_x)/300)^(factor)*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x)-Ph2D(ii_y,ii_x-1))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*Ex*(T2D(ii_y,ii_x)/300)^(factor)*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
             
                            Ey = (U2D(ii_y,ii_x)-U2D(ii_y-1,ii_x))/del_y;
                            Negry = (Ne2D(ii_y,ii_x)-Ne2D(ii_y-1,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y,ii_x)-Ph2D(ii_y-1,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        else
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x-1))/(2*del_x);
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x-1))/(2*del_x);
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x-1))/(2*del_x);
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
             
                            Ey = (U2D(ii_y,ii_x)-U2D(ii_y-1,ii_x))/del_y;
                            Negry = (Ne2D(ii_y,ii_x)-Ne2D(ii_y-1,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y,ii_x)-Ph2D(ii_y-1,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        end
                    elseif ii_y == Ny - channel_segment
                        if ii_x == 1
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
             
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y,ii_x))/del_y;
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        elseif ii_x == Nx
                            Ex = (U2D(ii_y,ii_x)-U2D(ii_y,ii_x-1))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x)-Ne2D(ii_y,ii_x-1))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x)-Ph2D(ii_y,ii_x-1))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
              
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y,ii_x))/del_y;
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        else
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x-1))/(2*del_x);
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x-1))/(2*del_x);
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x-1))/(2*del_x);
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
              
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y,ii_x))/del_y;
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y,ii_x))/del_y;
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y,ii_x))/del_y;
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        end
                    else
                        if ii_x == 1
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
            
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y-1,ii_x))/(2*del_y);
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y-1,ii_x))/(2*del_y);
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y-1,ii_x))/(2*del_y);
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        elseif ii_x == Nx
                            Ex = (U2D(ii_y,ii_x)-U2D(ii_y,ii_x-1))/del_x;
                            Negrx = (Ne2D(ii_y,ii_x)-Ne2D(ii_y,ii_x-1))/del_x;
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x)-Ph2D(ii_y,ii_x-1))/del_x;
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
             
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y-1,ii_x))/(2*del_y);
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y-1,ii_x))/(2*del_y);
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y-1,ii_x))/(2*del_y);
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        else
                            Ex = (U2D(ii_y,ii_x+1)-U2D(ii_y,ii_x-1))/(2*del_x);
                            Negrx = (Ne2D(ii_y,ii_x+1)-Ne2D(ii_y,ii_x-1))/(2*del_x);
                            Jnx = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negrx;
                            Phgrx = (Ph2D(ii_y,ii_x+1)-Ph2D(ii_y,ii_x-1))/(2*del_x);
                            Jpx = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ex*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgrx;
                            Jhx = (Jnx+Jpx)*Ex;
             
                            Ey = (U2D(ii_y+1,ii_x)-U2D(ii_y-1,ii_x))/(2*del_y);
                            Negry = (Ne2D(ii_y+1,ii_x)-Ne2D(ii_y-1,ii_x))/(2*del_y);
                            Jny = Currentobj.Device.q*Currentobj.Device.mu*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ne2D(ii_y,ii_x)+Currentobj.Device.q*Dn*(T2D(ii_y,ii_x)/300)^(factor)*Negry;
                            Phgry = (Ph2D(ii_y+1,ii_x)-Ph2D(ii_y-1,ii_x))/(2*del_y);
                            Jpy = Currentobj.Device.q*Currentobj.Device.muv*(T2D(ii_y,ii_x)/300)^(factor)*Ey*Ph2D(ii_y,ii_x)-Currentobj.Device.q*Dnv*(T2D(ii_y,ii_x)/300)^(factor)*Phgry;
                            Jhy = (Jny+Jpy)*Ey;
                            Joule_heat(ii_y,ii_x)=Jhx+Jhy;
                        end
                    end
                end
            end
        end %end of current_heat
       
    end %end of method
    
end %end of class

