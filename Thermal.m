classdef Thermal < handle 
    %Thermal Summary of this class goes here
    %   Detailed explanation goes here
    %Thermal is used to solve the Poisson equation.
    %The Thermal object will do the matrix assembly once it is created. The
    %member function 'ThermalSolver' is used to solve the Thermal equation
    %with approriate input arguments. 
    properties
        K
        B
        Device
        FEMgrid
    end
    
    methods
        function Thermalobj = Thermal(device,FEMgrid)
            Thermalobj.Device = device;
            Thermalobj.FEMgrid = FEMgrid;
            MatrixAssembly(Thermalobj,device,FEMgrid);
        end %end of constructor
        
        function MatrixAssembly(Thermalobj,device,FEMgrid)
            
            Thermalobj.K = sparse(FEMgrid.N_n,FEMgrid.N_n);
            Thermalobj.B = sparse(FEMgrid.N_n,1);
            
            for ii_e = 1 : FEMgrid.N_e
               
                if FEMgrid.Element(ii_e).material == device.channel_material
                    themcond = device.themcond_channel;
                elseif FEMgrid.Element(ii_e).material == device.oxide_material
                    themcond = device.themcond_oxide;
                end
                
                n1 = FEMgrid.Element(ii_e).n1;
                n2 = FEMgrid.Element(ii_e).n2;
                n3 = FEMgrid.Element(ii_e).n3;
                
                De = [themcond 0;0 themcond];
                Ktte = FEMgrid.Element(ii_e).Be.'*De*FEMgrid.Element(ii_e).Be*FEMgrid.Element(ii_e).A;
                
                Thermalobj.K(n1,n1)=Thermalobj.K(n1,n1)+Ktte(1,1);
                Thermalobj.K(n1,n2)=Thermalobj.K(n1,n2)+Ktte(1,2);
                Thermalobj.K(n1,n3)=Thermalobj.K(n1,n3)+Ktte(1,3);
                Thermalobj.K(n2,n1)=Thermalobj.K(n2,n1)+Ktte(2,1);
                Thermalobj.K(n2,n2)=Thermalobj.K(n2,n2)+Ktte(2,2);
                Thermalobj.K(n2,n3)=Thermalobj.K(n2,n3)+Ktte(2,3);
                Thermalobj.K(n3,n1)=Thermalobj.K(n3,n1)+Ktte(3,1);
                Thermalobj.K(n3,n2)=Thermalobj.K(n3,n2)+Ktte(3,2);
                Thermalobj.K(n3,n3)=Thermalobj.K(n3,n3)+Ktte(3,3);
            end
            
        end %end of MatrixAssembly
        
        function [T] = ThermalSolver(Thermalobj,U,Ne,Ph,Fn,Fp,T_old)
            
            [Heat_Source] = Joule_Heating(Thermalobj,U,Ne,Ph,Fn,Fp,T_old);
            
            for ii_n = 1 : Thermalobj.FEMgrid.N_n
                if Thermalobj.FEMgrid.Node(ii_n).bm == Thermalobj.Device.gate_bm 
                    Thermalobj.K(ii_n,:) = sparse(1,Thermalobj.FEMgrid.N_n); 
                    Thermalobj.K(ii_n,ii_n) = 1;
                    Heat_Source(ii_n,1)=0;
                    Thermalobj.B(ii_n,1) = Thermalobj.Device.tem_boundary;
                elseif Thermalobj.FEMgrid.Node(ii_n).bm == Thermalobj.Device.drain_bm 
                    Thermalobj.K(ii_n,:) = sparse(1,Thermalobj.FEMgrid.N_n);
                    Thermalobj.K(ii_n,ii_n) = 1;
                    Heat_Source(ii_n,1)=0;
                    Thermalobj.B(ii_n,1) = Thermalobj.Device.tem_boundary;
                elseif Thermalobj.FEMgrid.Node(ii_n).bm == Thermalobj.Device.source_bm 
                    Thermalobj.K(ii_n,:) = sparse(1,Thermalobj.FEMgrid.N_n);
                    Thermalobj.K(ii_n,ii_n) = 1;  
                    Heat_Source(ii_n,1)=0;
                    Thermalobj.B(ii_n,1) = Thermalobj.Device.tem_boundary;
                end
            end

            T = Thermalobj.K\(Heat_Source+Thermalobj.B);
        end %end of ThermalSolver
        
        function [Heat_Source] = Joule_Heating(Thermalobj,U,Ne,Ph,Fn,Fp,T_old) 
            
            factor = Thermalobj.Device.md_factor;
            
            Heat_Source = sparse(Thermalobj.FEMgrid.N_n,1);
            
            for ii_e = 1 : Thermalobj.FEMgrid.N_e
                
                n1 = Thermalobj.FEMgrid.Element(ii_e).n1;
                n2 = Thermalobj.FEMgrid.Element(ii_e).n2;
                n3 = Thermalobj.FEMgrid.Element(ii_e).n3;
                
                Ue = [U(n1) U(n2) U(n3)].';
                Nee = [Ne(n1) Ne(n2) Ne(n3)].';
                Phe = [Ph(n1) Ph(n2) Ph(n3)].';
                Fne = [Fn(n1) Fn(n2) Fn(n3)].';
                Fpe = [Fp(n1) Fp(n2) Fp(n3)].';
                Te = [T_old(n1) T_old(n2) T_old(n3)].';
                
                
                tem = (Te(1)+Te(2)+Te(3))/3;
                Fe = Thermalobj.Device.q*Thermalobj.Device.mu*(tem/300)^(factor)*Ue.'*Thermalobj.FEMgrid.Element(ii_e).Be.'*(-Thermalobj.FEMgrid.Element(ii_e).Be)*Fne*Thermalobj.FEMgrid.Element(ii_e).Se*Nee + ...
                     Thermalobj.Device.q*Thermalobj.Device.muv*(tem/300)^(factor)*Ue.'*Thermalobj.FEMgrid.Element(ii_e).Be.'*(-Thermalobj.FEMgrid.Element(ii_e).Be)*Fpe*Thermalobj.FEMgrid.Element(ii_e).Se*Phe;
                
                Heat_Source(n1,1)=Heat_Source(n1,1)+Fe(1,1);
                Heat_Source(n2,1)=Heat_Source(n2,1)+Fe(2,1);
                Heat_Source(n3,1)=Heat_Source(n3,1)+Fe(3,1);
            end
        end %end of Joule_Heating
        
%         function MatrixAssembly(Thermalobj,device,FEMgrid)
%             
%             N_n = length(FEMgrid.Node);
%             N_e = length(FEMgrid.Element);
%             
%             Thermalobj.K = sparse(N_n,N_n);
%             Thermalobj.S = sparse(N_n,N_n);
%             Thermalobj.B = sparse(N_n,1);
%             rbar = 1;
%             
%             for ii_e = 1 : N_e
%                 
%                 if FEMgrid.Element(ii_e).material == device.channel_material
%                     themcond = device.themcond_channel;
%                 elseif FEMgrid.Element(ii_e).material == device.oxide_material
%                     themcond = device.themcond_oxide;
%                 end
%                 
%                 if FEMgrid.Element(ii_e).material ~= device.oxide_material
%                     %%%%%%%% compute the overlap matrix %%%%%%%%%%%
%                     Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)=Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/6;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)=Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)=Thermalobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)=Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)=Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/6;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)=Thermalobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)=Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/12; 
%                     Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)=Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/12; 
%                     Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)=Thermalobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/6;     
%                 end
%                 %%%%%%%% compute the kinetic energy matrix %%%%%%%%%%%
%                 Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1) = Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x3)^2+...
%                                                (FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)^2)*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2) = Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3) = Thermalobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1) = Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2) = Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)^2+...
%                                                (FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y3)^2)*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3) = Thermalobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1) = Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2) = Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3) = Thermalobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)+(themcond/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x2)^2+...
%                                                (FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2)^2)*rbar;   
%             end
%             %Thermalobj.K = -Thermalobj.K;
%             
%             for ii_n = 1 : N_n
%                 if FEMgrid.Node(ii_n).bm == device.gate_bm 
%                     Thermalobj.K(ii_n,:) = sparse(1,N_n); 
%                     Thermalobj.S(ii_n,:) = sparse(1,N_n);
%                     %Thermalobj.Sad(ii_n) = 0;
%                     Thermalobj.K(ii_n,ii_n) = 1;
%                     Thermalobj.B(ii_n,1) = device.tem_boundary;
%                 elseif FEMgrid.Node(ii_n).bm == device.drain_bm 
%                     Thermalobj.K(ii_n,:) = sparse(1,N_n);
%                     Thermalobj.S(ii_n,:) = sparse(1,N_n);
%                     %Thermalobj.Sad(ii_n) = 0;
%                     Thermalobj.K(ii_n,ii_n) = 1;
%                     Thermalobj.B(ii_n,1) = device.tem_boundary;
%                 elseif FEMgrid.Node(ii_n).bm == device.source_bm 
%                     Thermalobj.K(ii_n,:) = sparse(1,N_n);
%                     Thermalobj.S(ii_n,:) = sparse(1,N_n);
%                     %Thermalobj.Sad(ii_n) = 0;
%                     Thermalobj.K(ii_n,ii_n) = 1;     
%                     Thermalobj.B(ii_n,1) = device.tem_boundary;
%                 end
%             end
%         end %end of MatrixAssembly

%         function [T] = ThermalSolver(Thermalobj,heat_source,Xlin,Ylin)
%             X = [Thermalobj.FEMgrid.Node.x];
%             Y = [Thermalobj.FEMgrid.Node.y]';
%             [xx,yy] = meshgrid(Xlin,Ylin);
%             heat_source = full(heat_source);
%             heat_mat = interp2(xx,yy,heat_source,X,Y);
%             heat_node = diag(heat_mat);
%             T = Thermalobj.K\(Thermalobj.S*heat_node+Thermalobj.B);
%         end %end of thermal_equation
        
    end %end of methods
    
end %end of classdef

