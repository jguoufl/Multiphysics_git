classdef Mechanical < handle & CommonFunctions
    
    properties
        K
        F
        B
        Them_Stress_T0
        Device
        FEMgrid
    end
    
    methods
        function Mechanicalobj = Mechanical(device, FEMgrid)
            
            Mechanicalobj = Mechanicalobj@CommonFunctions();
            
            Mechanicalobj.Device = device;
            Mechanicalobj.FEMgrid = FEMgrid;
            
            MatrixAssembly(Mechanicalobj,device,FEMgrid);
        end %end of constructor
        
        function MatrixAssembly(Mechanicalobj,device,FEMgrid)
            
            Mechanicalobj.K = sparse(2*FEMgrid.N_n,2*FEMgrid.N_n);
            Mechanicalobj.F = sparse(2*FEMgrid.N_n,1);
            Mechanicalobj.Them_Stress_T0 = sparse(2*FEMgrid.N_n,1);
            
            T0 = 300;
            for ii_e = 1 : FEMgrid.N_e
                
                if FEMgrid.Element(ii_e).material == device.channel_material
                    elast = device.elast_channel;
                    themexpan = device.themexpan_channel;
                elseif FEMgrid.Element(ii_e).material == device.oxide_material
                    elast = device.elast_oxide;
                    themexpan = device.themexpan_oxide;
                end
                
                Kdde = FEMgrid.Element(ii_e).Be_M.'*elast*FEMgrid.Element(ii_e).Be_M*FEMgrid.Element(ii_e).A;
                Them_Stress_T0e = FEMgrid.Element(ii_e).Be_M.'*elast*themexpan*T0*FEMgrid.Element(ii_e).A;
                
                n1 = FEMgrid.Element(ii_e).n1;
                n2 = FEMgrid.Element(ii_e).n2;
                n3 = FEMgrid.Element(ii_e).n3;
                Mechanicalobj.K(2*n1-1:2*n1,2*n1-1:2*n1) = Mechanicalobj.K(2*n1-1:2*n1,2*n1-1:2*n1)+[Kdde(1,1) Kdde(1,2);Kdde(2,1) Kdde(2,2)];
                Mechanicalobj.K(2*n1-1:2*n1,2*n2-1:2*n2) = Mechanicalobj.K(2*n1-1:2*n1,2*n2-1:2*n2)+[Kdde(1,3) Kdde(1,4);Kdde(2,3) Kdde(2,4)];
                Mechanicalobj.K(2*n1-1:2*n1,2*n3-1:2*n3) = Mechanicalobj.K(2*n1-1:2*n1,2*n3-1:2*n3)+[Kdde(1,5) Kdde(1,6);Kdde(2,5) Kdde(2,6)];
                Mechanicalobj.K(2*n2-1:2*n2,2*n1-1:2*n1) = Mechanicalobj.K(2*n2-1:2*n2,2*n1-1:2*n1)+[Kdde(3,1) Kdde(3,2);Kdde(4,1) Kdde(4,2)];
                Mechanicalobj.K(2*n2-1:2*n2,2*n2-1:2*n2) = Mechanicalobj.K(2*n2-1:2*n2,2*n2-1:2*n2)+[Kdde(3,3) Kdde(3,4);Kdde(4,3) Kdde(4,4)];
                Mechanicalobj.K(2*n2-1:2*n2,2*n3-1:2*n3) = Mechanicalobj.K(2*n2-1:2*n2,2*n3-1:2*n3)+[Kdde(3,5) Kdde(3,6);Kdde(4,5) Kdde(4,6)];
                Mechanicalobj.K(2*n3-1:2*n3,2*n1-1:2*n1) = Mechanicalobj.K(2*n3-1:2*n3,2*n1-1:2*n1)+[Kdde(5,1) Kdde(5,2);Kdde(6,1) Kdde(6,2)];
                Mechanicalobj.K(2*n3-1:2*n3,2*n2-1:2*n2) = Mechanicalobj.K(2*n3-1:2*n3,2*n2-1:2*n2)+[Kdde(5,3) Kdde(5,4);Kdde(6,3) Kdde(6,4)];
                Mechanicalobj.K(2*n3-1:2*n3,2*n3-1:2*n3) = Mechanicalobj.K(2*n3-1:2*n3,2*n3-1:2*n3)+[Kdde(5,5) Kdde(5,6);Kdde(6,5) Kdde(6,6)];
                
                Mechanicalobj.Them_Stress_T0(2*n1-1:2*n1,:)=Mechanicalobj.Them_Stress_T0(2*n1-1:2*n1,:)+Them_Stress_T0e(1:2,:);
                Mechanicalobj.Them_Stress_T0(2*n2-1:2*n2,:)=Mechanicalobj.Them_Stress_T0(2*n2-1:2*n2,:)+Them_Stress_T0e(3:4,:);
                Mechanicalobj.Them_Stress_T0(2*n3-1:2*n3,:)=Mechanicalobj.Them_Stress_T0(2*n3-1:2*n3,:)+Them_Stress_T0e(5:6,:);
            end
            
        end %end of MatrixAssembly
        
        function [displacement] = MechanicalSolver(Mechanicalobj,T,U)
            
            [Piezo_Stress, Them_Stress] = Piezo_Them(Mechanicalobj,T,U);
            displacement = Mechanicalobj.K\(Them_Stress-Mechanicalobj.Them_Stress_T0+Piezo_Stress+Mechanicalobj.F);
        end %end of MechanicalSolver
        
        function [Piezo_Stress, Them_Stress] = Piezo_Them(Mechanicalobj,T,U)
            
            Piezo_Stress = sparse(2*Mechanicalobj.FEMgrid.N_n,1);
            Them_Stress = sparse(2*Mechanicalobj.FEMgrid.N_n,1);
            
            for ii_e = 1 : Mechanicalobj.FEMgrid.N_e
                
                if Mechanicalobj.FEMgrid.Element(ii_e).material == Mechanicalobj.Device.channel_material
                    elast = Mechanicalobj.Device.elast_channel;
                    themexpan = Mechanicalobj.Device.themexpan_channel;
                    piezo = Mechanicalobj.Device.piezo_channel;
                elseif Mechanicalobj.FEMgrid.Element(ii_e).material == Mechanicalobj.Device.oxide_material
                    elast = Mechanicalobj.Device.elast_oxide;
                    themexpan = Mechanicalobj.Device.themexpan_oxide;
                    piezo = Mechanicalobj.Device.piezo_oxide;
                end
                
                n1 = Mechanicalobj.FEMgrid.Element(ii_e).n1;
                n2 = Mechanicalobj.FEMgrid.Element(ii_e).n2;
                n3 = Mechanicalobj.FEMgrid.Element(ii_e).n3;
                
                Ue = [U(n1) U(n2) U(n3)].';
                Te = [T(n1) T(n2) T(n3)].';
                
                Them_Stress_e = Mechanicalobj.FEMgrid.Element(ii_e).Be_M.'*elast*themexpan*Mechanicalobj.FEMgrid.Element(ii_e).Ge*Te;
                Piezo_Stress_e = Mechanicalobj.FEMgrid.Element(ii_e).Be_M.'*piezo.'*Mechanicalobj.FEMgrid.Element(ii_e).Be*Ue*Mechanicalobj.FEMgrid.Element(ii_e).A;
                
                Them_Stress(2*n1-1:2*n1,:)=Them_Stress(2*n1-1:2*n1,:)+Them_Stress_e(1:2,:);
                Them_Stress(2*n2-1:2*n2,:)=Them_Stress(2*n2-1:2*n2,:)+Them_Stress_e(3:4,:);
                Them_Stress(2*n3-1:2*n3,:)=Them_Stress(2*n3-1:2*n3,:)+Them_Stress_e(5:6,:);
                
                Piezo_Stress(2*n1-1:2*n1,:)=Piezo_Stress(2*n1-1:2*n1,:)+Piezo_Stress_e(1:2,:);
                Piezo_Stress(2*n2-1:2*n2,:)=Piezo_Stress(2*n2-1:2*n2,:)+Piezo_Stress_e(3:4,:);
                Piezo_Stress(2*n3-1:2*n3,:)=Piezo_Stress(2*n3-1:2*n3,:)+Piezo_Stress_e(5:6,:);
                
            end
            
        end %end of Piezo_Them
%         function MatrixAssembly(Mechanicalobj,device,FEMgrid)
%             
%             N_n = length(FEMgrid.Node);
%             N_e = length(FEMgrid.Element);
%             
%             Mechanicalobj.K = sparse(2*N_n,2*N_n);
%             Mechanicalobj.R = sparse(2*N_n,1);
%             t = 1;
%             
%             for ii_e = 1 : N_e
%                 yjk = FEMgrid.Element(ii_e).y2 - FEMgrid.Element(ii_e).y3;
%                 yki = FEMgrid.Element(ii_e).y3 - FEMgrid.Element(ii_e).y1;
%                 yij = FEMgrid.Element(ii_e).y1 - FEMgrid.Element(ii_e).y2;
%                 xkj = FEMgrid.Element(ii_e).x3 - FEMgrid.Element(ii_e).x2;
%                 xik = FEMgrid.Element(ii_e).x1 - FEMgrid.Element(ii_e).x3;
%                 xji = FEMgrid.Element(ii_e).x2 - FEMgrid.Element(ii_e).x1;
%    
%                 M = [yjk yki yij;0 0 0;xkj xik xji];
%                 N = [0 0 0;xkj xik xji;yjk yki yij];
%    
%                 Kele = t/(4*Ele(ii_e).A)*[M.'*D*M M.'*D*N;N.'*D*M N.'*D*N];
%    
%                 k11 = [Kele(1,1) Kele(1,4);Kele(4,1) Kele(4,4)];
%                 k12 = [Kele(1,2) Kele(1,5);Kele(4,2) Kele(4,5)];
%                 k13 = [Kele(1,3) Kele(1,6);Kele(4,3) Kele(4,6)];
%                 k21 = [Kele(2,1) Kele(2,4);Kele(5,1) Kele(5,4)];
%                 k22 = [Kele(2,2) Kele(2,5);Kele(5,2) Kele(5,5)];
%                 k23 = [Kele(2,3) Kele(2,6);Kele(5,3) Kele(5,6)];
%                 k31 = [Kele(3,1) Kele(3,4);Kele(6,1) Kele(6,4)];
%                 k32 = [Kele(3,2) Kele(3,5);Kele(6,2) Kele(6,5)];
%                 k33 = [Kele(3,3) Kele(3,6);Kele(6,3) Kele(6,6)];
%    
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1)+k11;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2)+k12;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3)+k13;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1)+k21;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2)+k22;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3)+k23;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n1-1:2*FEMgrid.Element(ii_e).n1)+k31;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n2-1:2*FEMgrid.Element(ii_e).n2)+k32;
%                 Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3) = Mechanicalobj.K(2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3,2*FEMgrid.Element(ii_e).n3-1:2*FEMgrid.Element(ii_e).n3)+k33;
%             end
%             
%             for ii_n = 1 : N_n
%                 if(device.Node(ii_n).bm == 2)
%                     Mechanicalobj.R(2*ii_n-1) = 1;
%                 end
%             end
%             
%             for ii_n = 1 : N_n
%                 if (device.Node(ii_n).bm == 3)
%                     Mechanicalobj.K(2*ii_n-1,:) = sparse(1,2*N_n);
%                     Mechanicalobj.K(2*ii_n,:) = sparse(1,2*N_n);
%                     Mechanicalobj.K(2*ii_n-1:2*ii_n,2*ii_n-1:2*ii_n)=[1 0;0 1];
%                 end
%             end
%             
% %             d = Mechanicalobj.K \ Mechanicalobj.R;
% %             Mechanicalobj.ux = d(1:2:2*N_n);
% %             Mechanicalobj.uy = d(2:2:2*N_n);
%         end %end of MatrixAssembly
        
%         function Displacement(Mechanicalobj)
%             d = Mechanicalobj.K \ Mechanicalobj.R;
%             Mechanicalobj.ux = d(1:2:2*N_n);
%             Mechanicalobj.uy = d(2:2:2*N_n);
%         end %end of displacement
        
    end %end of methods
    
end %end of classdef

