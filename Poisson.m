classdef Poisson < handle & CommonFunctions
    %Poisson Summary of this class goes here
    %   Detailed explanation goes here
    %Poisson is used to solve the Poisson equation.
    %The Poisson object will do the matrix assembly once it is created. The
    %member function 'PoissonSolver' is used to solve the Poisson equation
    %with approriate input arguments. It will update the boudary conditions
    %at each call.
    properties
        K %the Laplace operator matrix
        S %the overlap matrix
        B %the boundary condition vector
        Device
        FEMgrid
    end
    
    methods
        
        function Poissonobj = Poisson(device, FEMgrid)
            
            %declare deviceobj is inherited from CommonFunctions
            Poissonobj = Poissonobj@CommonFunctions();
            
            Poissonobj.Device = device;
            Poissonobj.FEMgrid = FEMgrid;
            %construct the K matrix for Poisson solver
            MatrixAssembly(Poissonobj, device, FEMgrid);
         
        end %end of constructor
        
        function MatrixAssembly(Poissonobj, device, FEMgrid)
            
            Poissonobj.K = sparse(FEMgrid.N_n,FEMgrid.N_n); 
            Poissonobj.S = sparse(FEMgrid.N_n,FEMgrid.N_n);
            
            for ii_e = 1 : FEMgrid.N_e
                
                if FEMgrid.Element(ii_e).material == device.channel_material
                    epsil = device.epsil_channel;
                elseif FEMgrid.Element(ii_e).material == device.oxide_material
                    epsil = device.epsil_oxide;
                end
                
                De = [epsil 0;0 epsil];
                Kvve = FEMgrid.Element(ii_e).Be.'*De*FEMgrid.Element(ii_e).Be*FEMgrid.Element(ii_e).A;
                
                n1 = FEMgrid.Element(ii_e).n1;
                n2 = FEMgrid.Element(ii_e).n2;
                n3 = FEMgrid.Element(ii_e).n3;
                
                Poissonobj.K(n1,n1)=Poissonobj.K(n1,n1)+Kvve(1,1);
                Poissonobj.K(n1,n2)=Poissonobj.K(n1,n2)+Kvve(1,2);
                Poissonobj.K(n1,n3)=Poissonobj.K(n1,n3)+Kvve(1,3);
                Poissonobj.K(n2,n1)=Poissonobj.K(n2,n1)+Kvve(2,1);
                Poissonobj.K(n2,n2)=Poissonobj.K(n2,n2)+Kvve(2,2);
                Poissonobj.K(n2,n3)=Poissonobj.K(n2,n3)+Kvve(2,3);
                Poissonobj.K(n3,n1)=Poissonobj.K(n3,n1)+Kvve(3,1);
                Poissonobj.K(n3,n2)=Poissonobj.K(n3,n2)+Kvve(3,2);
                Poissonobj.K(n3,n3)=Poissonobj.K(n3,n3)+Kvve(3,3);
                
                Poissonobj.S(n1,n1)=Poissonobj.S(n1,n1)+FEMgrid.Element(ii_e).Se(1,1);
                Poissonobj.S(n1,n2)=Poissonobj.S(n1,n2)+FEMgrid.Element(ii_e).Se(1,2);
                Poissonobj.S(n1,n3)=Poissonobj.S(n1,n3)+FEMgrid.Element(ii_e).Se(1,3);
                Poissonobj.S(n2,n1)=Poissonobj.S(n2,n1)+FEMgrid.Element(ii_e).Se(2,1);
                Poissonobj.S(n2,n2)=Poissonobj.S(n2,n2)+FEMgrid.Element(ii_e).Se(2,2);
                Poissonobj.S(n2,n3)=Poissonobj.S(n2,n3)+FEMgrid.Element(ii_e).Se(2,3);
                Poissonobj.S(n3,n1)=Poissonobj.S(n3,n1)+FEMgrid.Element(ii_e).Se(3,1);
                Poissonobj.S(n3,n2)=Poissonobj.S(n3,n2)+FEMgrid.Element(ii_e).Se(3,2);
                Poissonobj.S(n3,n3)=Poissonobj.S(n3,n3)+FEMgrid.Element(ii_e).Se(3,3);
                
            end
            Poissonobj.K = -Poissonobj.K;
            %modify the K, S, B according to the boundary conditions
            for ii_n = 1 : FEMgrid.N_n
                if FEMgrid.Node(ii_n).bm == device.gate_bm 
                    Poissonobj.K(ii_n,:) = sparse(1,FEMgrid.N_n); 
                    Poissonobj.S(ii_n,:) = sparse(1,FEMgrid.N_n);
                    Poissonobj.K(ii_n,ii_n) = 1;
                elseif FEMgrid.Node(ii_n).bm == device.drain_bm 
                    Poissonobj.K(ii_n,:) = sparse(1,FEMgrid.N_n);
                    Poissonobj.S(ii_n,:) = sparse(1,FEMgrid.N_n);
                    Poissonobj.K(ii_n,ii_n) = 1;
                    Poissonobj.B(ii_n,1) = device.Udrain;
                elseif FEMgrid.Node(ii_n).bm == device.source_bm 
                    Poissonobj.K(ii_n,:) = sparse(1,FEMgrid.N_n);
                    Poissonobj.S(ii_n,:) = sparse(1,FEMgrid.N_n);
                    Poissonobj.K(ii_n,ii_n) = 1;     
                    Poissonobj.B(ii_n,1) = device.Usource;
                end
            end
            
        end %end of MatrixAssembly
        
%         function MatrixAssembly(Poissonobj, device, FEMgrid)
%             
%             N_n = length(FEMgrid.Node);
%             N_e = length(FEMgrid.Element);
%             
%             Poissonobj.K = sparse(N_n,N_n); 
%             Poissonobj.S = sparse(N_n,N_n); 
%             Poissonobj.Sad = zeros(N_n,1); 
%             Poissonobj.B = sparse(N_n,1); 
%         
%             rbar = 1; %rbar=1 for the planar structure
%             
%             for ii_e = 1 : N_e
%                 
%                 if FEMgrid.Element(ii_e).material == device.channel_material
%                     epsil = device.epsil_channel;
%                 elseif FEMgrid.Element(ii_e).material == device.oxide_material
%                     epsil = device.epsil_oxide;
%                 end
%                 
%                 if FEMgrid.Element(ii_e).material ~= device.oxide_material
%                     %%%%%%%% compute the overlap matrix %%%%%%%%%%%
%                     Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)=Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/6;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)=Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)=Poissonobj.S(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)=Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)=Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/6;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)=Poissonobj.S(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/12;
%                     Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)=Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)+FEMgrid.Element(ii_e).A*rbar/12; 
%                     Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)=Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)+FEMgrid.Element(ii_e).A*rbar/12; 
%                     Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)=Poissonobj.S(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)+FEMgrid.Element(ii_e).A*rbar/6;       
%        
%                     Poissonobj.Sad(FEMgrid.Element(ii_e).n1)=Poissonobj.Sad(FEMgrid.Element(ii_e).n1)+1/4*(FEMgrid.Element(ii_e).l1*FEMgrid.Element(ii_e).d1+FEMgrid.Element(ii_e).l3*FEMgrid.Element(ii_e).d3);
%                     Poissonobj.Sad(FEMgrid.Element(ii_e).n2)=Poissonobj.Sad(FEMgrid.Element(ii_e).n2)+1/4*(FEMgrid.Element(ii_e).l1*FEMgrid.Element(ii_e).d1+FEMgrid.Element(ii_e).l2*FEMgrid.Element(ii_e).d2);
%                     Poissonobj.Sad(FEMgrid.Element(ii_e).n3)=Poissonobj.Sad(FEMgrid.Element(ii_e).n3)+1/4*(FEMgrid.Element(ii_e).l2*FEMgrid.Element(ii_e).d2+FEMgrid.Element(ii_e).l3*FEMgrid.Element(ii_e).d3);
%                 end
%                 %%%%%%%% compute the kinetic energy matrix %%%%%%%%%%%
%                 Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1) = Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n1)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x3)^2+...
%                                                (FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)^2)*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2) = Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n2)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3) = Poissonobj.K(FEMgrid.Element(ii_e).n1,FEMgrid.Element(ii_e).n3)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1) = Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n1)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2) = Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n2)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)^2+...
%                                                (FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y3)^2)*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3) = Poissonobj.K(FEMgrid.Element(ii_e).n2,FEMgrid.Element(ii_e).n3)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1) = Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n1)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x3-FEMgrid.Element(ii_e).x2)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y2-FEMgrid.Element(ii_e).y3)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2) = Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n2)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x3)*...
%                                                (FEMgrid.Element(ii_e).x2-FEMgrid.Element(ii_e).x1)+(FEMgrid.Element(ii_e).y3-FEMgrid.Element(ii_e).y1)*(FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2))*rbar;
%                 Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3) = Poissonobj.K(FEMgrid.Element(ii_e).n3,FEMgrid.Element(ii_e).n3)+(-epsil/4/FEMgrid.Element(ii_e).A)*((FEMgrid.Element(ii_e).x1-FEMgrid.Element(ii_e).x2)^2+...
%                                                (FEMgrid.Element(ii_e).y1-FEMgrid.Element(ii_e).y2)^2)*rbar;   
%             end
%             
%             %modify the K, S, Sad, B according to the boundary conditions
%             for ii_n = 1 : N_n
%                 if FEMgrid.Node(ii_n).bm == device.gate_bm 
%                     Poissonobj.K(ii_n,:) = sparse(1,N_n); 
%                     Poissonobj.S(ii_n,:) = sparse(1,N_n);
%                     Poissonobj.Sad(ii_n) = 0;
%                     Poissonobj.K(ii_n,ii_n) = 1;
%                 elseif FEMgrid.Node(ii_n).bm == device.drain_bm 
%                     Poissonobj.K(ii_n,:) = sparse(1,N_n);
%                     Poissonobj.S(ii_n,:) = sparse(1,N_n);
%                     Poissonobj.Sad(ii_n) = 0;
%                     Poissonobj.K(ii_n,ii_n) = 1;
%                     Poissonobj.B(ii_n,1) = device.Udrain;
%                 elseif FEMgrid.Node(ii_n).bm == device.source_bm 
%                     Poissonobj.K(ii_n,:) = sparse(1,N_n);
%                     Poissonobj.S(ii_n,:) = sparse(1,N_n);
%                     Poissonobj.Sad(ii_n) = 0;
%                     Poissonobj.K(ii_n,ii_n) = 1;     
%                     Poissonobj.B(ii_n,1) = device.Usource;
%                 end
%             end
%             
%             %%%%% modified area distribution matrix for discritizing the RG rate
%             Poissonobj.Sad = spdiags(Poissonobj.Sad,0,N_n,N_n);
%         end %end of MatrixAssembly
                
        function [U_bias] = PoissonSolver(Poissonobj,Fn_bias,Fp_bias,U_old,Laplace_flag,Vg,Vd,displacement)
            
            for ii_n = 1 : Poissonobj.FEMgrid.N_n
                    if Poissonobj.FEMgrid.Node(ii_n).bm == Poissonobj.Device.gate_bm
                        Poissonobj.B(ii_n, 1) = -Vg;
                    elseif Poissonobj.FEMgrid.Node(ii_n).bm == Poissonobj.Device.drain_bm
                        Poissonobj.B(ii_n, 1) = Poissonobj.Device.Udrain - Vd;
                    end
            end
            
            if Laplace_flag==1
               %%%%%% initial guess as Laplace solution
               U_bias=real(Poissonobj.K\Poissonobj.B);     % a direct method to solve AX=B
            else
                
                [Piezo_Charge]=Piezoelectric(Poissonobj,displacement);
                Np=length(Poissonobj.K);
                
                error_inner=1;
                criterion_inner=1e-3;
                while error_inner>criterion_inner
                    dummy_ro=zeros(Np,1); dummy_ro_prime=zeros(Np,1);
        
                    %%% compute the total charge density (electron as positive)
                    dummy_ro=-Poissonobj.Device.Ndav;
                    dummy_ro=dummy_ro+dummy(Poissonobj,U_old,Fn_bias,Fp_bias,Poissonobj.Device.ind_channel,0,Poissonobj.Device.Neff_e,Poissonobj.Device.Neff_h,Poissonobj.Device.Eg);
                    %%% compute the derivitive of the toatl charge density.
                    dummy_ro_prime=dummy(Poissonobj,U_old,Fn_bias,Fp_bias,Poissonobj.Device.ind_channel,1,Poissonobj.Device.Neff_e,Poissonobj.Device.Neff_h,Poissonobj.Device.Eg);

                    %%%%% Newton-Ralphson solution
                    Res=Poissonobj.K*U_old+(Poissonobj.Device.q/Poissonobj.Device.epsil_0)*Poissonobj.S*(dummy_ro)-Poissonobj.B-Piezo_Charge/Poissonobj.Device.epsil_0;
                    Jm=Poissonobj.K+(Poissonobj.Device.q/Poissonobj.Device.epsil_0)*Poissonobj.S*spdiags(dummy_ro_prime,0,Np,Np);
                    delt_U=-sparse(Jm)\Res;

                    %% bound delt_U to eliminate non-physical values
                    for i_node=1:length(delt_U)
                        if(abs(delt_U(i_node))<=1)
                            delt_U(i_node)=delt_U(i_node);
                        elseif(1<abs(delt_U(i_node)) && abs(delt_U(i_node)) <3.7)
                            delt_U(i_node)=sign(delt_U(i_node))*power(abs(delt_U(i_node)),1/5);
                        elseif(abs(delt_U(i_node))>=3.7)
                            delt_U(i_node)=sign(delt_U(i_node))*log(abs(delt_U(i_node)));
                        end
                    end
        
                    %%% update the parameters
                    error_inner=max(abs(full(delt_U)));
                    U_bias=U_old+delt_U;
                    U_old=U_bias;        
                end
            end
        end %end of PoissonSolver

        function [Piezo_Charge] = Piezoelectric(Poissonobj, displacement)
            
            Piezo_Charge = sparse(Poissonobj.FEMgrid.N_n,1);
            
            for ii_e = 1 : Poissonobj.FEMgrid.N_e
                
                if Poissonobj.FEMgrid.Element(ii_e).material == Poissonobj.Device.channel_material
                    piezo = Poissonobj.Device.piezo_channel;
                elseif Poissonobj.FEMgrid.Element(ii_e).material == Poissonobj.Device.oxide_material
                    piezo = Poissonobj.Device.piezo_oxide;
                end
                
                n1 = Poissonobj.FEMgrid.Element(ii_e).n1;
                n2 = Poissonobj.FEMgrid.Element(ii_e).n2;
                n3 = Poissonobj.FEMgrid.Element(ii_e).n3;
                de = [displacement(n1*2-1,1);
                      displacement(n1*2,1);
                      displacement(n2*2-1,1);
                      displacement(n2*2,1);
                      displacement(n3*2-1,1);
                      displacement(n3*2,1)];
                Kvde = Poissonobj.FEMgrid.Element(ii_e).Be.'*piezo*Poissonobj.FEMgrid.Element(ii_e).Be_M*de*Poissonobj.FEMgrid.Element(ii_e).A;
                
                Piezo_Charge(n1,1)=Piezo_Charge(n1,1)+Kvde(1,1);
                Piezo_Charge(n2,1)=Piezo_Charge(n2,1)+Kvde(2,1);
                Piezo_Charge(n3,1)=Piezo_Charge(n3,1)+Kvde(3,1);
                
            end
        end %end of Piezoelectric
        
    end %end of methods
    
end %end of classdef

