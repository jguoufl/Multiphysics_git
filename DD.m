classdef DD < handle & CommonFunctions
    %DD Summary of this class goes here
    %   Detailed explanation goes here
    %DD is used to solve the contiuity equation of Drift-Diffusion model.
    %When a DD object is created, it will get the information of the device
    %and do nothing else. Matrix assembly will be done when the member
    %function 'DDSolver' is called. And it will solve the continuity
    %equation with finite volume Scharfetter-Gummel scheme.
    properties
        KD
        BD
        KDv
        BDv
        Device
    end
    
    methods
        
        function DDobj = DD(device)
            DDobj = DDobj@CommonFunctions();
            DDobj.Device = device;
        end %end of constructor
        
        function [Ne,Ph,Fn_bias,Fp_bias] = DDSolver(DDobj,FEMgrid,Fn_bias_old,Fp_bias_old,U,T)
            MatrixAssembly(DDobj,DDobj.Device,FEMgrid,Fn_bias_old,Fp_bias_old,U, T);
            [Ne,Ph,Fn_bias,Fp_bias] = Solver(DDobj,U);
        end
        
        function MatrixAssembly(DDobj,device,FEMgrid,Fn_bias_old,Fp_bias_old,U,T)
            
            factor = DDobj.Device.md_factor;
            
            %%%%%%% set up the matrices elements for the DD Eqn.
            DDobj.KD=sparse(FEMgrid.N_n,FEMgrid.N_n);   % the Laplace operator matrix
            DDobj.BD=sparse(FEMgrid.N_n,1);     % the boundary condition vector
            DDobj.KDv=sparse(FEMgrid.N_n,FEMgrid.N_n);   % the Laplace operator matrix for holes 
            DDobj.BDv=sparse(FEMgrid.N_n,1);     % the boundary condition vector for holes
            
            zetan=(Fn_bias_old-U)./device.kBT;
            zetap=-(Fp_bias_old-(U-device.Eg))./device.kBT;
            
            deg_fac_c=fermi(DDobj,zetan,1,1/2)./fermi(DDobj,zetan,1,-1/2);    % treat degeneratcy and non-parabolic band structure.
            deg_fac_p=fermi(DDobj,zetap,1,1/2)./fermi(DDobj,zetap,1,-1/2);    % treat degeneratcy and non-parabolic band structure.
            
            for ii_e=1:FEMgrid.N_e
                
                n1 = FEMgrid.Element(ii_e).n1;
                n2 = FEMgrid.Element(ii_e).n2;
                n3 = FEMgrid.Element(ii_e).n3;
                
                if FEMgrid.Element(ii_e).material~= device.oxide_material   % exclude the oxide region for transport
                    %%%%%%%% compute the kinetic energy matrix %%%%%%%%%%%
                    %%% the 1st edge
                    tem = (T(n1)+T(n2))/2;
                    deg_fac=(deg_fac_c(n1)+deg_fac_c(n2))/2;
                    delt=(U(n2)-U(n1))/(device.kBT*deg_fac);
                    c1=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d1/FEMgrid.Element(ii_e).l1*Bern(DDobj,delt);    % coefficient in Schaff-Gummel, for electron flux * (-1)
                    c2=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d1/FEMgrid.Element(ii_e).l1*Bern(DDobj,-delt);    % coefficient in Schaff_gunnel 
                    deg_fac=(deg_fac_p(n1)+deg_fac_p(n2))/2;
                    deltv=(U(n2)-U(n1))/(device.kBT*deg_fac);
                    v1=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d1/FEMgrid.Element(ii_e).l1*Bern(DDobj,-deltv);    % coefficient in Schaff-Gummel, for hole flux * (-1)
                    v2=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d1/FEMgrid.Element(ii_e).l1*Bern(DDobj,deltv);    % coefficient in Schaff_gunnel 
    
                    DDobj.KD(n1,n1)=DDobj.KD(n1,n1)-c1;
                    DDobj.KD(n1,n2)=DDobj.KD(n1,n2)+c2;
                    DDobj.KD(n2,n2)=DDobj.KD(n2,n2)-c2;
                    DDobj.KD(n2,n1)=DDobj.KD(n2,n1)+c1;
    
                    DDobj.KDv(n1,n1)=DDobj.KDv(n1,n1)-v1;
                    DDobj.KDv(n1,n2)=DDobj.KDv(n1,n2)+v2;
                    DDobj.KDv(n2,n2)=DDobj.KDv(n2,n2)-v2;
                    DDobj.KDv(n2,n1)=DDobj.KDv(n2,n1)+v1;
    
                    %%% the 2nd edge
                    tem = (T(n2)+T(n3))/2;
                    deg_fac=(deg_fac_c(n2)+deg_fac_c(n3))/2;
                    delt=(U(n3)-U(n2))/(device.kBT*deg_fac);    
                    c1=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d2/FEMgrid.Element(ii_e).l2*Bern(DDobj,delt);    % coefficient in Schaff-Gummel
                    c2=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d2/FEMgrid.Element(ii_e).l2*Bern(DDobj,-delt);    % coefficient in Schaff_gunnel 
                    deg_fac=(deg_fac_p(n2)+deg_fac_p(n3))/2;
                    deltv=(U(n3)-U(n2))/(device.kBT*deg_fac); 
                    v1=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d2/FEMgrid.Element(ii_e).l2*Bern(DDobj,-deltv);    % coefficient in Schaff-Gummel
                    v2=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d2/FEMgrid.Element(ii_e).l2*Bern(DDobj,deltv);    % coefficient in Schaff_gunnel  
    
                    DDobj.KD(n2,n2)=DDobj.KD(n2,n2)-c1;
                    DDobj.KD(n2,n3)=DDobj.KD(n2,n3)+c2;
                    DDobj.KD(n3,n3)=DDobj.KD(n3,n3)-c2;
                    DDobj.KD(n3,n2)=DDobj.KD(n3,n2)+c1;
    
                    DDobj.KDv(n2,n2)=DDobj.KDv(n2,n2)-v1;
                    DDobj.KDv(n2,n3)=DDobj.KDv(n2,n3)+v2;
                    DDobj.KDv(n3,n3)=DDobj.KDv(n3,n3)-v2;
                    DDobj.KDv(n3,n2)=DDobj.KDv(n3,n2)+v1;
    
                    %%% the 3rd edge
                    tem = (T(n1)+T(n3))/2;
                    deg_fac=(deg_fac_c(n3)+deg_fac_c(n1))/2;
                    delt=(U(n1)-U(n3))/(device.kBT*deg_fac);
                    c1=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d3/FEMgrid.Element(ii_e).l3*Bern(DDobj,delt);    % coefficient in Schaff-Gummel
                    c2=device.mu*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d3/FEMgrid.Element(ii_e).l3*Bern(DDobj,-delt);    % coefficient in Schaff_gunnel 
                    deg_fac=(deg_fac_p(n3)+deg_fac_p(n1))/2;
                    deltv=(U(n1)-U(n3))/(device.kBT*deg_fac); 
                    v1=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d3/FEMgrid.Element(ii_e).l3*Bern(DDobj,-deltv);    % coefficient in Schaff-Gummel
                    v2=device.muv*(tem/300)^(factor)*(device.kBT*deg_fac)*FEMgrid.Element(ii_e).d3/FEMgrid.Element(ii_e).l3*Bern(DDobj,deltv);    % coefficient in Schaff_gunnel
    
                    DDobj.KD(n3,n3)=DDobj.KD(n3,n3)-c1;
                    DDobj.KD(n3,n1)=DDobj.KD(n3,n1)+c2;
                    DDobj.KD(n1,n1)=DDobj.KD(n1,n1)-c2;
                    DDobj.KD(n1,n3)=DDobj.KD(n1,n3)+c1; 
    
                    DDobj.KDv(n3,n3)=DDobj.KDv(n3,n3)-v1;
                    DDobj.KDv(n3,n1)=DDobj.KDv(n3,n1)+v2;
                    DDobj.KDv(n1,n1)=DDobj.KDv(n1,n1)-v2;
                    DDobj.KDv(n1,n3)=DDobj.KDv(n1,n3)+v1; 
                end 
            end

            CHAROX=1;       % the dummy charge density in oxide to avoid singularity

            %% set the boundary conditions for the gate
            for ii_n=1:FEMgrid.N_n
                if FEMgrid.Node(ii_n).bm==device.source_bm     
                    DDobj.KD(ii_n,:)=sparse(1,FEMgrid.N_n);    % correct DDobj.KD and SD matrices for the boundary nodes
                    DDobj.KD(ii_n,ii_n)=1;
                    DDobj.BD(ii_n,1)=device.Nss;  % potential of the Al contact
        
                    DDobj.KDv(ii_n,:)=sparse(1,FEMgrid.N_n);    % correct DDobj.KD and SD matrices for the boundary nodes
                    DDobj.KDv(ii_n,ii_n)=1;
                    DDobj.BDv(ii_n,1)=device.Pss;  % potential of the Al contact
                elseif FEMgrid.Node(ii_n).bm==device.drain_bm     
                    DDobj.KD(ii_n,:)=sparse(1,FEMgrid.N_n);
                    DDobj.KD(ii_n,ii_n)=1;
                    DDobj.BD(ii_n,1)=device.Ndd;
        
                    DDobj.KDv(ii_n,:)=sparse(1,FEMgrid.N_n);
                    DDobj.KDv(ii_n,ii_n)=1;
                    DDobj.BDv(ii_n,1)=device.Pdd;
                end
            end
            
            for ii_ox=1:length(device.ind_oxide)       % no charge in the oxide nodes
                ii_n=device.ind_oxide(ii_ox);
                DDobj.KD(ii_n,:)=sparse(1,FEMgrid.N_n);    % correct DDobj.KD and SD matrices for the boundary nodes
                DDobj.KD(ii_n,ii_n)=1;
                DDobj.BD(ii_n,1)=CHAROX;  % charge density in oxide is zero
    
                DDobj.KDv(ii_n,:)=sparse(1,FEMgrid.N_n);    % correct DDobj.KD and SD matrices for the boundary nodes
                DDobj.KDv(ii_n,ii_n)=1;
                DDobj.BDv(ii_n,1)=CHAROX;  % charge density in oxide is zero 
            end
        end %end of MatrixAssembly
        
        function [Ne,Ph,Fn_bias,Fp_bias] = Solver(DDobj,U)
            %%% solve for the electron density
            Ne=DDobj.KD\DDobj.BD;
            Ph=DDobj.KDv\DDobj.BDv;
            % Ne_old=KD\BD;
            % Ph_old=KDv\BDv;
            % Gnp=zeros(2*N_n,1);
               
            % %%% set up the matrices
            % np_old=[Ne_old; Ph_old];
            % Ldd=sparse([KD sparse(N_n,N_n); sparse(N_n,N_n) KDv]);
            % Bv=sparse([BD; BDv]);
            % Sdd=sparse([Sad sparse(N_n,N_n); sparse(N_n,N_n) Sad]);
 
            % [np_nb]=ddNR(Ldd,Sdd,Bv,Gnp,np_old,ni);
            % Ne=np_nb(1:N_n);
            % Ph=np_nb(N_n+1:2*N_n);

            %%% compute the quasi Fermi levels
            N_n = length(U);
            Fn_bias=zeros(N_n,1);
            Fp_bias=zeros(N_n,1);
            zetan = sparse(N_n,1);
            zetap = sparse(N_n,1);
            zetan(DDobj.Device.ind_channel)=anti_dummy(DDobj,Ne(DDobj.Device.ind_channel)./DDobj.Device.Neff_e,1/2,1); 
            Fn_bias(DDobj.Device.ind_channel)=(DDobj.Device.kBT*zetan(DDobj.Device.ind_channel)+U(DDobj.Device.ind_channel));   % The electron quasi-Fermi level
    
            zetap(DDobj.Device.ind_channel)=anti_dummy(DDobj,Ph(DDobj.Device.ind_channel)./DDobj.Device.Neff_h,1/2,1); 
            Fp_bias(DDobj.Device.ind_channel)=(U(DDobj.Device.ind_channel)-DDobj.Device.Eg)-DDobj.Device.kBT*zetap(DDobj.Device.ind_channel);    % The hole quasi-Fermi level
        end %end of Solver
            
    end %end of methods
    
end %end of classdef

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           