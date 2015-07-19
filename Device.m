classdef Device < handle & CommonFunctions
    %Device Summary of this class goes here
    %   Detailed explanation goes here
    %Device is used to specify device paramters, as well as device boundary
    %conditons that well be used for Poisson, DD, Thermal and Mechinical
    %classes.
    properties(SetAccess = private)
        %Physical constants
        kBT = 0.0259; %unit in eV
        hbar = 1.055e-34; %unit in m^2*kg/s
        epsil_0 = 8.854e-12; %unit in F/m
        q = 1.6e-19; %unit in C
        m0 = 9.11e-31; %unit in kg
        
        %Device Geometry
        LScale = 1e-9; %length scale, unit in m

        LChannel = 100; %unit in nm
        LSource = 10; %unit in nm
        LDrain = 10; %unit in nm
        TOxide = 5; %unit in nm
        TChannel = 5; %unit in nm
        
%         %Device Material
%         %InAs membrane is used as the channel material
%         Eg = 0.36; %band gap, unit in eV
%         me_channel = 0.023; %effective mass of electron 
%         mh_channel = 0.40; %effecive mass of hole
%         Neff_e %effecitve density of state of electron
%         Neff_h %effecitve density of state of hole
%         ni = 1e21; %intrinsic carrier concentration, unit in /m^3
%         mu = 4; %electron mobility, unit in m^2/(V*s)
%         muv = 0.05; %hole mobility, unit in m^2/(V*s)
%         md_factor = -0.3; %mobility degradition factor
%         epsil_oxide = 3.9; %relative dielectric constant of oxide
%         epsil_channel = 15; %relative dielectric constant of channel
%         themcond_channel = 27; %thermal conductivity of channel, unit in W/(mK)
%         themcond_oxide = 1.4; %thermal conductivity of oxide, unit in W/(mK)
%         
%         %elastic constants
%         C_channel = [8.3e10 4.5e10 0;4.5e10 8.3e10 0;0 0 4.0e10]; %unit in Pa

        %Device Material
        %GaN is used as the channel material
        %SiO2 is used as the gate dielectric
        Eg = 3.47;
        me_channel = 0.20;
        mh_channel = 0.80;
        Neff_e = 2.2343e24;
        Neff_h = 4.6246e25;
        ni
        mu = 0.1;
        muv = 0.02;
        md_factor = -1.5;
        epsil_channel = 8.9;
        epsil_oxide = 3.9;
        themcond_channel = 130;
        themcond_oxide = 1.4;
        %elastic constants
        elast_channel = [374e9 70e9 0;70e9 379e9 0;0 0 101e9]; %unit in Pa
        elast_oxide = [87.26e9 11.95e9 0;11.95e9 105.8e9 0;0 0 57.15];
        %piezoelectric matrix
        piezo_channel = [0 0 -0.3;-0.32 0.63 0]; %unit in C/m2
        piezo_oxide = [0.171 0 0;0 0 0];
        %thermal expansion
        themexpan_channel = [5.59e-6;3.17e-6;0]; %unit in K
        themexpan_oxide = [0.6e-6;0.6e-6;0];
        
        %charge density boundary conditions
        Usource = 0; %unit in V
        Udrain = 0; %unit in V
        Nss
        Ndd
        Pss
        Pdd
        
        %channel doping
        Nd_channel = 0; %n-type doping density in channel, unit in m^3
        hole_flag = 1; %0: holes are not considered; 1: holes are considered
        
        %temperature boundary condition
        tem_boundary = 300; %temperature boundary at source/drain and gate
        
        %FEM grid boundary markers
        gate_bm = 1;
        drain_bm = 4;
        source_bm = 6;
        channel_material = 2;
        oxide_material = 1;
        
        ind_material
        ind_channel
        ind_oxide
        Ndav
    end
    
    methods
        
        function Deviceobj = Device(FEMgrid)
            
            %declare deviceobj is inherited from CommonFunctions
            Deviceobj = Deviceobj@CommonFunctions();
            %initilize the device parameters, most of the parameters can be
            %initilized in the properties directly.
            Deviceobj.me_channel = Deviceobj.me_channel * Deviceobj.m0;
            Deviceobj.mh_channel = Deviceobj.mh_channel * Deviceobj.m0;
            %Deviceobj.Neff_e = 2*(2*pi*Deviceobj.me_channel*Deviceobj.kBT*Deviceobj.q/(2*pi*Deviceobj.hbar)^2)^(3/2);
            %Deviceobj.Neff_h = 2*(2*pi*Deviceobj.mh_channel*Deviceobj.kBT*Deviceobj.q/(2*pi*Deviceobj.hbar)^2)^(3/2);
            Deviceobj.ni = (Deviceobj.Neff_e*Deviceobj.Neff_h)^(0.5)*exp(-Deviceobj.Eg/(2*Deviceobj.kBT));
            N_n = length(FEMgrid.Node);
            N_e = length(FEMgrid.Element);
            Deviceobj.ind_material = sparse(N_n,1);
            for ii_e = 1 : N_e
                if FEMgrid.Element(ii_e).material == Deviceobj.channel_material
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n1) = Deviceobj.channel_material;
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n2) = Deviceobj.channel_material;
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n3) = Deviceobj.channel_material;
                elseif FEMgrid.Element(ii_e).material == Deviceobj.oxide_material
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n1) = Deviceobj.oxide_material;
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n2) = Deviceobj.oxide_material;
                    Deviceobj.ind_material(FEMgrid.Element(ii_e).n3) = Deviceobj.oxide_material;
                end
            end
            
            Deviceobj.ind_channel = find(Deviceobj.ind_material == Deviceobj.channel_material);
            Deviceobj.ind_oxide = find(Deviceobj.ind_material == Deviceobj.oxide_material);
            Deviceobj.Ndav = sparse(N_n,1);
            Deviceobj.Ndav(Deviceobj.ind_channel) = Deviceobj.Nd_channel;
            
            Deviceobj.Nss = Deviceobj.Neff_e*fermi(Deviceobj,-Deviceobj.Usource/Deviceobj.kBT,1,1/2);
            Deviceobj.Ndd = Deviceobj.Neff_e*fermi(Deviceobj,-Deviceobj.Udrain/Deviceobj.kBT,1,1/2);
            Deviceobj.Pss = Deviceobj.Neff_h*fermi(Deviceobj,(Deviceobj.Usource-Deviceobj.Eg)/Deviceobj.kBT,1,1/2);
            Deviceobj.Pdd = Deviceobj.Neff_h*fermi(Deviceobj,(Deviceobj.Udrain-Deviceobj.Eg)/Deviceobj.kBT,1,1/2);
        end %end of constructor
        
    end %end of methods
    
end %end of classdef

