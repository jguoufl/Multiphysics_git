clear all;
close all;

%read FEM grid to create Node and Element
length_scale = 1e-9; %unit in m
FEMgrid = FEMGrid('QM.n', 'QM.e', length_scale);
N_n = length(FEMgrid.Node);
N_e = length(FEMgrid.Element);

%construct the device geometry structure and initialize device parameters
device = Device(FEMgrid);
%create and initilize Poisson, DDSolver, Current, Thermal object.
poissonobj = Poisson(device, FEMgrid);
ddobj = DD(device);
%currentobj = Current(device);
thermalobj = Thermal(device,FEMgrid);
mechanicalobj = Mechanical(device,FEMgrid);

N_Vg = 1;
Vg = linspace(0.5,0.5,N_Vg);
N_Vd = 1;
Vd = linspace(0.5,0.5,N_Vd);

for ii_vg = 1 : N_Vg
   for ii_vd = 1 : N_Vd
      
       Fn_bias_old = zeros(N_n,1);
       Fp_bias_old = -Vd(ii_vd)*ones(N_n,1);
       T_old = 300*ones(N_n,1);
       dis_old = zeros(2*N_n,1);
       U_old = zeros(N_n,1);
       
       [U] = poissonobj.PoissonSolver(Fn_bias_old,Fp_bias_old,U_old,1,Vg(ii_vg), Vd(ii_vd),dis_old);
       
       dis_error = 10;
       dis_criterion = 0.01;
       while dis_error > dis_criterion
           
           tem_error = 10;
           tem_criterion = 0.01;
           while tem_error > tem_criterion
               
               Fn_error = 10;
               Fn_criterion = 0.01;
               while Fn_error > Fn_criterion
                   U_old = U;
                   [U] = poissonobj.PoissonSolver(Fn_bias_old,Fp_bias_old,U_old,0,Vg(ii_vg), Vd(ii_vd),dis_old);
                   [Ne_bias,Ph_bias,Fn_bias,Fp_bias] = ddobj.DDSolver(FEMgrid,Fn_bias_old,Fp_bias_old,U,T_old);
                   Fn_error = max(abs(Fn_bias_old(device.ind_channel)-Fn_bias(device.ind_channel)))/max(abs(Fn_bias(device.ind_channel)))
                   Fn_bias_old = Fn_bias;
                   Fp_bias_old = Fp_bias;
               end %end of electrical
               
               [T] = thermalobj.ThermalSolver(U,Ne_bias,Ph_bias,Fn_bias,Fp_bias,T_old);
               tem_error = max(abs(T-T_old))/max(abs(T))
               T_old = T;
           end %end of thermal
           
           [dis] = mechanicalobj.MechanicalSolver(T,U);
           dis_error = max(abs(dis-dis_old))/max(abs(dis))
           dis_old = dis;
       end %end of Mechancial
       
       Ne(ii_vg,ii_vd,:) = full(Ne_bias);
       Ph(ii_vg,ii_vd,:) = full(Ph_bias);
   end
end

