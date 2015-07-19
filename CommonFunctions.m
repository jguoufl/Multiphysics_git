classdef CommonFunctions < handle
    %CommonFunctions Summary of this class goes here
    %   Detailed explanation goes here
    %CommonFunctions defines a set of common functions, e.g. fermi, dummy
    %and so on. If objects of other classes want to use the functions
    %defined here, those classes must inherit from this class, and the
    %objects need to be declared as objects of this class as well in the
    %constructor.
    
    properties 
    end
    
    methods
        
        function [y]=fermi(CommonFunctionsobj,x,fermi_flag,fermi_order)
            if fermi_order==1/2
                if fermi_flag==1
                    exp_fac=exp(-0.17*(x+1.0).^2);
                    nu=x.^4+50.0+33.6*x.*(1.0-0.68*exp_fac);
                    zeta=3.0*sqrt(pi)./(4.0*nu.^0.375);
                    y=exp(x)./(1.0+zeta.*exp(x));	
                elseif fermi_flag==0
                    y=exp(x);
                end
            elseif fermi_order==0
                if fermi_flag==1
                    y=log(1+exp(x));
                elseif fermi_flag==0
                    y=exp(x);
                end
            elseif fermi_order==-1/2
                if fermi_flag==1
                    exp_fac=exp(-0.17*(x+1.0).^2);
                    nu=x.^4+50.0+33.6*x.*(1.0-0.68*exp_fac);
                    zeta=3.0*sqrt(pi)./(4.0*nu.^0.375);
                    nu_prime=4*x.^3+33.6-22.848*exp_fac.*(1-0.34*(x+x.^2));
                    zeta_prime=-(9*sqrt(pi)/32)*nu.^(-11/8).*nu_prime;
                    y=(exp(-x)-zeta_prime)./(exp(-x)+zeta).^2;
                elseif fermi_flag==0
                    y=exp(x);
                end
            end
        end %end of fermi
        
        function [y]=dummy(CommonFunctionsobj,U,Fn_bias,Fp_bias,ind_q,derivitive_flag, Neff_e, Neff_h, Eg)
           
            kBT = 0.0259;
            Np=length(U);   % total number of grid points
            y=zeros(Np,1);

            zetan=(Fn_bias-U)./kBT; % Normalized relative Fermi level position
            zetap=(Fp_bias-(U-Eg))./kBT; % Normalized relative Fermi level position


            if derivitive_flag~=1   % charge  
                %     y(ind_wire)=Nc_wire*fermi(zetan(ind_wire),1,1/2)-...
                %         Nc_wire*fermi(-zetap(ind_wire),1,1/2);     % wire part
                %     y(ind_fill)=Nc_fill*fermi(zetan(ind_fill),1,1/2)-...
                %         Nc_fill*fermi(-zetap(ind_fill),1,1/2);     % fill material part  
                y(ind_q) = Neff_e*fermi(CommonFunctionsobj,zetan(ind_q),1,1/2) - Neff_h*fermi(CommonFunctionsobj,-zetap(ind_q),1,1/2);
            else    % derivitive of charge
                %     y(ind_wire)=(-Nc_wire/kBT)*fermi(zetan(ind_wire),1,-1/2)+...
                %         (-Nc_wire/kBT)*fermi(-zetap(ind_wire),1,-1/2);
                %     y(ind_fill)=(-Nc_fill/kBT)*fermi(zetan(ind_fill),1,-1/2)+...
                %         (-Nc_fill/kBT)*fermi(-zetap(ind_fill),1,-1/2);
                y(ind_q) = (-Neff_e/kBT)*fermi(CommonFunctionsobj,zetan(ind_q),1,-1/2) + (-Neff_h/kBT)*fermi(CommonFunctionsobj,-zetap(ind_q),1,-1/2);
            end
        end %end of dummy
        
        function [y]=anti_dummy(CommonFunctionsobj,x,dummy_flag,fermi_flag)
            if dummy_flag==0
                if fermi_flag==0
                    y=log(x);
                elseif fermi_flag==1
                    y=log(exp(x)-1);
                end
            elseif dummy_flag==1/2
                if fermi_flag==0
                    y=log(x);
                elseif fermi_flag==1
                    y=log(x)+3.53553e-1*x-4.95009e-3*x.^2+1.48386e-4*x.^3-4.42563e-6*x.^4;
                end
            end
        end %end of anti_dummy
        
        function [y]=Bern(CommonFunctionsobj,x) % evaluate Bernoull'function
            y=zeros(length(x),1);
            for ii=1:length(x)
                if abs(x(ii))>1e-5
                    y(ii)=x(ii)/(exp(x(ii))-1);
                else
                    y(ii)=1;
                end
            end 
        end %end of Bern
        
        function [xlin, ylin, Z]=uniform_interpolation(CommonFunctionsobj,n,lengthscale)
            % convert to unit in lengthscale
            x=n(:,1)/lengthscale;
            y=n(:,2)/lengthscale;
            
            z=n(:,3);   
            
            nombx=501;  
            nomby=501;  
            xlin = linspace(min(x),max(x),nombx);
            ylin = linspace(min(y),max(y),nomby);
            [X,Y] = meshgrid(xlin,ylin);
            Z = griddata(x,y,z,X,Y);
            
            xlin=xlin*lengthscale;     % convert to m
            ylin=ylin*lengthscale;     % convert to m
        end
    end %end of methods
    
end %end of classdef

