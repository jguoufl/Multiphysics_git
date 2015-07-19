% generate 3D figures
function [xlin, ylin, Z]=pr(n)

%x=n(:,1)/1e-9;y=n(:,2)/1e-9;z=n(:,3);% convert to nm
x = n(:,1);y = n(:,2);z = n(:,3);
nombx=501;
nomby=501; 
xlin = linspace(min(x),max(x),nombx);
ylin = linspace(min(y),max(y),nomby);
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(x,y,z,X,Y);
viewmode=2;
if viewmode==1
    surf(X,Y,Z)
    shading interp;
    axis tight;
    %view(0,90)
    colorbar
    colormap(0.9*jet+0.1*flag)
elseif viewmode == 2
    gmapja(X,Y,Z);
    h_xlabel=get(gca, 'xlabel');    h_ylabel=get(gca, 'ylabel');
    set(h_xlabel,'string','x [m]','fontsize',[16]);
    set(h_ylabel,'string','y [m]','fontsize',[16]);
else
    %Z = exp(Z);
    pcolor(Z);
    shading interp;
end
xlin=xlin*1e-9;     % convert to m
ylin=ylin*1e-9;     % convert to m


