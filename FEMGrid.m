classdef FEMGrid < handle
    %FEMGrid Summary of this class goes here
    %   Detailed explanation goes here
    %FEMGrid is used to read FEM grid information from input files to
    %construct the Node and Element. Node stores the position and boundary
    %marker of the nodes in the grid. And Elements stores the indices of
    %nodes construting an element, the nodes' positions, the element area,
    %edge lengths and so on.
    properties(SetAccess = private)
        Node
        Element
        N_n
        N_e
    end
    
    methods
        function FEMobj = FEMGrid(NodeFile, ElementFile, LengthScale)
            Nmtr = load(NodeFile);
            Emtr = load(ElementFile);
            
            FEMobj.N_n = length(Nmtr);
            for ii_n = 1 : FEMobj.N_n
               FEMobj.Node(ii_n).x = LengthScale * Nmtr(ii_n,2);
               FEMobj.Node(ii_n).y = LengthScale * Nmtr(ii_n,3);
               FEMobj.Node(ii_n).bm = Nmtr(ii_n,4);
            end
            
            FEMobj.N_e = length(Emtr);
            for ii_e = 1 : FEMobj.N_e
                FEMobj.Element(ii_e).n1 = Emtr(ii_e,2) + 1;
                FEMobj.Element(ii_e).n2 = Emtr(ii_e,3) + 1;
                FEMobj.Element(ii_e).n3 = Emtr(ii_e,4) + 1;
                FEMobj.Element(ii_e).material = Emtr(ii_e,13);
                
                FEMobj.Element(ii_e).x1 = FEMobj.Node(FEMobj.Element(ii_e).n1).x; 
                FEMobj.Element(ii_e).y1 = FEMobj.Node(FEMobj.Element(ii_e).n1).y; 
                FEMobj.Element(ii_e).x2 = FEMobj.Node(FEMobj.Element(ii_e).n2).x;    
                FEMobj.Element(ii_e).y2 = FEMobj.Node(FEMobj.Element(ii_e).n2).y;    
                FEMobj.Element(ii_e).x3 = FEMobj.Node(FEMobj.Element(ii_e).n3).x;    
                FEMobj.Element(ii_e).y3 = FEMobj.Node(FEMobj.Element(ii_e).n3).y;    
                FEMobj.Element(ii_e).A = 0.5*det([1 FEMobj.Element(ii_e).x1 FEMobj.Element(ii_e).y1; 1 FEMobj.Element(ii_e).x2 FEMobj.Element(ii_e).y2; 1 FEMobj.Element(ii_e).x3 FEMobj.Element(ii_e).y3]); 
                FEMobj.Element(ii_e).l1 = sqrt((FEMobj.Element(ii_e).x1-FEMobj.Element(ii_e).x2)^2+(FEMobj.Element(ii_e).y1-FEMobj.Element(ii_e).y2)^2); 
                FEMobj.Element(ii_e).l2 = sqrt((FEMobj.Element(ii_e).x2-FEMobj.Element(ii_e).x3)^2+(FEMobj.Element(ii_e).y2-FEMobj.Element(ii_e).y3)^2);
                FEMobj.Element(ii_e).l3 = sqrt((FEMobj.Element(ii_e).x3-FEMobj.Element(ii_e).x1)^2+(FEMobj.Element(ii_e).y3-FEMobj.Element(ii_e).y1)^2);
    
                xctmp = LengthScale * Emtr(ii_e,11);    yctmp = LengthScale * Emtr(ii_e,12);  
                FEMobj.Element(ii_e).d1 = sqrt((xctmp-(FEMobj.Element(ii_e).x1+FEMobj.Element(ii_e).x2)/2)^2+(yctmp-(FEMobj.Element(ii_e).y1+FEMobj.Element(ii_e).y2)/2)^2); 
                FEMobj.Element(ii_e).d2 = sqrt((xctmp-(FEMobj.Element(ii_e).x2+FEMobj.Element(ii_e).x3)/2)^2+(yctmp-(FEMobj.Element(ii_e).y2+FEMobj.Element(ii_e).y3)/2)^2); 
                FEMobj.Element(ii_e).d3 = sqrt((xctmp-(FEMobj.Element(ii_e).x3+FEMobj.Element(ii_e).x1)/2)^2+(yctmp-(FEMobj.Element(ii_e).y3+FEMobj.Element(ii_e).y1)/2)^2); 
                
                %basis function in the element Ni=ai+bix+ciy
                a1 = (FEMobj.Element(ii_e).x2*FEMobj.Element(ii_e).y3 - FEMobj.Element(ii_e).x3*FEMobj.Element(ii_e).y2)/(2*FEMobj.Element(ii_e).A);
                b1 = (FEMobj.Element(ii_e).y2 - FEMobj.Element(ii_e).y3)/(2*FEMobj.Element(ii_e).A);
                c1 = (FEMobj.Element(ii_e).x3 - FEMobj.Element(ii_e).x2)/(2*FEMobj.Element(ii_e).A);
                a2 = (FEMobj.Element(ii_e).x3*FEMobj.Element(ii_e).y1 - FEMobj.Element(ii_e).x1*FEMobj.Element(ii_e).y3)/(2*FEMobj.Element(ii_e).A);
                b2 = (FEMobj.Element(ii_e).y3 - FEMobj.Element(ii_e).y1)/(2*FEMobj.Element(ii_e).A);
                c2 = (FEMobj.Element(ii_e).x1 - FEMobj.Element(ii_e).x3)/(2*FEMobj.Element(ii_e).A);
                a3 = (FEMobj.Element(ii_e).x1*FEMobj.Element(ii_e).y2 - FEMobj.Element(ii_e).x2*FEMobj.Element(ii_e).y1)/(2*FEMobj.Element(ii_e).A);
                b3 = (FEMobj.Element(ii_e).y1 - FEMobj.Element(ii_e).y2)/(2*FEMobj.Element(ii_e).A);
                c3 = (FEMobj.Element(ii_e).x2 - FEMobj.Element(ii_e).x1)/(2*FEMobj.Element(ii_e).A);
                FEMobj.Element(ii_e).Be = [b1 b2 b3;c1 c2 c3];%first derivative on basis functions for Poisson, thermal, carrier density.
                %first derivative on basis functions for mechanical part to compute strain and stress.
                FEMobj.Element(ii_e).Be_M = [b1 0 b2 0 b3 0;0 c1 0 c2 0 c3;c1 b1 c2 b2 c3 b3];
           
                %intergral of N'N over the element
                FEMobj.Element(ii_e).Se = [2 1 1;1 2 1;1 1 2]*FEMobj.Element(ii_e).A/12;
                %intergral of N over the element
                FEMobj.Element(ii_e).Ge = [1 1 1]*FEMobj.Element(ii_e).A/3;
            end
        end %end of constructor
        
    end %end of methods
    
end %end of classdef

