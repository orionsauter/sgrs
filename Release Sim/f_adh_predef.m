function [ X1, X2, X3 ] = f_adh_predef( Ad_set )

% Outputs adhesive force parameters based on the 3 cases of the paper "Prediction of the LISA-Pathfinder release mechanism in-flight performance"
% Ad_set = 0 , 1 , 2 , 3

if Ad_set==1 % SET 1 (red)
    X1=0.4e6; 
    X2=1.5191e+21;
    X3=3.4404;
elseif Ad_set==2 % SET 2 (green)
    X1=0.4e6; 
    X2=6.8008e+12;
    X3=2.1107;
elseif Ad_set==3 % SET 3 (blue)
    X1=0.4e6; 
    X2=1.0866e+07;
    X3=1.1400;
elseif Ad_set==0 % No adhesion
    X1=0;
    X2=0;
    X3=0;
else
    disp('No adhesion force profile selected')
end

end

