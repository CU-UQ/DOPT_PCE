%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_sol  = le_duffing_call(inputs)
x_sol = zeros(4,1);
t0 = 0;
tend = 5;
tstep = 500;
[T,x] = ode45(@(t,y)le_duffing(t,y, inputs),linspace(t0,tend,tstep),[1; 0]');

for i=1:4
    x_sol(i) = interp1(T,x(:,1),i);
end

end

function dxdt = le_duffing(~, x, inputs)
epsilon = inputs(1);
omega = inputs(2);
xi = inputs(3);

    dxdt = [x(2); -2*omega*xi*x(2) - omega^2*(x(1) - epsilon*x(1)^3)];
end



