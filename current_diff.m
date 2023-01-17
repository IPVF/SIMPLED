% Diffusion current
function f=current_diff(x,t,u,DuDx,parameters)
f=[parameters.D*DuDx(1)];
end