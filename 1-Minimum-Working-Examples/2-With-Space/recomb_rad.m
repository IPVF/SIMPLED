
% Radiative recomb
function s=recomb_rad(x,t,u,DuDx,parameters)
s=[-(u(1)*(u(1)+parameters.Ndop)-parameters.ni^2)*parameters.k2];
end
