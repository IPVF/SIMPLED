
% Radiative recomb
function s=recomb_rad(x,t,u,DuDx,parameters)
s=[-(u(1)^2-parameters.ni^2)*parameters.k2];
end
