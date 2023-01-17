
function [sol]=solveSystemEquations(parameters)
    % solved the partial differential drift diffusion equation
    [sol] = pdepe(0,@(z,t,n,DnDz) PDESolveDiffRad(z,t,n,DnDz,parameters),...
        @(z) PDEic(z,parameters),@(zl,ul,zr,ur,t) PDEbc(zl,ul,zr,ur,t,parameters),parameters.zvec,parameters.tvec);
    
end

function [c,f,s]= PDESolveDiffRad(x,t,u,DuDx,parameters)
c=[1];

% currents
f = current_diff(x,t,u,DuDx,parameters);

% RECOMBINATION
s = recomb_rad(x,t,u,DuDx,parameters);
s = s + recomb_srh(x,t,u,DuDx,parameters); % takes the generation into account as well
s = s + [parameters.generationRate(x,t,parameters)];

end

function u0=PDEic(x,parameters)
u0=[parameters.initialElectrons(x,parameters)];
end

function [pl,ql,pr,qr]=PDEbc(zl,ul,zr,ur,t,parameters)
% BOUNDARY CONDITIONS : NO CURRENT FLOW AT THE BORDERS

srhTermTop = -1 * parameters.Stop*(ul(1)^2-parameters.ni^2)/(ul(1) + parameters.Ntop); 
srhTermBot = 1 * parameters.Sbot*(ur(1)^2-parameters.ni^2)/(ur(1) + parameters.Nbot);

pl = [srhTermTop];%ul-ur;
ql = [1];
pr = [srhTermBot];%ul-ur;
qr = [1];
end

