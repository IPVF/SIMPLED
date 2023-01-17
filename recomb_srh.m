
% Trap assisted recomb
function s=recomb_srh(x,t,u,DuDx,parameters)

if parameters.srhOn 
   
    srhTerm = -parameters.k1 * (u(1)^2-parameters.ni^2)/(u(1)+parameters.Nbulk);
    s=[srhTerm];

end

end     
    