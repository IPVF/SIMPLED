% Model of recombination with no space dimension
tau = 50e-9; %s
k2 = 1e-10;% cm^3/s
Ndop = 1e14; % cm^-3
pseudo_thickness = 0.5e-4; %cm
ni = 1e4; % cm^-3

tspan = [0 1e-6]; % Time span in s 
n0 = [1e18]; % Initial density of carrier

 opts = odeset('Refine',5);
[tvec,deltan] = ode15s(@(t,y) odefcn(t,y,tau,k2,Ndop), tspan, n0,opts);

PLsignal = deltan.*(deltan+Ndop);
maxiPL = max(PLsignal);

%% PLOTTING THE DECAYS
% this figure contains the n and PL decays
f=figure(58)
yyaxis left
[p] = polyfit(tvec(end-20:end),log(PLsignal(end-20:end)),1);
tau_PL = -1/p(1);
semilogy(tvec,PLsignal/maxiPL,'DisplayName',sprintf("$\\tau_{PL}$ = %0.2fns",tau_PL*1e9))
hold on;
ylim([1e-15 10])
ylabel("$I_\textrm{PL}$ [arb.u.]")
yyaxis right
[p] = polyfit(tvec(end-20:end),log(deltan(end-20:end)),1);
tau_dyn = -1/p(1);
semilogy(tvec,deltan,'DisplayName',sprintf("$\\tau_{dyn}$ = %0.2fns",tau_dyn*1e9))
hold on;
ylim([1e-15 10]*max(deltan))
set(gca,'YScale','log')
legend()
ylabel("$\Delta n$ [cm$^{-3}$]")
xlabel("Time [s]")
axis square
xlim([-0.1 1]*1e-6)
set(f,'Position',[50 50 400 300])

%% 
dndt = diff(deltan)./diff(tvec);
Jsc_max = -dndt;
IPL = log(pseudo_thickness*k2*deltan.*(Ndop+deltan));%0.026 *  /ni^2
Voc_max = 0.026*log(deltan.*(Ndop+deltan)/ni^2);
[indSwitchPL] = findCloserIndexInList(deltan,Ndop);

f=figure(78);
set(gca,'Box','on')
semilogy(IPL(1:end-1),Jsc_max,'DisplayName','tr-no space');
hold on;
xlabel("$\textrm{log}\ I_{\rm PL}$ (ph.cm$^{-2}$.s$^{-1}$)")
ylabel("$-dn/dt$ (cm$^{-3}$.s$^{-1}$)")
ylim([1e-6 1]/(1.6e-19 * pseudo_thickness) )
xlim(IPL(findCloserIndexInList(Voc_max,[1 1.5])))
set(f,'Position',[50 50 300 300])
axis square;

%% COMPUTING THE CONTRIBUTIONS of dPL/dt
% The derivative of n being defined we can compute
% Y = dPL/dt normalized

% Color Palette
colorPalette=[  [38, 70, 83]/255,	
                [42, 157, 143]/255,	          	
                [233, 196, 106]/255,	          	
                [244, 162, 97]/255,	          	
                [231, 111, 81]/255,	          
                [144, 44, 20]/255 ]; 

dPLdt = diff(PLsignal)./diff(tvec);

Yrel = zeros(length(tvec)-1,5);
Yrel(:,2)=-(k2*(2*deltan(1:end-1)+Ndop).*deltan(1:end-1).*(deltan(1:end-1)+Ndop))./dPLdt;
Yrel(:,3)=-(1*deltan(1:end-1).*(2*deltan(1:end-1)+Ndop)/tau)./dPLdt;

labels = ["Diffusion" "Radiative Rec." "Bulk Rec." "Top Surface Rec." "Bot. Surface Rec."];

f=figure(10);
set(f,'Position',[500 500 400 450])
customArea(tvec(1:end-1)',Yrel,labels,0.9,colorPalette)
ylim([0 1.05])
title("Loss Analysis dPL/dt")
xlabel("Time [s]")
ylabel("Proportion of $$dPL$$/$$dt$$ caused by")
set(gca,'Box','on')
lgd=legend('Location','Southoutside');
lgd.NumColumns = 2;

%% COMPUTING THE CONTRIBUTIONS of dPL/dt
% The derivative of n being defined we can compute
% Y =  dn/dt normalized

% Color Palette
colorPalette=[  
                [42, 157, 143]/255,	          	
                [233, 196, 106]/255,	          	
                [244, 162, 97]/255,	          	
                [231, 111, 81]/255,	          
                [144, 44, 20]/255 ]; 

dndt = diff(deltan)./diff(tvec);

Yrel = zeros(length(tvec)-1,4);
Yrel(:,1)=-(k2*deltan(1:end-1).*(deltan(1:end-1)+Ndop))./dndt;
Yrel(:,2)=-(deltan(1:end-1)/tau)./dndt;

labels = ["Radiative Rec." "Bulk Rec." "Top Surface Rec." "Bot. Surface Rec."];

f=figure(11);
set(f,'Position',[500 500 400 450])
customArea(tvec(1:end-1)',Yrel,labels,0.9,colorPalette)
ylim([0 1.05])
title("Loss Analysis dn/dt")
xlabel("Time [s]")
ylabel("Proportion of $$dPL$$/$$dt$$ caused by")
set(gca,'Box','on')
lgd=legend('Location','Southoutside');
lgd.NumColumns = 2;


function dndt = odefcn(t,n,tau,k2,Ndopant)
dndt = - n/tau - k2*n*(n+Ndopant);
end
