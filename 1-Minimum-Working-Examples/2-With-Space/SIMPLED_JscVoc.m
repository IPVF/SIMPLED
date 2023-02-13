% MAIN_ComputeContributionsDecay
%
% This code aims at simulating a decay with a simple drift diffusion model
% and assigning the decay to its origins.
%
% INPUTS:
%   - parameters of the drift diffusion model
%   - duration of the simulation after the pulse
%   - laser fluence

% For the plots:
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% DEFINTION OF THE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
durationOfSimulation = 1000e-9; % in seconds
laserFluence = 5e13; % in absorbed ph/cm²

% Color Palette
colorPalette=[  [38, 70, 83]/255,
    [42, 157, 143]/255,
    [233, 196, 106]/255,
    [244, 162, 97]/255,
    [231, 111, 81]/255,
    [144, 44, 20]/255 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the parameter struct
parameters=struct();

%%  parameters
parameters.thickness = 0.5e-4;%in cm
parameters.zvec = linspace(0,parameters.thickness,200); %in cm
parameters.tvec = logspace(-14,log10(durationOfSimulation),10000);

% TRANSPORT PROPERTIES
parameters.D = 5e-3; % cm^2/s
parameters.ni = 1e4; % Carrier density

% RECOMBINATION : RADIATIVE
parameters.k2 = 1e-10; %cm^3/s

% DOPING FOR THE RADIATIVE REC
parameters.Ndop = 1e14;

%  RECOMBINATION : SRH
parameters.srhOn = 1; % 1 if on, 0 if off.
parameters.Nbulk = 0;% per cc
parameters.k1 = 1/50e-9;% bulk recombination rate in s^-1

% TOP SRH
parameters.Stop = 0;% Top surface rec. vel. in cm/s
parameters.Ntop = 0; % per cc

% Bottom SRH
parameters.Sbot = 0;% Bot. surface rec. vel. in cm/s
parameters.Nbot = 0; % per cc

% INITIAL CONDITION
parameters.alpha = 1e5; % absorption coefficient at the laser wavelength cm^-1
parameters.ngamma = laserFluence; % Fluence of the laser ph.cm^-2 ;
parameters.initialElectrons = @(x,parameters) parameters.ni; % Initial carrier density

% Temporal pulse parameters: (the pulse is temporally Gaussian with a width
% sigmaPulse and a mean time position tPulse
parameters.sigmaPulse = 5e-12;
parameters.tPulse = 50e-12;

% GENERATION RATE
parameters.factorFluence = 1;
parameters.generationRate =  @(z,t,modelParameterZ) 1*modelParameterZ.factorFluence*modelParameterZ.alpha*exp(-modelParameterZ.alpha*z)*(modelParameterZ.ngamma)*(1/(sqrt(2*pi)*modelParameterZ.sigmaPulse)).*exp(-(t-modelParameterZ.tPulse).^2/(2*modelParameterZ.sigmaPulse^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parametersContinuous = parameters;
parametersContinuous.cw_flux = 1e17; % in ph/cm²/s
parametersContinuous.tvec = logspace(-14,log10(10e-6),3000);
parametersContinuous.uniformG = false;
parametersContinuous.generationRate =   @(z,t,modelParameterZ) modelParameterZ.cw_flux*modelParameterZ.alpha*(modelParameterZ.uniformG + (1-modelParameterZ.uniformG)*exp(-modelParameterZ.alpha*z));









%% TRANSIENT SOLUTION OF THE PROBLEM
% Solving the equations
[sol] = solveSystemEquations(parameters);

%% CONTINUOUS SOLUTION OF THE PROBLEM
[sol_cw] = solveSystemEquations(parametersContinuous);

%% Delta n(t)
deltan = (1/parameters.thickness)*squeeze(trapz(parameters.zvec,sol(:,:,1),2));


%% Delta n(cw)
deltan_cw = (1/parameters.thickness)*squeeze(trapz(parameters.zvec,sol_cw(:,:,1),2));

figure()
semilogy(parametersContinuous.tvec,deltan_cw)
xlabel("Time (s)")
ylabel("$\Delta n$ [cm$^{-3}$]")
title("Reaching the cw case")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting dn/dt vs. log(IPL)
%% PL(t)
PLsignal = squeeze(trapz(parameters.zvec,sol(:,:,1).*(sol(:,:,1)+parameters.Ndop),2));
maxiPL = max(PLsignal);

dndt = diff(deltan)./diff(parameters.tvec') - squeeze(sol(1:end-1,1,1)*parameters.Stop- sol(1:end-1,end,1)*parameters.Sbot);
JscMax = 1.6e-19 * dndt * parameters.thickness;
VocMax =0.026*log(PLsignal/parameters.thickness/(parameters.ni^2));

% figure(78);
% [val,indTBegPlot] = min(abs( parameters.tvec - 200e-12));
% %semilogy(VocMax(indTBegPlot:end-1),-JscMax(indTBegPlot:end),'DisplayName','tr-with space','LineStyle','--','LineWidth',2);
% hold on;
% xlim([0 1.5])
% ylim([1e-6 1])

fluxesForCW = logspace(10,20,30);
Jsc = zeros(size(fluxesForCW));
Voc = zeros(size(fluxesForCW));
PLdots = zeros(size(fluxesForCW));

for uniformGeneration = [true]
    parametersContinuous.uniformG = uniformGeneration;
    for k = 1:length(fluxesForCW)

        parametersContinuous.cw_flux = fluxesForCW(k);
        [sol_cw] = solveSystemEquations(parametersContinuous);

        integralOfR_vs_time = (1/parameters.thickness)*squeeze(trapz(parameters.zvec,sol_cw(:,:,1)*parametersContinuous.k1 +  (sol_cw(:,:,1).*(sol_cw(:,:,1)+parametersContinuous.Ndop)-parametersContinuous.ni^2)*parametersContinuous.k2   ,2));
        Jsc(k) = 1.6e-19 * parameters.thickness * (integralOfR_vs_time(end) - squeeze(sol_cw(end,1,1)*parameters.Stop+ sol_cw(end,end,1)*parameters.Sbot)');

        % in units of PL n(top)*p(bot)
        PLdots(k) = sol_cw(end,end)*(sol_cw(end,1)+parameters.Ndop);

        Voc(k) = 0.026*log(sol_cw(end,1)*(sol_cw(end,end)+parameters.Ndop)/(parameters.ni^2));
        

    end
    if uniformGeneration
        labelString = 'cw-uniform G';
    else
        labelString = 'cw-beer-lambert';
    end
    f=figure(79);
    yyaxis right;
    l=semilogy(Voc,Jsc,'Marker','o','DisplayName',labelString,'MarkerFaceColor','auto','LineStyle','none');
    l.MarkerFaceColor = l.Color;
    xlabel("$V_{\rm oc}$ (V)")
    %kT log(I_{PL})$,
    %ylabel("-dn/dt")
    ylabel("$J_{\rm gen}$ (A.cm$^{-2}$)")
    ylim([1e-6 1])
    xlim([1 1.5])
    axis square
    set(f,'Position',[50 50 400 300])
end



% figure(58)
% yyaxis left
% [p] = polyfit(parameters.tvec(end-20:end),log(PLsignal(end-20:end)),1);
% tau_PL = -1/p(1);
% semilogy(parameters.tvec,PLsignal/maxiPL,'DisplayName',sprintf("With Transport $\\tau_{PL}$ = %0.2fns",tau_PL*1e9))
% ylim([1e-15 1])
% ylabel("$I_\textrm{PL}$ [arb.u.]")
% yyaxis right
% [p] = polyfit(parameters.tvec(end-20:end),log(deltan(end-20:end)),1);
% tau_dyn = -1/p(1);
% semilogy(parameters.tvec,deltan,'DisplayName',sprintf("$\\tau_{dyn}$ = %0.2fns",tau_dyn*1e9))
% ylim([1e-15 1]*max(deltan))
% legend()
% ylabel("$\Delta n$ [cm$^{-3}$]")
% xlabel("Time [s]")


%% Comparing the spatial distributions

deltan_of_z_end_tr = squeeze(sol(end,:,1));
deltan_of_z_end_cw = squeeze(sol_cw(end,:,1));
meanDeltan_cw = mean(deltan_of_z_end_cw);

figure(59)
semilogy(parameters.zvec,deltan_of_z_end_tr/max(deltan_of_z_end_tr),'DisplayName','tr')
hold on;
semilogy(parameters.zvec,deltan_of_z_end_cw/max(deltan_of_z_end_cw),'DisplayName','cw')
legend()
ylabel("$\Delta n$ [norm.]")
xlabel("z-position in thickness [cm]")

%% Operational time
Gofz =  squeeze(parametersContinuous.generationRate(parametersContinuous.zvec,1,parametersContinuous));
tau_op = 1e9 * deltan_of_z_end_cw ./ Gofz;
figure(590)
semilogy(parameters.zvec,tau_op)
ylabel("$\tau_{op}$ [ns]")
xlabel("z-position in thickness [cm]")

tau_op_mean = meanDeltan_cw / mean(Gofz);
disp(sprintf("With the cw spatial distribution, tau_op = <n>/<G> = %0.2f ns.",tau_op_mean*1e9))


%% Perform measurement
result = struct();
result.sol = sol;
result.PL = PLsignal;
result.PL_rel = PLsignal / maxiPL;
result.parameters = parameters;
result.tvec = parameters.tvec;
result.zvec = parameters.zvec;


elecDens = result.sol;
sizeLateral = 2;
elecDensNormTime = zeros(sizeLateral,size(elecDens,2),size(elecDens,1));

radiativeRec = zeros(1,size(elecDens,1));
defectRec =  zeros(1,size(elecDens,1));
topSurfaceTerm =  zeros(1,size(elecDens,1));
bottomSurfaceTerm =  zeros(1,size(elecDens,1));
diffusionTerm =  zeros(1,size(elecDens,1));
PL =  zeros(1,size(elecDens,1));

rateRecTopSurfacePerSurfaceArea = zeros(1,size(elecDens,1));
rateRecBottomSurfacePerSurfaceArea = zeros(1,size(elecDens,1));
rateRecBulkSRHPerSurfaceArea = zeros(1,size(elecDens,1));
rateRecBulkRadPerSurfaceArea =  zeros(1,size(elecDens,1));
% Number absolute of rec per second at surface : surface * S * deltaN
% Number absolute of rec per second in bulk : thickness * S * bulk rec


% Loop over time:
for indT=1:size(elecDens,1)
    % Loop over space
    for x=1:sizeLateral
        elecDensNormTime(x,:,indT) = elecDens(indT,:)/max(elecDens(indT,:));
    end

    % Get the excess carrier density
    deltaN = squeeze(elecDens(indT,:));

    % Compute the raditive and bulk defect terms
    radiativeRec(indT) = trapz(result.parameters.zvec,result.parameters.k2 * (deltaN .* (deltaN+parameters.Ndop)-parameters.ni^2).*(2*deltaN+parameters.Ndop));
    defectRec(indT) = trapz(result.parameters.zvec,result.parameters.k1 *(deltaN.^2 ./(deltaN + result.parameters.Nbulk) ).*(2*deltaN+parameters.Ndop));
    topSurfaceTerm(indT) = result.parameters.Stop * deltaN(1)^2 / (deltaN(1) + result.parameters.Ntop).*(2*deltaN(1)+parameters.Ndop);
    bottomSurfaceTerm(indT) = result.parameters.Sbot * deltaN(end)^2 / (deltaN(end) + result.parameters.Nbot).*(2*deltaN(end)+parameters.Ndop);

    % Compute the spatial derivative of DeltaN for the diffusion term
    derivSpat = diff(deltaN)./diff(result.parameters.zvec);
    diffusionTerm(indT) = result.parameters.D*trapz(result.parameters.zvec(1:end-1),derivSpat.^2);

    % Compute the PL intensity
    PL(indT) = trapz(result.parameters.zvec,deltaN.*(deltaN+parameters.Ndop)-parameters.ni^2);

    % Compute the number of rec directly happening at the interfaces and
    % in the bulk
    rateRecBulkRadPerSurfaceArea(indT) = trapz(result.parameters.zvec,result.parameters.k2 * deltaN.^2);
    rateRecBulkSRHPerSurfaceArea(indT) =  trapz(result.parameters.zvec,result.parameters.k1 *deltaN.^2 ./(deltaN + result.parameters.Nbulk));
    rateRecTopSurfacePerSurfaceArea(indT) = result.parameters.Stop * deltaN(1)^2 / (deltaN(1) + result.parameters.Ntop);
    rateRecBottomSurfacePerSurfaceArea(indT) = result.parameters.Sbot * deltaN(end)^2 / (deltaN(end) + result.parameters.Nbot);
end

% Compute the sum of all contributions
totalVar = -4*(diffusionTerm + topSurfaceTerm + bottomSurfaceTerm + defectRec + radiativeRec);
% Compute the derivative of the PL "independantly" from the contrib
% (explains not normalized exactly to 1 values) but allows to detect non
% detection of some recombination
dPLdt = diff(PL)./diff(result.parameters.tvec);

% Compare the two function to check in a plot: they should be very close one to
% the other.
% UNCOMMENT IF DEBUG
% figure();
% plot(result.tvec(1:end-1),dPLdt);
% hold on;
% plot(result.tvec,totalVar);
% figure();
% semilogy(result.tvec,PL);

Yrel = zeros(length(result.tvec)-1,5);
Yrel(:,1)=-2*diffusionTerm(1:end-1)./dPLdt;
Yrel(:,2)=-radiativeRec(1:end-1)./dPLdt;
Yrel(:,3)=-defectRec(1:end-1)./dPLdt;
Yrel(:,4)=-topSurfaceTerm(1:end-1)./dPLdt;
Yrel(:,5)=-bottomSurfaceTerm(1:end-1)./dPLdt;

labels = ["Diffusion" "Radiative Rec." "Bulk Rec." "Top Surface Rec." "Bot. Surface Rec."];

f=figure(10);
set(f,'Position',[500 500 400 450])
[val,indTBegPlot] = min(abs( result.tvec - 200e-12));
%indTEndPlot = indTBegPlot+findCloserIndexInList(-diff(PL(indTBegPlot:end))./PL(indTBegPlot:end-1),2e-3); % If the PL decay has ended (no relative change of more than 0.01%).
customArea(result.tvec(indTBegPlot:end-1),Yrel(indTBegPlot:end,1:5),labels,0.9,colorPalette)
ylim([0 1.05])
title("Loss Analysis")
xlabel("Time [s]")
ylabel("Proportion of $$dPL$$/$$dt$$ caused by")
set(gca,'Box','on')
lgd=legend('Location','Southoutside');
lgd.NumColumns = 2;

% Computation of the total number of recombination squared
totalDiffusion =  trapz(result.tvec(indTBegPlot:end),diffusionTerm(indTBegPlot:end));
totalTopSurface =  trapz(result.tvec(indTBegPlot:end),topSurfaceTerm(indTBegPlot:end));
totalBotSurface =  trapz(result.tvec(indTBegPlot:end),bottomSurfaceTerm(indTBegPlot:end));
totalBulk =  trapz(result.tvec(indTBegPlot:end),defectRec(indTBegPlot:end));
totalRad =  trapz(result.tvec(indTBegPlot:end),radiativeRec(indTBegPlot:end));
total = totalDiffusion+totalTopSurface+totalBulk+totalRad+totalBotSurface;
dataBar = [totalDiffusion/total, totalRad/total, totalBulk/total, totalTopSurface/total, totalBotSurface/total];

names={'Diffusion','Rad. Rec.','Bulk Rec.','Top Rec.','Bot. Rec.'};

f=figure(30);
set(f,'Position',[100 100 300 300])
for i=1:5
    X = categorical(names(i));
    b=bar(X,dataBar(i));
    set(b,'FaceColor',colorPalette(i,:));
    text(X,dataBar(i)+0.05,sprintf('%0.2f',dataBar(i)),'HorizontalAlignment','center');
    hold on;
end
ylabel("Share of total derivative for PL")
ylim([0 1])
title("Share of PL Decay")

totalRecPerSurfaceArea_bulk_rad = trapz(result.tvec,rateRecBulkRadPerSurfaceArea);
totalRecPerSurfaceArea_bulk_srh = trapz(result.tvec,rateRecBulkSRHPerSurfaceArea);
totalRecPerSurfaceArea_topSurface = trapz(result.tvec,rateRecTopSurfacePerSurfaceArea);
totalRecPerSurfaceArea_botSurface = trapz(result.tvec,rateRecBottomSurfacePerSurfaceArea);
totalRecPerSurfaceArea = totalRecPerSurfaceArea_topSurface+totalRecPerSurfaceArea_bulk_srh+totalRecPerSurfaceArea_bulk_rad+totalRecPerSurfaceArea_botSurface;
dataBar2 = [totalRecPerSurfaceArea_bulk_rad,totalRecPerSurfaceArea_bulk_srh,totalRecPerSurfaceArea_topSurface,totalRecPerSurfaceArea_botSurface]/totalRecPerSurfaceArea;
names2={'Rad. Rec.','Bulk Rec.','Top Rec.','Bot. Rec.'};

f=figure(32);
set(f,'Position',[100 100 300 330])
for i=1:4
    X2 = categorical(names2(i));
    b=bar(X2,dataBar2(i));
    set(b,'FaceColor',colorPalette(i+1,:));
    text(X2,dataBar2(i)+0.05,sprintf('%0.2f',dataBar2(i)),'HorizontalAlignment','center');
    hold on;
end
ylabel("Share of total Recombination")
title({"Share of Recombination",sprintf("After %0.2f ns",result.parameters.tvec(end)*1e9),strcat(sprintf("%0.3f",100*totalRecPerSurfaceArea/result.parameters.ngamma),"\% of carriers recombined")})
ylim([0 1])

% To export the data:
%writematrix(YrelExport','testExportCsv5.xlsx')
