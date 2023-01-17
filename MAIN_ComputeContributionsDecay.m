% MAIN_ComputeContributionsDecay
%
% This code aims at simulating a decay with a simple drift diffusion model
% and assigning the decay to its origins. 
%
% INPUTS:
%   - parameters of the drift diffusion model
%   - duration of the simulation after the pulse
%   - laser fluence


% DEFINTION OF THE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
durationOfSimulation = 2000e-9;
laserFluence = 1e12; % in ph/cm²

% Color Palette
colorPalette=[  [38, 70, 83]/255,	          	
                [42, 157, 143]/255,	          	
                [233, 196, 106]/255,	          	
                [244, 162, 97]/255,	          	
                [231, 111, 81]/255,	          
                [144, 44, 20]/255       ];          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the parameter struct
parameters=struct();

%%  parameters
parameters.thickness = 0.5e-4;%in cm
parameters.zvec = linspace(0,parameters.thickness,100); %in cm
parameters.tvec = logspace(-14,log10(durationOfSimulation),3000);

% TRANSPORT PROPERTIES
parameters.D = 4.5e-3; % cm^2/s
parameters.ni = 1e4;

% RECOMBINATION : RADIATIVE
parameters.k2 = 2.9e-11; %cm^3/s

%  RECOMBINATION : SRH
parameters.srhOn = 1;
parameters.Nbulk = 0;% per cc
parameters.k1 = 3.9e5;% in s^-1

% TOP SRH
parameters.Stop = 50;%cm/s
parameters.Ntop = 0; % per cc

% Bottom SRH
parameters.Sbot = 25;%cm/s
parameters.Nbot = 0; % per cc

% INITIAL CONDITION
parameters.alpha = 2*6.5e4; %cm^-1
parameters.ngamma = laserFluence; %ph.cm^-2 ; corresponds to ngamma * alpha(Laser wavelength)
parameters.initialElectrons = @(x,parameters) parameters.ni;

parameters.sigmaPulse = 5e-12;
parameters.tPulse = 50e-12;

% GENERATION RATE
parameters.factorFluence = 1;
parameters.generationRate =  @(z,t,modelParameterZ) 1*modelParameterZ.factorFluence*modelParameterZ.alpha*exp(-modelParameterZ.alpha*z)*(modelParameterZ.ngamma)*(1/(sqrt(2*pi)*modelParameterZ.sigmaPulse)).*exp(-(t-modelParameterZ.tPulse).^2/(2*modelParameterZ.sigmaPulse^2));%modelParameterZ.alpha*exp(-modelParameterZ.alpha*z)

% Post-Treatment
parameters.postTreatment = "gating";
parameters.delta = 1.3411e-9;% delay between laser and camera 0.586e-9 for reguar and 1.3411e-9 for superResolved
parameters.gateWidth = 3e-9; %s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOLUTION OF THE PROBLEM
% Solving the equations
[sol] = solveSystemEquations(parameters);

%% PL(t)
PLsignal = squeeze(trapz(parameters.zvec,sol(:,:,1).^2,2));
maxiPL = max(PLsignal);

figure(58)
semilogy(parameters.tvec,PLsignal/maxiPL)
hold on;
xlabel("Time [s]")
ylabel("PL [a.u.]")

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
    radiativeRec(indT) = trapz(result.parameters.zvec,result.parameters.k2 * deltaN.^3);
    defectRec(indT) = trapz(result.parameters.zvec,result.parameters.k1 *deltaN.^3 ./(deltaN + result.parameters.Nbulk));
    topSurfaceTerm(indT) = result.parameters.Stop * deltaN(1)^3 / (deltaN(1) + result.parameters.Ntop);
    bottomSurfaceTerm(indT) = result.parameters.Sbot * deltaN(end)^3 / (deltaN(end) + result.parameters.Nbot);
    
    % Compute the spatial derivative of DeltaN for the diffusion term
    derivSpat = diff(deltaN)./diff(result.parameters.zvec);
    diffusionTerm(indT) = result.parameters.D*trapz(result.parameters.zvec(1:end-1),derivSpat.^2);
    
    % Compute the PL intensity
    PL(indT) = trapz(result.parameters.zvec,deltaN.^2);
    
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
Yrel(:,2)=-2*radiativeRec(1:end-1)./dPLdt;
Yrel(:,3)=-2*defectRec(1:end-1)./dPLdt;
Yrel(:,4)=-2*topSurfaceTerm(1:end-1)./dPLdt;
Yrel(:,5)=-2*bottomSurfaceTerm(1:end-1)./dPLdt;

labels = ["Diffusion" "Radiative Rec." "Bulk Rec." "Top Surface Rec." "Bot. Surface Rec."];

f=figure(10);
set(f,'Position',[500 500 400 450])
[val,indTBegPlot] = min(abs( result.tvec - 200e-12));
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


%writematrix(YrelExport','testExportCsv5.xlsx')
