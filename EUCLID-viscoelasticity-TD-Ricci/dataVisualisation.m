clc;
clear all;
close all;

%% import data
load Id_TD_q1_EL22_2ME_dt0.500_Tsim300.000.mat;

%% variables
time = data.time;
plotTimeVisualisation = [3 5 length(time)];

deformationTime = data.deformation_time;                                                                  % 2*nodes x timesteps (before X rows, than Y rows)

sigmaCoordinate = 2;     % sigma11 = 1 , sigma22 = 2 , sigma33 = 3 , sigma12 = 4
sigmaCoordinateText = ["11" "22" "33" "12"];

elementStress = data.stress_time;                                                                         % nSigma x elements x timesteps 

% for CST I'll do a mean for values for each node
nodalStress = zeros(length(elementStress(:,1,1)), length(coord(:,1)), length(plotTimeVisualisation));     % nSigma x nNodes x plottedTimesteps
nodalContribution = zeros(length(coord(:,1)),1);                                                          % nNodes
for t = 1:length(plotTimeVisualisation)
    for i = 1:length(elementStress(:,1,1))
        for j = 1:length(conne(:,1))
            for k = 1:length(conne(1,:))
                nodalStress(i,conne(j,k),t) = nodalStress(i,conne(j,k),t) + elementStress(i,j,plotTimeVisualisation(t));    
                nodalContribution(conne(j,k)) = nodalContribution(conne(j,k)) + 1; 
            end
        end
    end
end
nodalContribution = nodalContribution / (length(elementStress(:,1,1))*length(plotTimeVisualisation));
nodalStress1 = nodalStress./nodalContribution';

%% mesh visualisation
nPlot = 1;
figure(nPlot);
hold on;
axis equal;
xlim([min(coord(:,2)) - 1 , max(coord(:,2)) + 1]);
ylim([min(coord(:,3)) - 1 , max(coord(:,3)) + 1]);
xlabel('x [dm]');
ylabel('y [dm]');
for i = 1:length(conne(:,1))
    plot([coord(conne(i,1),2) coord(conne(i,2),2)],[coord(conne(i,1),3) coord(conne(i,2),3)], 'k');
    plot([coord(conne(i,2),2) coord(conne(i,3),2)],[coord(conne(i,2),3) coord(conne(i,3),3)], 'k');
    plot([coord(conne(i,3),2) coord(conne(i,1),2)],[coord(conne(i,3),3) coord(conne(i,1),3)], 'k');
end
plot(coord(:,2),coord(:,3), "r.");
title("Mesh visualisation");
nPlot = nPlot +1;

%% nodal deformation
zoom = 1000;
for t = 1:length(plotTimeVisualisation)
    figure(nPlot);
    hold on;
    axis equal;
    xlim([min(coord(:,2)) - 1 , max(coord(:,2)) + 1]);
    ylim([min(coord(:,3)) - 1 , max(coord(:,3)) + 1]);
    xlabel('x [dm]');
    ylabel('y [dm]');
    % initial
    for i = 1:length(conne(:,1))
        plot([coord(conne(i,1),2) coord(conne(i,2),2)],[coord(conne(i,1),3) coord(conne(i,2),3)], ':k','LineWidth',0.5);
        plot([coord(conne(i,2),2) coord(conne(i,3),2)],[coord(conne(i,2),3) coord(conne(i,3),3)], ':k','LineWidth',0.5);
        plot([coord(conne(i,3),2) coord(conne(i,1),2)],[coord(conne(i,3),3) coord(conne(i,1),3)], ':k','LineWidth',0.5);
    end
    plot(coord(:,2),coord(:,3), "b.");
    % deformated
    for i = 1:length(conne(:,1))
        plot([coord(conne(i,1),2)+zoom*deformationTime(conne(i,1),plotTimeVisualisation(t)) coord(conne(i,2),2)+zoom*deformationTime(conne(i,2),plotTimeVisualisation(t))],[coord(conne(i,1),3)+zoom*deformationTime(conne(i,1)+length(deformationTime(:,1))/2,plotTimeVisualisation(t)) coord(conne(i,2),3)+zoom*deformationTime(conne(i,2)+length(deformationTime(:,1))/2,plotTimeVisualisation(t))], 'k','LineWidth',1.2);
        plot([coord(conne(i,2),2)+zoom*deformationTime(conne(i,2),plotTimeVisualisation(t)) coord(conne(i,3),2)+zoom*deformationTime(conne(i,3),plotTimeVisualisation(t))],[coord(conne(i,2),3)+zoom*deformationTime(conne(i,2)+length(deformationTime(:,1))/2,plotTimeVisualisation(t)) coord(conne(i,3),3)+zoom*deformationTime(conne(i,3)+length(deformationTime(:,1))/2,plotTimeVisualisation(t))], 'k','LineWidth',1.2);
        plot([coord(conne(i,3),2)+zoom*deformationTime(conne(i,3),plotTimeVisualisation(t)) coord(conne(i,1),2)+zoom*deformationTime(conne(i,1),plotTimeVisualisation(t))],[coord(conne(i,3),3)+zoom*deformationTime(conne(i,3)+length(deformationTime(:,1))/2,plotTimeVisualisation(t)) coord(conne(i,1),3)+zoom*deformationTime(conne(i,1)+length(deformationTime(:,1))/2,plotTimeVisualisation(t))], 'k','LineWidth',1.2);
    end
    plot(coord(:,2) + zoom*deformationTime(1:length(deformationTime(:,1))/2,plotTimeVisualisation(t)), ...
        coord(:,3) + zoom*deformationTime(1+length(deformationTime(:,1))/2:end,plotTimeVisualisation(t)), "r.",'MarkerSize',10);
    title({"Nodal deformation"},{"timestep " + num2str(plotTimeVisualisation(t)) + " --> t = " + num2str(time(plotTimeVisualisation(t))) + "s"});
    nPlot = nPlot +1;
end

%% nodal stress
% the constant strain triangle CST has constant strain so constant stress
% inside the element. It is possible to plot the stress value for each
% element in a single plot containing the whole mesh
lowLim = min(min(nodalStress1(sigmaCoordinate,:,:)))-10;
highLim = max(max(nodalStress1(sigmaCoordinate,:,:)))+10;

for t = 1:length(plotTimeVisualisation)
    figure(nPlot);
    hold on;
    axis equal;
    xlim([min(coord(:,2)) - 1 , max(coord(:,2)) + 1]);
    ylim([min(coord(:,3)) - 1 , max(coord(:,3)) + 1]);
    p=patch('faces',conne,'vertices',[coord(:,2) coord(:,3)],'facevertexcdata',nodalStress1(sigmaCoordinate,:,t)','FaceColor','interp');
    set(p,'EdgeColor','flat');
    colorbar;
    clim([lowLim highLim]);
    %colormap autumn
    title({"Nodal \sigma" + sigmaCoordinateText(sigmaCoordinate)},{"timestep " + num2str(plotTimeVisualisation(t)) + " --> t = " + num2str(time(plotTimeVisualisation(t))) + "s"});
    xlabel('x [dm]');
    ylabel('y [dm]');
    nPlot = nPlot +1;
end

%% element stress
lowLim = min(min(elementStress(sigmaCoordinate,:,plotTimeVisualisation)))-10;
highLim = max(max(elementStress(sigmaCoordinate,:,plotTimeVisualisation)))+10;
for t = 1:length(plotTimeVisualisation)
    figure(nPlot)
    hold on;
    axis equal;
    xlim([min(coord(:,2)) - 1 , max(coord(:,2)) + 1]);
    ylim([min(coord(:,3)) - 1 , max(coord(:,3)) + 1]);

    X = zeros(2,2);
    Y = zeros(2,2);
    Z = zeros(2,2);

    for j = 1:length(conne(:,1))
        X(1,1) = coord(conne(j,1),2);
        X(1,2) = coord(conne(j,2),2);
        X(2,1) = coord(conne(j,1),2);
        X(2,2) = coord(conne(j,3),2);

        Y(1,1) = coord(conne(j,1),3);
        Y(1,2) = coord(conne(j,2),3);
        Y(2,1) = coord(conne(j,1),3);
        Y(2,2) = coord(conne(j,3),3);

        Z(1,1) = elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0001;
        Z(1,2) = elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0000;
        Z(2,1) = elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0001;
        Z(2,2) = elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) - 0.0001;

        contourf(X,Y,Z,'EdgeColor','none');
    end
    for i = 1:length(conne(:,1))
        plot([coord(conne(i,1),2) coord(conne(i,2),2)],[coord(conne(i,1),3) coord(conne(i,2),3)], 'k');
        plot([coord(conne(i,2),2) coord(conne(i,3),2)],[coord(conne(i,2),3) coord(conne(i,3),3)], 'k');
        plot([coord(conne(i,3),2) coord(conne(i,1),2)],[coord(conne(i,3),3) coord(conne(i,1),3)], 'k');
    end
    title({"Element \sigma" + sigmaCoordinateText(sigmaCoordinate)},{"timestep " + num2str(plotTimeVisualisation(t)) + " --> t = " + num2str(time(plotTimeVisualisation(t))) + "s"});
    xlabel('x [dm]');
    ylabel('y [dm]');
    cl = colorbar;
    clim([lowLim highLim]);
    colormap('jet');
    nPlot = nPlot +1;
end

%% stress over time
element = 1;
node = 1;
sigmaTime = zeros(length(time),1);
for t = 1:length(time)
    sigmaTime(t) = elementStress(sigmaCoordinate,element,t);
end

figure(nPlot);
plot(time,sigmaTime,"b");
title("\sigma" + sigmaCoordinateText(sigmaCoordinate) + "(t), element " + num2str(element));
xlabel('time [s]');
ylabel('\sigma' + sigmaCoordinateText(sigmaCoordinate) + ' [N/dm^2]');
nPlot = nPlot +1;

%% Y deformation over time
node = length(deformationTime(:,1))/2;
nodalDeformationTime = deformationTime(node+length(deformationTime(:,1))/2,:);

figure(nPlot);
plot(time,nodalDeformationTime,"b");
title("Y displacement node " + num2str(node));
xlabel('time [s]');
ylabel('Y displacement [dm]');
nPlot = nPlot +1;

%% stress vs Y deformation
figure(nPlot);
set(gca,'FontSize',12);
colororder({'b','r'});
yyaxis left
ax = gca;
ax.YColor = 'b';
plot(time,sigmaTime,"b");
hold on;
ylabel('\sigma' + sigmaCoordinateText(sigmaCoordinate) + ' [N/dm^2]');

yyaxis right
plot(time,nodalDeformationTime,"r");
ylabel('Y displacement [dm]');

xlabel('time [s]');
nPlot = nPlot +1;