%{
EUCLIDpreprocessing.m goal:
prepare input data for EUCLID from the raw DIC data

input:
force, time, position of nodes in the initial configuration and their displacements over time

output:
-conne.txt
    connectivity matrix
    for triangular elements the matrix consists of (number of elements) rows and 4 columns. 
    The 4 columns contain the index of the element and then the index of its component nodes.
    example
    [
    1, 1, 38, 7
    2, 5, 6, 39
    3, 8, 9, 42
    4, 3, 4, 40
    5, 12, 13, 41
    6, 10, 11, 43
    7, 14, 15, 44
    ...
    ]

-coord.csv
    The matrix contains information about the mesh nodes.
    It consists of (number of nodes) rows and 5 columns.
    The 5 columns contain the node index, initial position along the x-axis, 
    initial position along the y-axis, Degrees of Freedom (DoFs) along the 
    vertical direction and DoFs in the horizontal direction.
    The DoFs along the vertical direction are either 1 or 2, 1 if the node is
    constrained in the upper side of the mesh and 2 if it is constrained in the lower side.
    The Dofs along the horizontal direction are either 3 or 4, 3 in the case 
    where the node is constrained in the right-hand side of the mesh and 4 if it is constrained in the left-hand side.
    example
    [
    id,x,y,verticalDofs,horizontalDofs
    1,-10.6284,28.9876,1,4
    2,-10.6211,27.7976,0,4
    3,-10.6120,26.6066,0,4
    4,-10.5983,24.8201,0,4
    5,-10.5847,23.0336,0,4
    6,-10.5801,22.4381,0,4
    7,-10.5756,21.8425,0,4
    ...
    ]

-F.csv
    F contains (number of DoFs) rows and (number of timesteps) columns.
    F is in kN.
    in our case:
    in the first row is the force acting on the upper edge in the vertical direction,
    in the second the force acting on the lower edge in the vertical direction,
    in the third the force acting on the right edge in the horizontal direction 
    and in the fourth the force acting on the left edge in the horizontal direction
    example
    [
    0.00459071,0.00472383,0.004679,0.00465854, ...
    -0.00459071,-0.00472383,-0.004679,-0.00465854, ...
    0,0,0,0, ...
    0,0,0,0, ...
    ]

-time.csv
    time contains the time instant at which the DIC data was acquired. 
    It's (number of timesteps) rows and 1 column
    
    example
    [
    0
    1.06508
    2.07184
    3.0716
    4.06575
    5.06411
    6.079
    ...
    ]

-U.csv
    U is (2 x number of nodes) rows and (number of timesteps) columns.
    the first (number of nodes) rows contain for each node the displacements in the horizontal direction
    while the second (number of nodes) rows contain the displacements along the vertical direction
    [
    0,0.001,0.0012,0.0017, ...
    0,0.0009,0.0012,0.0017, ...
    0,0.0008,0.0011,0.0019, ...
    0,0.0005,0.0012,0.0016, ...
    ...
    ]

----  PROCEDURE  ----

the matrices F.csv and time.csv can be created at any time because we use the input data.
in my matlab code they are prepared at the beginning in %% import data and then the force is arranged in the matrix in %% force

1- Identification of the nodes that are never lost during the entire test. The mesh will be created with respect to them.
2- identification of the order of the nodes in the input file from the DIC.
	These two steps are done before the %% reduced mesh section

the variable smallSpecimenRegion was previously used to have meshes of reduced size, 
now it is always negative so you can jump directly to the else %% total mesh

3- finding the nodes on the edges of the mesh (you could collapse point 2 and 3 into the 
   same function if you create a code that can handle everything on its own)
4- finding the nodes that will compose the mesh (you could start by meshing all the nodes in the DIC, 
   the code is in "if fullNodes"). The next "else" allows you to find a mesh with fewer nodes.
5- composition of the mesh (I did the procedure with Rhino, etc.). You should find some library that does this in Phyton). 
   When creating the mesh, you must generate the matrices conne.txt and coord.csv

6- the U matrix must be created. I go back into the code once I have finished the Rhino, etc procedure to create 
   the mesh in "if coordPart". At this point you just need to take the displacements obtained from the DIC 
   for each node and each timestep and arrange them correctly in the matrix
%}

clc;
clear all;
close all;

%% configuration
firstRUN = true;
meshTop = 29;  % mm
height = 59.2;   % mm
% boundary nodes
dx = 0.2;
dy = 0.2;
left = -10.2+dx;
right = 12.3-dx;
top = meshTop-dy;
bottom = meshTop-height+dy;

smallSpecimenRegion = false;
equilateralElements = true;
plots = true;
export = false;

fullNodes = true;
coordPart = false;
exportCoord = false;

% experiment number
experimentNumber = 6;
% 1 first test
% 2 laser cut test
% 3 ellipsoid hole
% 30 ellipsoid hole complete mesh
% 4 3 holes
% 40 3 holes complete mesh
% 5 bells
% 50 bells complete mesh
% 6 ellipsoid hole with reference


%% import data
% FN3000_50kN [kN] , Wegaufnehmer_WA10 [mm] , Kraft 10kN [kN] , DMS_Axial [µm/m] , Kraft 50kN [kN] , Maschinenweg [mm] , MX840B_CH 7 [V] , Trigger [Volt] , Zeit [s] , D1 min. [mm] , D2 min. [mm] , DeltaL 1 [mm] , DeltaL 2 [mm] , Radius 1 [mm] , Radius 2 [mm]
% experiment 1 has 2 separated csv files for X and Y displacement
if experimentNumber == 1
    force = importdata("experiments\1\EUCLID_lve_pp_rec_1.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folderX = "experiments\1\EUCLID_lve_pp_rec_1_dX\";
    dummy = importdata(folderX + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    % displacement y
    folderY = "experiments\1\EUCLID_lve_pp_rec_1_dY\";

    for t = 2:length(time)
        dummy = importdata(folderX + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4) = dummy.data(:,5);
        dummy = importdata(folderY + fileName1 + num2str(t-1) + fileName2);
        nodes{t}(:,5) = dummy.data(:,5);
    end

elseif experimentNumber == 2
    force = importdata("experiments\2\20241213_EUCLID_lve_pa6_rec_LaserCut_Sand_1Hz.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\2\APA6LCRectSand\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 3
    force = importdata("experiments\3\PA6_WJ_ellipsoid.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\3\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 4
    force = importdata("experiments\4\PA6_WJ_three_holes.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\4\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 5
    force = importdata("experiments\5\PA6_WJ_two_bells.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\5\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 30
    force = importdata("experiments\30\PA6_WJ_ellipsoid.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\30\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 40
    force = importdata("experiments\40\PA6_WJ_three_holes.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\40\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 50
    force = importdata("experiments\50\PA6_WJ_two_bells.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "Flächenkomponente 1_";
    fileName2 = ".000 s.csv";
    % displacement x
    folder = "experiments\50\";
    dummy = importdata(folder + fileName1 + num2str(0) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

    elseif experimentNumber == 6
    force = importdata("experiments\3\PA6_WJ_ellipsoid.csv");
    time = force.data(:,9);  % Zeit [s]
    force = force.data(:,3); % Kraft 10kN [kN]

    nodes = cell(length(time),1);

    fileName1 = "EUCLID_WJ_1_hole_000";
    fileName2 = ".csv";
    % displacement x
    folder = "experiments\PA6_WJ_1_hole_with_reference\EUCLID_WJ_1_hole_DIC\";
    dummy = importdata(folder + fileName1 + num2str(1) + fileName2);
    nodes{1} = dummy.data(:,1:3);  % [ index , X, Y ]
    nodes{1}(:,4:5) = 0;           % [ index , X, Y , deltaX , deltaY ]
    reference = nodes{1};

    for t = 2:length(time)
        if t < 11
            fileName1 = "EUCLID_WJ_1_hole_000";
        elseif t < 101
            fileName1 = "EUCLID_WJ_1_hole_00";
        elseif t < 1001
            fileName1 = "EUCLID_WJ_1_hole_0";
        else
            fileName1 = "EUCLID_WJ_1_hole_";
        end
        dummy = importdata(folder + fileName1 + num2str(t-1) + fileName2);
        nodes{t} = dummy.data(:,1:3);
        nodes{t}(:,4:5) = dummy.data(:,5:6);
    end

    dt = zeros(length(time)-1,1);
    for t = 2:length(time)
        dt(t-1) = time(t)-time(t-1);
    end

end

%% permanent nodes
missedNodes = [];
lastNode = zeros(length(time)-1,1);
for t = 2:length(time)
    lastNode(t-1) = nodes{t}(end,1);
    
    missedNodes = clearNodes(nodes{t}(:,1),missedNodes);
end
missedNodes = unique(missedNodes);

if all(lastNode == lastNode(1))
    disp("last node ok");
else
    disp("last node problem");
end

finalNodes = deleteNodes(reference,missedNodes);
finalNodes = finalNodes(:,1:3);

if firstRUN
    figure(1)
    hold on;
    plot(reference(:,2),reference(:,3),'.b');
    plot(reference(1,2),reference(1,3),'.r');
    plot(reference(2,2),reference(2,3),'xk');
    plot(reference(end,2),reference(end,3),'ob');
    legend("dic","first","second","last");
    axis equal;

    figure(2)
    hold on;
    plot(reference(:,2),reference(:,3),'.b');
    plot(finalNodes(:,2),finalNodes(:,3),'or');
    for i = 1:length(finalNodes(:,1))-1
        plot([finalNodes(i,2) finalNodes(i+1,2)],[finalNodes(i,3) finalNodes(i+1,3)],'k');
    end
    legend("dic","permanent nodes");
    axis equal;

    return
end




%% reduced mesh
if smallSpecimenRegion
% meshTop = min(finalNodes(:,3)) + (max(finalNodes(:,3)) - min(finalNodes(:,3)))/2;
while true
    for i = 1:length(missedNodes)
        if nodes{1}(missedNodes(i),3) > meshTop || nodes{1}(missedNodes(i),3) < meshTop - height

        else
            meshTop = meshTop - height;
            continue
        end
    end
    
    if meshTop < min(finalNodes(:,3))
        disp("impossible to find a region");
        break
    end
    break
end

reducedMesh = [];
for i = 1:length(finalNodes(:,1))
    if finalNodes(i,3) < meshTop && finalNodes(i,3) > meshTop-height
        reducedMesh(end+1,:) = finalNodes(i,:); 
    end
end

reducedTest = cell(t,1);
for i = 1:t
    reducedTest{i} = zeros(length(reducedMesh(:,1)),5);
    for j = 1:length(reducedMesh(:,1))
        reducedTest{i}(j,:) = nodes{i}(nodes{i}(:,1) == reducedMesh(j,1),:);
    end
    reducedTest{i}(:,1) = (1:1:j)';
end

% bottom left
for i = 2:length(reducedMesh(:,1))
    if reducedMesh(i,1) - reducedMesh(i-1,1) > 1
        artificialFix = i;
        break
    end
end

artificialFix = 1;

U = zeros(2*length(reducedMesh(:,1)),t);
for i = 1:t
    reducedTest{i}(:,4) = reducedTest{i}(:,4) - reducedTest{i}(artificialFix,4);
    reducedTest{i}(:,5) = reducedTest{i}(:,5) - reducedTest{i}(artificialFix,5);
    U(1:length(reducedMesh(:,1)),i) = reducedTest{i}(:,4);
    U(length(reducedMesh(:,1))+1:2*length(reducedMesh(:,1)),i) = reducedTest{i}(:,5);
end

% 
rhinoInput = reducedMesh(:,2:3);
% 1 up , 2 down , 3 right , 4 left
% 1 - 3 positive force , 2 - 4 negative force
coord = zeros(length(reducedMesh(:,1)),5);
for i = 1:length(reducedMesh(:,1))
    coord(i,1:3) = [i reducedMesh(i,2) reducedMesh(i,3)];
    if reducedMesh(i,3) > top
        coord(i,4) = 1;
    elseif reducedMesh(i,3) < bottom
        coord(i,4) = 2;
    end
    if reducedMesh(i,2) > right
        coord(i,5) = 3;
    elseif reducedMesh(i,2) < left
        coord(i,5) = 4;
    end
end

%% force
forces = zeros(4,t);
for i = 1:t
    forces(:,i) = [force(i) -force(i) 0 0]';
end

%% plots
if plots
nfigure = 1;

figure(nfigure)
hold on;
axis equal;
% ylim([-35 20]);
plot(nodes{1}(:,2),nodes{1}(:,3),"k.");
for i = 1:length(missedNodes)
    plot(reference(reference(:,1) == missedNodes(i),2),reference(reference(:,1) == missedNodes(i),3),"rx");
end
plot([min(nodes{1}(:,2)) max(nodes{1}(:,2))],[meshTop meshTop],"g");
plot([min(nodes{1}(:,2)) max(nodes{1}(:,2))],[meshTop-height meshTop-height],"g");
plot([min(nodes{1}(:,2)) min(nodes{1}(:,2))],[meshTop meshTop-height],"g");
plot([max(nodes{1}(:,2)) max(nodes{1}(:,2))],[meshTop meshTop-height],"g");
plot(reducedMesh(:,2),reducedMesh(:,3),"og");
nfigure = nfigure +1;

figure(nfigure)
hold on;
axis equal;
plot(nodes{1}(:,2),nodes{1}(:,3),".r");
plot(nodes{100}(:,2),nodes{100}(:,3),".b");
% plot(nodes{500}(:,2),nodes{500}(:,3),".k");
for i = 1:length(missedNodes)
    plot(nodes{1}(nodes{1}(:,1) == missedNodes(i),2),nodes{1}(nodes{1}(:,1) == missedNodes(i),3),"xk");
end
nfigure = nfigure +1;

figure(nfigure)
hold on;
axis equal;
plot(nodes{1}(:,2),nodes{1}(:,3),".r");
plot(nodes{end}(:,2),nodes{end}(:,3),".b");
for i = 1:length(missedNodes)
    plot(nodes{1}(nodes{1}(:,1) == missedNodes(i),2),nodes{1}(nodes{1}(:,1) == missedNodes(i),3),"xr");
end
legend({"reference configuration","last configuration","missed nodes"},'Location','northwest','FontSize',10);

%%
nfigure = nfigure +1;
figure(nfigure)
hold on;
axis equal;
plot(reducedTest{1}(:,2),reducedTest{1}(:,3),".r");
plot(reducedTest{1}(:,2) + reducedTest{end}(:,4),reducedTest{1}(:,3) + reducedTest{end}(:,5),".b");
plot(reducedTest{1}(artificialFix,2),reducedTest{1}(artificialFix,3),"ok");
legend({"reference configuration","last configuration"},'Location','northwest','FontSize',10);
nfigure = nfigure +1;

nfigure = nfigure +1;
figure(nfigure)
hold on;
axis equal;
for i = 1:length(coord(:,1))
    plot(coord(i,2),coord(i,3),".k");
    if coord(i,4) == 1 && coord(i,5) == 3
        plot(coord(i,2),coord(i,3),"sr");
    elseif coord(i,4) == 1 && coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"sm");
    elseif coord(i,4) == 2 && coord(i,5) == 3
        plot(coord(i,2),coord(i,3),"sb");
    elseif coord(i,4) == 2 && coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"sk");
    elseif coord(i,4) == 1 
        plot(coord(i,2),coord(i,3),"^r");
    elseif coord(i,4) == 2
        plot(coord(i,2),coord(i,3),"vr");
    elseif coord(i,5) == 3
        plot(coord(i,2),coord(i,3),">r");
    elseif coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"<r");
    end
end
nfigure = nfigure +1;

end

%% export
if export
    nameFile = "experiments/" + num2str(experimentNumber) + "/mesh.txt";
    fID = fopen(nameFile,'w');
    fprintf(fID,'%.4f,%.4f\n',reducedMesh(:,2:3)');  % matrice
    fclose(fID);

    nameDir = "experiments/" + experimentNumber + "/EUCLIDdata/";
    nomi = ["id" "x" "y" "verticalDoFs" "horizontalDoFs"];
    nameFile = strcat(nameDir,'coord.csv');
    fID = fopen(nameFile,'w');
    fprintf(fID,'%s,%s,%s,%s,%s\n',nomi');  % nomi
    fprintf(fID,'%i,%.4f,%.4f,%i,%i\n',coord');  % matrice
    fclose(fID);

    nameFile = strcat(nameDir + 'U.csv');
    writematrix(U,nameFile,"Delimiter",",");

    % nomi = ["id" "Ux" "Uy"];
    % for i = 1:t
    %     nameFile = strcat(nameDir + 'U_' + num2str(i) + '.csv');
    %     fID = fopen(nameFile,'w');
    %     fprintf(fID,'%s,%s,%s\n',nomi');  % nomi
    %     fprintf(fID,'%i,%.4f,%.4f\n',[(1:1:j)' reducedTest{i}(:,4:5)]');  % matrice
    %     fclose(fID);
    % end

    nameFile = strcat(nameDir + 'F.csv');
    writematrix(forces,nameFile,"Delimiter",",");

    nameFile = strcat(nameDir + 'time.csv');
    writematrix(time,nameFile,"Delimiter",",");

    nameFile = strcat(nameDir + 'dt.csv');
    writematrix(dt,nameFile,"Delimiter",",");
end
else
    %% total mesh
    % experiment number 2
    %           x          y        -     x          y
    % right     8.8489    32.5881   -     8.7984   -31.902
    % bottom    8.5058   -32.3032   -   -14.7488   -32.0349
    % left    -14.874    -31.471    -   -15.1047    32.4226
    % top     -14.7246    32.7489   -     7.9659    32.9518   -   last (8.563     32.9572)

    % experiment number 3
    %           x          y        -     x          y
    % right    10.8687    28.6384   -    11.135    -28.5143
    % bottom   11.1169   -29.1073   -   -11.1475   -29.3294
    % left    -11.1559   -28.7341   -   -11.8295    28.4152
    % top     -11.8365    29.0105   -    10.7085    29.1831

    % experiment number 4
    %           x          y        -     x          y
    % right    11.8633    28.755    -    12.2061   -29.0096
    % bottom   12.181    -29.5989   -   -10.1578   -29.9654
    % left    -10.1842   -29.3708   -   -10.6256    28.3931
    % top     -10.6284    28.9876   -    11.8598    29.3505   -  last (10.7402  29.4906) 
    
    % experiment number 5
    %           x          y        -     x          y
    % right    10.2343    28.0422   -    10.5439   -29.0959
    % bottom   10.3994   -29.5822   -   -11.9067   -29.571
    % left    -11.9102   -28.9758   -   -12.2475    28.1633
    % top     -12.1878    28.7026   -    10.231     28.6385   -  last (9.2358    28.9518)

    % experiment number 5
    %           x          y        -     x          y
    % right    11.4954    29.9098   -    11.5372   -30.0948    
    % bottom   11.8766   -30.6151   -   -11.98     -30.7468
    % left    -11.5491   -30.3475   -   -12.3488    29.923
    % top     -12.0842    30.4717   -    11.1763    30.4478   

    % put x positions
    if experimentNumber == 2
        rightEdge = finalNodes(2:find(finalNodes(:,2) == 8.7984),:);
        bottomEdge = finalNodes(find(finalNodes(:,2) == 8.7984)+1:find(finalNodes(:,2) == -14.7488),:);
        leftEdge = finalNodes(find(finalNodes(:,2) == -14.7488)+1:find(finalNodes(:,2) == -15.1047),:);
        topEdge = finalNodes(find(finalNodes(:,2) == -15.1047)+1:find(finalNodes(:,2) == 7.9659),:);
        topEdge = [topEdge; finalNodes(1,:)];
        topBottomNodes = 17;   %  5 ,  9 , 17
        leftRightNodes = 37;  % 10 , 19 , 37
    elseif experimentNumber == 3 || experimentNumber == 30
        topEdge = finalNodes(find(finalNodes(:,2) == -11.8365):find(finalNodes(:,2) == 10.7085),:);
        rightEdge = finalNodes(find(finalNodes(:,2) == 10.8687):find(finalNodes(:,2) == 11.135),:);
        bottomEdge = finalNodes(find(finalNodes(:,2) == 11.1169):find(finalNodes(:,2) == -11.1475),:);
        leftEdge = finalNodes(find(finalNodes(:,2) == -11.1559):find(finalNodes(:,2) == -11.8295),:);
        
        topBottomNodes = 13;   %  5 ,  9 , 17
        leftRightNodes = 33;  % 10 , 19 , 37
    elseif experimentNumber == 4 || experimentNumber == 40
        topEdge = finalNodes(1:find(finalNodes(:,2) == 11.8598),:);
        rightEdge = finalNodes(find(finalNodes(:,2) == 11.8633):find(finalNodes(:,2) == 12.2061),:);
        bottomEdge = finalNodes(find(finalNodes(:,2) == 12.181):find(finalNodes(:,2) == -10.1578),:);
        leftEdge = finalNodes(find(finalNodes(:,2) == -10.1842):find(finalNodes(:,2) == -10.6256),:);
        topEdge = [finalNodes(find(finalNodes(:,2) == -10.6284):find(finalNodes(:,2) == 10.7402),:); topEdge];
        
        topBottomNodes = 17;   %  5 ,  9 , 17
        leftRightNodes = 37;  % 10 , 19 , 37
    elseif experimentNumber == 5 || experimentNumber == 50
        topEdge = finalNodes(1:find(finalNodes(:,2) == 10.231),:);
        rightEdge = finalNodes(find(finalNodes(:,2) == 10.2343):find(finalNodes(:,2) == 10.5439),:);
        bottomEdge = finalNodes(find(finalNodes(:,2) == 10.3994):find(finalNodes(:,2) == -11.9067),:);
        leftEdge = finalNodes(find(finalNodes(:,2) == -11.9102):find(finalNodes(:,2) == -12.2475),:);
        topEdge = [finalNodes(find(finalNodes(:,3) == 28.7026):find(finalNodes(:,2) == 9.2358),:); topEdge];
        
        topBottomNodes = 17;   %  5 ,  9 , 17
        leftRightNodes = 37;  % 10 , 19 , 37
    elseif experimentNumber == 6
        topEdge = finalNodes(1:find(finalNodes(:,2) == 11.1763),:);
        rightEdge = finalNodes(find(finalNodes(:,2) == 11.4954):find(finalNodes(:,2) == 11.5372),:);
        bottomEdge = finalNodes(find(finalNodes(:,2) == 11.8766):find(finalNodes(:,2) == -11.98),:);
        leftEdge = finalNodes(find(finalNodes(:,2) == -11.5491):find(finalNodes(:,2) == -12.3488),:);
        % topEdge = [finalNodes(find(finalNodes(:,3) == ):find(finalNodes(:,2) == ),:); topEdge];
        
        topBottomNodes = 17;   %  5 ,  9 , 17
        leftRightNodes = 37;  % 10 , 19 , 37
    end

    if fullNodes
        topMesh = topEdge;
        rightMesh = rightEdge; 
        bottomMesh = bottomEdge;
        leftMesh = leftEdge;
        internalMesh = finalNodes;
        internalMesh(ismember(internalMesh(:,1),topEdge(:,1)),:) = [];
        internalMesh(ismember(internalMesh(:,1),rightEdge(:,1)),:) = [];
        internalMesh(ismember(internalMesh(:,1),leftEdge(:,1)),:) = [];
        internalMesh(ismember(internalMesh(:,1),bottomEdge(:,1)),:) = [];

        idealFineMesh = [topMesh; rightMesh; bottomMesh; leftMesh; internalMesh];

    else
        topMesh = zeros(topBottomNodes,length(topEdge(1,:)));
        for i = 1:topBottomNodes
            if i == 1
                topMesh(i,:) = topEdge(1,:);
            elseif i == topBottomNodes
                topMesh(i,:) = topEdge(end,:);
            else
                topMesh(i,:) = topEdge(floor((length(topEdge(:,1))/(topBottomNodes-1))*(i-1)),:);
            end
        end

        bottomMesh = zeros(topBottomNodes,length(bottomEdge(1,:)));
        for i = 1:topBottomNodes
            if i == 1
                bottomMesh(i,:) = bottomEdge(1,:);
            elseif i == topBottomNodes
                bottomMesh(i,:) = bottomEdge(end,:);
            else
                bottomMesh(i,:) = bottomEdge(floor((length(bottomEdge(:,1))/(topBottomNodes-1))*(i-1)),:);
            end
        end

        leftMesh = zeros(leftRightNodes-2,length(leftEdge(1,:)));
        for i = 1:leftRightNodes
            if i == 1 || i == leftRightNodes
            else
                leftMesh(i-1,:) = leftEdge(floor((length(leftEdge(:,1))/(leftRightNodes-1))*(i-1)),:);
            end
        end

        rightMesh = zeros(leftRightNodes-2,length(rightEdge(1,:)));
        for i = 1:leftRightNodes
            if i == 1 || i == leftRightNodes
            else
                rightMesh(i-1,:) = rightEdge(floor((length(rightEdge(:,1))/(leftRightNodes-1))*(i-1)),:);
            end
        end

        internalMesh = [];

        for i = 1:(leftRightNodes-2)
            for j = 1:(topBottomNodes-2)
                [x,y] = linesIntersection(bottomMesh(end-(j),2),bottomMesh(end-(j),3),...
                    bottomMesh(end-(j),2),bottomMesh(end-(j),3),...
                    topMesh(j+1,2),topMesh(j+1,3),...
                    topMesh(j+1,2),topMesh(j+1,3),...
                    leftMesh(i,2),leftMesh(i,3),...
                    leftMesh(i,2),leftMesh(i,3),...
                    rightMesh(end-(i-1),2),rightMesh(end-(i-1),3),...
                    rightMesh(end-(i-1),2),rightMesh(end-(i-1),3));
                internalMesh(end+1,:) = [(i-1)*(leftRightNodes-2)+j x y];
                % internalMesh(end+1,:) = [(i-1)*(leftRightNodes-2)+j (topMesh(1+j,2)+bottomMesh(end-j,2))/2 (rightMesh(i,3)+leftMesh(end-(i-1)))/2];
            end
        end

        if equilateralElements
            levels = leftRightNodes-1;
            levelNodes = topBottomNodes-1;
            for i = 1:levels
                for j = 1:levelNodes
                    if i == 1
                        % linesIntersection(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8)
                        [x,y] = linesIntersection(bottomMesh(end-(j-1),2),bottomMesh(end-(j-1),3),...
                            bottomMesh(end-(j)  ,2),bottomMesh(end-(j)  ,3),...
                            topMesh(j  ,2),topMesh(j  ,3),...
                            topMesh(j+1,2),topMesh(j+1,3),...
                            bottomMesh(end,2),bottomMesh(end,3),...
                            leftMesh(1,2),leftMesh(1,3),...
                            bottomMesh(1,2),bottomMesh(1,3),...
                            rightMesh(end,2),rightMesh(end,3));
                        % internalMesh(end+1,:) = [1 ((bottomMesh(j,2)+bottomMesh(j+1,2))/2) bottomMesh(end,3)+abs((leftMesh(i,3)-bottomMesh(end,3))/2)];
                    elseif i == levels
                        [x,y] = linesIntersection(bottomMesh(end-(j-1),2),bottomMesh(end-(j-1),3),...
                            bottomMesh(end-(j)  ,2),bottomMesh(end-(j)  ,3),...
                            topMesh(j  ,2),topMesh(j  ,3),...
                            topMesh(j+1,2),topMesh(j+1,3),...
                            leftMesh(end,2),leftMesh(end,3),...
                            topMesh(1,2),topMesh(1,3),...
                            rightMesh(1,2),rightMesh(1,3),...
                            topMesh(end,2),topMesh(end,3));
                        % internalMesh(end+1,:) = [1 ((bottomMesh(j,2)+bottomMesh(j+1,2))/2) leftMesh(i-1,3)+abs((leftMesh(i-1,3)-topMesh(end,3))/2)];
                    else
                        [x,y] = linesIntersection(bottomMesh(end-(j-1),2),bottomMesh(end-(j-1),3),...
                            bottomMesh(end-(j)  ,2),bottomMesh(end-(j)  ,3),...
                            topMesh(j  ,2),topMesh(j  ,3),...
                            topMesh(j+1,2),topMesh(j+1,3),...
                            leftMesh(i-1,2),leftMesh(i-1,3),...
                            leftMesh(i  ,2),leftMesh(i  ,3),...
                            rightMesh(end-(i-2),2),rightMesh(end-(i-2),3),...
                            rightMesh(end-(i-1),2),rightMesh(end-(i-1),3));
                        % internalMesh(end+1,:) = [1 ((bottomMesh(j,2)+bottomMesh(j+1,2))/2) leftMesh(i-1,3)+abs((leftMesh(i-1,3)-leftMesh(i,3))/2)];
                    end
                    internalMesh(end+1,:) = [1 x y];
                end
            end
        end

        idealFineMesh = [topMesh; rightMesh; bottomMesh; leftMesh; internalMesh];

    end

    disp('finding realFineMesh...');
    realFineMesh = zeros(length(idealFineMesh(:,1)),length(idealFineMesh(1,:)));
    for i = 1:length(realFineMesh(:,1))
        bestNode = finalNodes(1,2:3);
        distanceBestNode = ((finalNodes(1,2)-idealFineMesh(i,2))^2+(finalNodes(1,3)-idealFineMesh(i,3))^2)^0.5;
        for j = 1:length(finalNodes(:,1))
            distance = ((finalNodes(j,2)-idealFineMesh(i,2))^2+(finalNodes(j,3)-idealFineMesh(i,3))^2)^0.5;
            if distance < distanceBestNode
                distanceBestNode = distance;
                bestNode = [finalNodes(j,2) finalNodes(j,3)];
            end
        end
        realFineMesh(i,:) = [i bestNode];
    end

    %% check same nodes
    realFineMesh = unique( realFineMesh(:,[2 3]), 'rows');

    %%
    disp('plot');
    figure(1)
    axis equal;
    hold on;
    plot(finalNodes(:,2),finalNodes(:,3),'.k');
    plot(topEdge(:,2),topEdge(:,3),'^r');
    plot(topMesh(:,2),topMesh(:,3),'ob');
    plot(rightEdge(:,2),rightEdge(:,3),'>r');
    plot(rightMesh(:,2),rightMesh(:,3),'ob');
    plot(bottomEdge(:,2),bottomEdge(:,3),'vr');
    plot(bottomMesh(:,2),bottomMesh(:,3),'ob');
    plot(leftEdge(:,2),leftEdge(:,3),'<r');
    plot(leftMesh(:,2),leftMesh(:,3),'ob');   
    plot(internalMesh(:,2),internalMesh(:,3),'om',MarkerSize=10);
    plot(realFineMesh(:,1),realFineMesh(:,2),'xb',MarkerSize=10);
    legend({"nodes","boundary","externalNodes","","","","","","","idealNodes","finalNodes"},'Location','northwest','FontSize',10);
end

if coordPart
nameFile = "experiments/" + num2str(experimentNumber) + "/nodes.txt";
mesh = importdata(nameFile);

index = zeros(length(mesh(:,1)),1);
for i = 1:length(mesh(:,1))
    for j = 1:length(nodes{1}(:,1))
        if mesh(i,2) == nodes{1}(j,2) && mesh(i,3) == nodes{1}(j,3)
            index(i) = nodes{1}(j,1);
        end
    end
end

U = zeros(2*length(mesh(:,1)),t);
for i = 1:t
    for j = 1:length(mesh(:,1))
        U(j,i) = nodes{i}(nodes{1}(:,1) == index(j),4) - nodes{i}(nodes{1}(:,1) == bottomMesh(end,1),4);
        U(length(mesh(:,1))+j,i) = nodes{i}(nodes{1}(:,1) == index(j),5) - nodes{i}(nodes{1}(:,1) == bottomMesh(end,1),5);
    end
end

%
% 1 up , 2 down , 3 right , 4 left
% 1 - 3 positive force , 2 - 4 negative force
leftMesh(end+1,:) = topMesh(1,:);
leftMesh(end+1,:) = bottomMesh(end,:);
rightMesh(end+1,:) = topMesh(end,:);
rightMesh(end+1,:) = bottomMesh(1,:);

leftEdge(end+1,:) = topEdge(1,:);
leftEdge(end+1,:) = bottomEdge(end,:);
rightEdge(end+1,:) = topEdge(end,:);
rightEdge(end+1,:) = bottomEdge(1,:);

% coord = zeros(length(mesh(:,1)),5);
% for i = 1:length(mesh(:,1))
%     coord(i,1:3) = [i mesh(i,2) mesh(i,3)];
%     if ismember(coord(i,2),topMesh(:,2)) && ismember(coord(i,3),topMesh(:,3))
%         coord(i,4) = 1;
%     end
%     if ismember(coord(i,2),bottomMesh(:,2)) && ismember(coord(i,3),bottomMesh(:,3))
%         coord(i,4) = 2;
%     end
%     if ismember(coord(i,2),rightMesh(:,2)) && ismember(coord(i,3),rightMesh(:,3))
%         coord(i,5) = 3;
%     end
%     if ismember(coord(i,2),leftMesh(:,2)) && ismember(coord(i,3),leftMesh(:,3))
%         coord(i,5) = 4;
%     end
% end

coord = zeros(length(mesh(:,1)),5);
for i = 1:length(mesh(:,1))
    coord(i,1:3) = [i mesh(i,2) mesh(i,3)];
    if ismember(coord(i,2),topEdge(:,2)) && ismember(coord(i,3),topEdge(:,3))
        coord(i,4) = 1;
    end
    if ismember(coord(i,2),bottomEdge(:,2)) && ismember(coord(i,3),bottomEdge(:,3))
        coord(i,4) = 2;
    end
    if ismember(coord(i,2),rightEdge(:,2)) && ismember(coord(i,3),rightEdge(:,3))
        coord(i,5) = 3;
    end
    if ismember(coord(i,2),leftEdge(:,2)) && ismember(coord(i,3),leftEdge(:,3))
        coord(i,5) = 4;
    end
end

%% force
forces = zeros(4,t);
for i = 1:t
    forces(:,i) = [force(i) -force(i) 0 0]';
end

%% plots
nfigure = 100;
%
figure(nfigure)
hold on;
axis equal;
plot(mesh(:,2),mesh(:,3),".r");
plot(mesh(:,2) + U(1:length(mesh(:,1)),50),mesh(:,3) + U(length(mesh(:,1))+1:end,50),".m");
plot(mesh(:,2) + U(1:length(mesh(:,1)),t),mesh(:,3) + U(length(mesh(:,1))+1:end,t),".b");
legend({"reference configuration","timestep 50 configuration","last configuration"},'Location','northwest','FontSize',10);
nfigure = nfigure +1;

figure(nfigure)
hold on;
axis equal;
for i = 1:length(coord(:,1))
    plot(coord(i,2),coord(i,3),".k");
    if coord(i,4) == 1 && coord(i,5) == 3
        plot(coord(i,2),coord(i,3),"sr");
    elseif coord(i,4) == 1 && coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"sm");
    elseif coord(i,4) == 2 && coord(i,5) == 3
        plot(coord(i,2),coord(i,3),"sb");
    elseif coord(i,4) == 2 && coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"sk");
    elseif coord(i,4) == 1 
        plot(coord(i,2),coord(i,3),"^r");
    elseif coord(i,4) == 2
        plot(coord(i,2),coord(i,3),"vr");
    elseif coord(i,5) == 3
        plot(coord(i,2),coord(i,3),">r");
    elseif coord(i,5) == 4
        plot(coord(i,2),coord(i,3),"<r");
    end
end
nfigure = nfigure +1;

end
if exportCoord
    % export
    nameFile = "experiments/" + num2str(experimentNumber) + "/mesh.txt";
    fID = fopen(nameFile,'w');
    fprintf(fID,'%.4f,%.4f\n',mesh(:,2:3)');  % matrice
    fclose(fID);

    nameDir = "experiments/" + experimentNumber + "/EUCLIDdata/";
    nomi = ["id" "x" "y" "verticalDoFs" "horizontalDoFs"];
    nameFile = strcat(nameDir,'coord.csv');
    fID = fopen(nameFile,'w');
    fprintf(fID,'%s,%s,%s,%s,%s\n',nomi');  % nomi
    fprintf(fID,'%i,%.4f,%.4f,%i,%i\n',coord');  % matrice
    fclose(fID);

    nameFile = strcat(nameDir + 'U.csv');
    writematrix(U,nameFile,"Delimiter",",");

    % nomi = ["id" "Ux" "Uy"];
    % for i = 1:t
    %     nameFile = strcat(nameDir + 'U_' + num2str(i) + '.csv');
    %     fID = fopen(nameFile,'w');
    %     fprintf(fID,'%s,%s,%s\n',nomi');  % nomi
    %     fprintf(fID,'%i,%.4f,%.4f\n',[(1:1:j)' reducedTest{i}(:,4:5)]');  % matrice
    %     fclose(fID);
    % end

    nameFile = strcat(nameDir + 'F.csv');
    writematrix(forces,nameFile,"Delimiter",",");

    nameFile = strcat(nameDir + 'time.csv');
    writematrix(time,nameFile,"Delimiter",",");

    % nameFile = strcat(nameDir + 'dt.csv');
    % writematrix(dt,nameFile,"Delimiter",",");
end

%%
figure(200)
hold on;
axis equal;
plot(finalNodes(:,2),finalNodes(:,3),'.k');
plot(idealFineMesh(:,2),idealFineMesh(:,3),".r");
xlabel("x [mm]");
ylabel("y [mm]");
legend({"DIC nodes"},'Location','northwest','FontSize',10);

figure(201)
hold on;
axis equal;
plot(realFineMesh(:,1),realFineMesh(:,2),'.k');
plot(topMesh(:,2),topMesh(:,3),".r");
plot(bottomMesh(:,2),bottomMesh(:,3),".r");
plot(leftMesh(:,2),leftMesh(:,3),"xb");
plot(rightMesh(:,2),rightMesh(:,3),"xb");
xlabel("x [mm]");
ylabel("y [mm]");
legend({"DIC nodes"},'Location','northwest','FontSize',10);

if coordPart
    conne = importdata("experiments/" + num2str(experimentNumber) + "/conne.txt");
    figure(202)
    hold on;
    axis equal;
    for i = 1:length(conne(:,1))
        plot([coord(conne(i,1+1),2) coord(conne(i,2+1),2)],[coord(conne(i,1+1),3) coord(conne(i,2+1),3)], 'k');
        plot([coord(conne(i,2+1),2) coord(conne(i,3+1),2)],[coord(conne(i,2+1),3) coord(conne(i,3+1),3)], 'k');
        plot([coord(conne(i,3+1),2) coord(conne(i,1+1),2)],[coord(conne(i,3+1),3) coord(conne(i,1+1),3)], 'k');
    end
    plot(realFineMesh(:,1),realFineMesh(:,2),'.r');
    plot(mesh(:,2),mesh(:,3),'ob');

    figure(203)
    hold on;
    axis equal;
    for i = 1:length(conne(:,1))
        plot([coord(conne(i,1+1),2) coord(conne(i,2+1),2)],[coord(conne(i,1+1),3) coord(conne(i,2+1),3)], 'k');
        plot([coord(conne(i,2+1),2) coord(conne(i,3+1),2)],[coord(conne(i,2+1),3) coord(conne(i,3+1),3)], 'k');
        plot([coord(conne(i,3+1),2) coord(conne(i,1+1),2)],[coord(conne(i,3+1),3) coord(conne(i,1+1),3)], 'k');
    end
    plot(realFineMesh(:,1),realFineMesh(:,2),'.r');
end