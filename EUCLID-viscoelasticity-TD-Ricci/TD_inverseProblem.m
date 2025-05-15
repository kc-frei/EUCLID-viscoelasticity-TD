%% clearing memory
clc;
disp('(0/8) - Clearing memory');
close all;
clear all;

% MAT05
% mm N s

%% configuration
F0=50;    % N
Nel=72;  % mesh elements 

thickness = 1;

tauMin = 0.01;                  % 10x dt
tauMax = 600;

% Maxwell elements    -->    1 + nG + 1 + nK
nG = 20;                    % deviatoric
nK = 20;                    % volumetric

EUCLIDdataPercentage = 100;%   percentage compared to the experiment
sampling = 1/1;              % remaining data (1/1) (1/2) ...
lassoProcedure = false;      % Lasso - least Square Non Negative
clusteringRange = 0.3;
exportData = false;
lambda_i = 1;
lambda_r = 1;  %10000000

planeStress = true;

realData = true;
experimentNumber = 50;
nameDir = "experiments/" + experimentNumber + "/EUCLIDdata/";

%% import data 
disp('(1/8) - Import data'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if realData
    nodes = importdata(nameDir + "coord.csv")';
    coord = nodes.data(:,1:3);
    conne = importdata(nameDir + "conne.txt");
    conne = conne(:,2:4);
    U = importdata(nameDir + "U.csv");
    Nel = length(conne(:,1));
    F = importdata(nameDir + "F.csv");
    F = F*1000;
else
    filename=sprintf('2ME/Id_TD_q%1.0f_EL%1.0f_2ME_dt0.500_Tsim300.000.mat',F0,Nel);
    load(filename);
end


%% variables
disp('(2/8) - Variables creation'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if realData
    time = importdata(nameDir + "time.csv")';
    
    dt = ones(length(time),1);
    for t = 2:length(time)
        dt(t) = time(t)-time(t-1);
    end

    boolTime = true(length(time),1);
    timesteps = sum(boolTime);
else
    time = data.time;
    time = time(1:ceil(length(time)*EUCLIDdataPercentage/100));
    boolTime = false(length(time),1);
    for t = 1:length(time)
        if t == 1
            boolTime(t) = 1;
            continue
        end
        if mod(t-1,(sampling)^-1) == 0
            boolTime(t) = true;
        end
    end
    time = time(boolTime);
    dt = ones(length(time),1);
    dt(:) = time(2) - time(1); % data.dt;                                <-- check
    timesteps = sum(boolTime);
end


tau_G = logspace(log10(tauMin),log10(tauMax),nG)'; %tau_G = [5;10];%
tau_K = logspace(log10(tauMin),log10(tauMax),nK)'; %[7;15];%

% tau_G = logspace(floor(log10(tauMin)),ceil(log10(tauMax)),nG)';
% tau_K = logspace(floor(log10(tauMin)),ceil(log10(tauMax)),nK)';

% tau_G = linspace(tauMin,tauMax,nG)';
% tau_K = linspace(tauMin,tauMax,nK)';

% betaDevElement and betaVolElement creaation
disp('(3/8) - Beta creation'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if realData
    if planeStress
        [betGnm1_, betKnm1_] = inverse_problem_input_realData_noEps33(nameDir,time,dt,coord,U,conne,Nel,nG,nK,tau_G,tau_K);
    else
        [betGnm1_, betKnm1_] = inverse_problem_input_realData(nameDir,time,dt,coord,U,conne,Nel,nG,nK,tau_G,tau_K);
    end
else
    [betGnm1_, betKnm1_] = inverse_problem_input(filename,F0,Nel,nG,nK,tau_G,tau_K);
end

disp('(4/8) - Other variables creation'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nElements = length(conne(:,1));
nNodes = length(coord(:,1));
nNodesElement = length(conne(1,:));

if planeStress
    m = [1 1 0]';
    Idev = eye(3)-(1/2)*(m*m');
    Ivol = (1/2)*(m*m');
    Dmu = [2 0 0;
           0 2 0;
           0 0 1];
else
    m = [1 1 1 0]';
    Idev = eye(4)-(1/3)*(m*m');
    Ivol = (1/3)*(m*m');
    Dmu = [2 0 0 0;
           0 2 0 0;
           0 0 2 0;
           0 0 0 1];
end


epsElement = cell(1,nElements);
epsDevElement = cell(1,nElements);                % Dev = deviatoric
epsDevVisElement = cell(nG,nElements);            % G --> Vis = viscous
epsVolElement = cell(1,nElements);                % Vol = volumetric
epsVolVisElement = cell(nK,nElements);            % K --> Vis = viscous
betaDevElement = cell(nG,nElements);
betaVolElement = cell(nK,nElements);
Ue = cell(1,nElements);
Be = cell(1,nElements);
Bd = cell(1,nElements);
detJ = cell(1,nElements);
b = cell(1,nElements);
aeG = cell(1,nElements);
aeK = cell(1,nElements);
ae = cell(1,nElements);
s = cell(1,nElements);
pm = cell(1,nElements);
sigma = cell(1,nElements);

for e = 1:nElements
    if planeStress
        epsElement{e} = zeros(3,timesteps);
        epsDevElement{e} = zeros(3,timesteps);
    else
        epsElement{e} = zeros(4,timesteps);
        epsDevElement{e} = zeros(4,timesteps);
    end

    epsVolElement{e} = zeros(1,timesteps);

    aeG{e} = cell(1,timesteps);
    aeK{e} = cell(1,timesteps);
    ae{e} = cell(1,timesteps);

    for gamma = 1:nG
        if planeStress
            epsDevVisElement{gamma,e} = zeros(3,timesteps);
            betaDevElement{gamma,e} = zeros(3,timesteps);
        else
            epsDevVisElement{gamma,e} = zeros(4,timesteps);
            betaDevElement{gamma,e} = zeros(4,timesteps);
        end
    end
    for gamma = 1:nK
        epsVolVisElement{gamma,e} = zeros(1,timesteps);
        betaVolElement{gamma,e} = zeros(1,timesteps);
    end

    Ue{e} = zeros(2*nNodesElement,timesteps);
    if planeStress
        s{e} = zeros(3,timesteps);
        pm{e} = zeros(3,timesteps);
        sigma{e} = zeros(3,timesteps);
    else
        s{e} = zeros(4,timesteps);
        pm{e} = zeros(4,timesteps);
        sigma{e} = zeros(4,timesteps);
    end
end

% variable filling --------------------------------------------------------
% deformations
% strain --> for each element e11 e22 e33 e12
% dataStrain_time = data.strain_time(:,:,boolTime);
% for e = 1:nElements
%     for t = 1:timesteps
%         epsElement{e}(:,t) = dataStrain_time(:,e,t);
%         epsDevElement{e}(:,t) = Idev*epsElement{e}(:,t);
%         epsVolElement{e}(t) = m'*epsElement{e}(:,t);
%     end
% end

% displacements
if realData
    dataDeformation_time = U;
else
    dataDeformation_time = data.deformation_time(:,boolTime);
end
for t = 1:timesteps
    for e = 1:nElements
        for i = 1:nNodesElement
            Ue{e}(1+2*(i-1),t) = dataDeformation_time(conne(e,i),t);
            Ue{e}(2+2*(i-1),t) = dataDeformation_time(conne(e,i)+nNodes,t);
        end
    end
end

% BGt , BGt_1 , BKt , BKt_1
BGt_1 = zeros(nG+1,timesteps);
for t = 1:timesteps
    for gamma = 2:nG+1
        BGt_1(gamma,t) = 1-(dt(t)/(2*tau_G(gamma-1)+dt(t)));
    end
end
BGt = BGt_1;
BGt(1,:) = 1;

BKt_1 = zeros(nK+1,timesteps);
for t = 1:timesteps
    for gamma = 2:nK+1
        BKt_1(gamma,t) = 1-(dt(t)/(2*tau_K(gamma-1)+dt(t)));
    end
end
BKt = BKt_1;
BKt(1,:) = 1;

% data evaluation in each element
betGnm1 = cell(timesteps,1);
betKnm1 = cell(timesteps,1);
index = 1;
for t = 1:length(boolTime)
    if boolTime(t) == 1
        betGnm1{index} = betGnm1_{t};
        betKnm1{index} = betKnm1_{t};
        index = index +1;
    end
end
for t = 1:timesteps
    for e = 1:nElements
        for gamma = 1:nG
            betaDevElement{gamma,e}(:,t) = betGnm1{t}(:,gamma,e);
        end

        for gamma = 1:nK
            betaVolElement{gamma,e}(t) = betKnm1{t}(1,gamma,e);
        end
    end
end

for e = 1:nElements
    vrtx = conne(e,:); % vertexes array
    X = [coord(vrtx(1),2),coord(vrtx(2),2),coord(vrtx(3),2)]';
    Y = [coord(vrtx(1),3),coord(vrtx(2),3),coord(vrtx(3),3)]';
    N_xi  = [-1 1 0]';
    N_eta = [-1 0 1]';
    Jmat(1,1) = N_xi'*X; Jmat(1,2) = N_xi'*Y;
    Jmat(2,1) = N_eta'*X; Jmat(2,2) = N_eta'*Y;

    w_gp = 0.5;
    detJ_ = det(Jmat);
    Jm1 = 1/detJ_*[Jmat(2,2), -Jmat(1,2); -Jmat(2,1), Jmat(1,1)];
    N1_x = Jm1(1,:)*[N_xi(1) N_eta(1)]';
    N1_y = Jm1(2,:)*[N_xi(1) N_eta(1)]';
    N2_x = Jm1(1,:)*[N_xi(2) N_eta(2)]';
    N2_y = Jm1(2,:)*[N_xi(2) N_eta(2)]';
    N3_x = Jm1(1,:)*[N_xi(3) N_eta(3)]';
    N3_y = Jm1(2,:)*[N_xi(3) N_eta(3)]';

    if planeStress
        Be{e} = [N1_x   0    N2_x   0    N3_x   0;
                  0    N1_y   0    N2_y   0    N3_y;
                 N1_y  N1_x  N2_y  N2_x  N3_y  N3_x];
    else
        Be{e} = [N1_x   0    N2_x   0    N3_x   0;
                  0    N1_y   0    N2_y   0    N3_y;
                  0     0     0     0     0     0;
                N1_y  N1_x  N2_y  N2_x  N3_y  N3_x];
    end

    Bd{e} = Idev*Be{e};
    b{e} = m'*Be{e};
    detJ{e} = detJ_;
end
%--------------------------------------------------------------------------
% ae matrix
disp('(5/8) - ae matrix creation'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for t = 1:timesteps
    for e = 1:nElements
        % evaluation of ae at the beginning
        if t == 1
            aeG{e}{t} = Be{e}'*Dmu*Bd{e}*Ue{e}(:,t)*w_gp*detJ{e}*BGt(:,t)';
            aeK{e}{t} = b{e}'*b{e}*Ue{e}(:,t)*w_gp*detJ{e}*BKt(:,t)';
        else 
            aeG{e}{t} = Be{e}'*Dmu*Bd{e}*Ue{e}(:,t)*w_gp*detJ{e}*BGt(:,t)';
            for gamma = 2:nG+1
                aeG{e}{t}(:,gamma) = aeG{e}{t}(:,gamma) - Be{e}'*Dmu*w_gp*detJ{e}*betaDevElement{gamma-1,e}(:,t-1);
            end

            aeK{e}{t} = b{e}'*b{e}*Ue{e}(:,t)*w_gp*detJ{e}*BKt(:,t)';
            for gamma = 2:nK+1
                aeK{e}{t}(:,gamma) = aeK{e}{t}(:,gamma) - b{e}'*w_gp*detJ{e}*betaVolElement{gamma-1,e}(t-1);
            end

        end

        ae{e}{t} = [aeG{e}{t},aeK{e}{t}];

    end
end

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   boundary   $$$$$
if realData
    % 1 up , 2 down , 3 right , 4 left
    dataTopY_time = F(1,boolTime);
    dataBottomY_time = F(2,boolTime);
    dataRightX_time = F(3,boolTime);
    dataLeftX_time = F(4,boolTime);
    dof_topX = nodes.data(nodes.data(:,4)==1,1);
    dof_topY = dof_topX;
    dof_topX = 2*(dof_topX-1)+1;
    dof_topY = 2*(dof_topY-1)+2;
    dof_bottomX = nodes.data(nodes.data(:,4)==2,1);
    dof_bottomY = dof_bottomX;
    dof_bottomX = 2*(dof_bottomX-1)+1;
    dof_bottomY = 2*(dof_bottomY-1)+2;
    dof_rightX = nodes.data(nodes.data(:,4)==0,:);
    dof_rightX = dof_rightX(dof_rightX(:,5)==3,1);
    dof_rightY = dof_rightX;
    dof_rightX = 2*(dof_rightX-1)+1;
    dof_rightY = 2*(dof_rightY-1)+2;
    dof_leftX = nodes.data(nodes.data(:,4)==0,:);
    dof_leftX = dof_leftX(dof_leftX(:,5)==4,1);
    dof_leftY = dof_leftX;
    dof_leftX = 2*(dof_leftX-1)+1;
    dof_leftY = 2*(dof_leftY-1)+2;
    edgeGDL = {dof_bottomX,dof_bottomY,dof_leftX,dof_leftY,dof_rightX,dof_rightY,dof_topX,dof_topY};
else
    dataRLeftx_time = data.RLeftx_time(:,boolTime);
    dataRBottomx_time = data.RBottomx_time(:,boolTime);
    edgeGDL = {data.dof_bottomX,data.dof_bottomY,data.dof_leftX,data.dof_leftY,data.dof_rightX,data.dof_rightY,data.dof_topX,data.dof_topY};
end

nEdges = length(edgeGDL);                        
nNodesExternalGDL = cellfun('length',edgeGDL);  
externalGDL = [];                               
for i = 1:nEdges                                 
    for j = 1:nNodesExternalGDL(i)   
        if any(externalGDL == edgeGDL{i}(j))

        else
            externalGDL(end+1) = edgeGDL{i}(j);
        end
    end                                          
end                                              
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

A_int = cell(1,timesteps);
R_int = zeros(2*nNodes-length(externalGDL),1);
A_bnd = cell(1,timesteps);
R_bnd = cell(1,timesteps);
index_vector = true(2*nNodes,1);
index_vector(externalGDL) = false;
for t = 1:timesteps
    A_int{t} = zeros(2*nNodes,length(BGt(:,1))+length(BKt(:,1)));
    A_bnd{t} = zeros(nEdges,length(BGt(:,1))+length(BKt(:,1)));
end

% assembly
% A_int , R_int
for t = 1:timesteps
    for e = 1:nElements
        for nN = 1:nNodesElement
            A_int{t}(2*(conne(e,nN)-1)+1:2*(conne(e,nN)-1)+2,:) = A_int{t}(2*(conne(e,nN)-1)+1:2*(conne(e,nN)-1)+2,:) + ae{e}{t}(2*(nN-1)+1,2*(nN-1)+2,:);
        end
    end
    A_int{t} = A_int{t}(index_vector,:);
end
% A_bnd , R_bnd
for t = 1:timesteps
    if realData
        % edgeGDL = {[bottom X],[bottom Y],[left X],[left Y],[right X],[right Y],[top X],[top Y]};
        % force direction: up --> positive, down --> negative
        R_bnd{t}  =  [0;         F(2,t);     0;       0;       0;        0;        0;      F(1,t)]; 
    else
        % edgeGDL = {[bottom X]                ,[bottom Y],[left X]                 ,[left Y],[right X],[right Y],[top X],[top Y]};
        R_bnd{t} = [sum(dataRBottomx_time(:,t));-F0/2;     sum(dataRLeftx_time(:,t));0;       0;        0;        0;      F0/2]; % <---
    end
    for e = 1:nElements
        for nN = 1:nNodesElement
            for r = 1:nEdges
                if(any(fix((1+edgeGDL{r})/2) == conne(e,nN)))
                    A_bnd{t}(r,:) = A_bnd{t}(r,:) + ae{e}{t}(2*(nN-1)+1+mod(1+edgeGDL{r}(1),2),:);
                end
            end
        end
    end
end


% A_exp , R_exp
disp('(6/8) - A_exp and R_exp assembly'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A_exp = [];
R_exp = [];
for t = 1:timesteps
    if  lambda_i == 0
    A_exp = [A_exp ; A_bnd{t} * lambda_r];
    R_exp = [R_exp ; R_bnd{t} * lambda_r];
    else
    A_exp = [A_exp ; A_int{t} * lambda_i];
    R_exp = [R_exp ; R_int * lambda_i];
    A_exp = [A_exp ; A_bnd{t} * lambda_r];
    R_exp = [R_exp ; R_bnd{t} * lambda_r];
    end
end

%% results
disp('(7/8) - Theta evaluation'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if lassoProcedure
    Nlam = 10;
    lammin = -6;
    lammax = -1;    %   -14     2    1000
    Lambda = logspace(lammin,lammax,Nlam);

    [result1,info] = lasso(A_exp,R_exp,'Lambda',Lambda,'MaxIter',10e3,'RelTol',1e-8,'Standardize',true);
else
    [result,resnorm,residual,exitflag,output,lambda] = lsqnonneg(A_exp,R_exp);
end

thetaResult = result;

tau = [0; tau_G; 0; tau_K];
tau(find(result == 0)) = [];
result(find(result == 0)) = [];
firstKindex = find(tau == 0);
firstKindex = firstKindex(end);

betaIndex = find(thetaResult);
betaIndex(betaIndex <= nG+1) = betaIndex(betaIndex <= nG+1) - 1;
betaIndex(betaIndex > nG+1) = betaIndex(betaIndex > nG+1) - nG - 2;

betaElementResult = cell(length(result),nElements);
for e = 1:nElements
    for gamma = 2:firstKindex-1
        betaElementResult{gamma,e}(:,:) = betaDevElement{betaIndex(gamma),e}(:,:);
    end
    for gamma = firstKindex+1:length(result)
        betaElementResult{gamma,e}(:,:) = betaVolElement{betaIndex(gamma),e}(:,:);
    end
end
%% plot preparation
disp('(8/8) - Plots'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% clustering
% Gs
for i = 1:firstKindex-1

    if i == 1 && tau(i) == 0
        if firstKindex == 2
            clusteredTauG = 0;
            clusteredGs = result(1);
            continue
        else
            clusteredTauG = [0;tau(2)];
            clusteredGs = [result(1);result(2)];
            continue
        end
    elseif i == 1
        clusteredTauG = [tau(1)];
        clusteredGs = [result(1)];
        continue
    end

    if abs(clusteredTauG(end) - tau(i)) < clusteringRange*tau(i) && abs(clusteredTauG(end) - tau(i)) > 0
        clusteredTauG(end) = clusteredGs(end)/(result(i)+clusteredGs(end))*clusteredTauG(end) + result(i)/(result(i)+clusteredGs(end))*tau(i);
        clusteredGs(end) = clusteredGs(end) + result(i);
        % disp([tau(i),i,j]);
    elseif abs(clusteredTauG(end) - tau(i)) > 0
        clusteredTauG(end+1) = tau(i);
        clusteredGs(end+1) = result(i);
    end

end

% Ks
for i = firstKindex:length(tau)

    if i == firstKindex && tau(i) == 0
        if firstKindex == length(tau)
            clusteredTauK = 0;
            clusteredKs = result(firstKindex);
            continue
        else
            clusteredTauK = [0;tau(firstKindex+1)];
            clusteredKs = [result(firstKindex);result(firstKindex+1)];
            continue
        end
    elseif i == firstKindex
        clusteredTauK = [tau(firstKindex)];
        clusteredKs = [result(firstKindex)];
        continue
    end

    if abs(clusteredTauK(end) - tau(i)) < clusteringRange*tau(i) && abs(clusteredTauK(end) - tau(i)) > 0
        clusteredTauK(end) = clusteredKs(end)/(result(i)+clusteredKs(end))*clusteredTauK(end) + result(i)/(result(i)+clusteredKs(end))*tau(i);
        clusteredKs(end) = clusteredKs(end) + result(i);
        % disp([tau(i),i,j]);
    elseif abs(clusteredTauK(end) - tau(i)) > 0
        clusteredTauK(end+1) = tau(i);
        clusteredKs(end+1) = result(i);
    end

end

%% plots
if realData
    figure(2)
    subplot(1,2,1);
    title("G");
    hold on;
    % ylim([0 max([target_mat.G;clusteredGs])+1000]);
    % xlim([-1 20]);
    % xlim([-1 max([target_mat.tauG;tau(1:firstKindex-1)])+5]);
    % plot(target_mat.tauG,target_mat.G,"squarer",MarkerSize=10,LineWidth=2);
    plot(tau(1:firstKindex-1),result(1:firstKindex-1),'.b',MarkerSize=5);
    % plot(clusteredTauG,clusteredGs,'.b',MarkerSize=20);
    % plot(0,target_mat.Ginf,'squarer',MarkerSize=10,LineWidth=2);
    legend({"result"},'Location','northwest','FontSize',10);

    subplot(1,2,2);
    title("K");
    hold on;
    % ylim([0 max([target_mat.K;result(firstKindex:end)])+1000]);
    % xlim([-1 20]);
    % xlim([-1 max([target_mat.tauK;tau(firstKindex:end)])+5]);
    % plot(target_mat.tauK,target_mat.K,"squarer",MarkerSize=10,LineWidth=2);
    plot(tau(firstKindex:end),result(firstKindex:end),'.b',MarkerSize=5);
    % plot(clusteredTauK,clusteredKs,'.b',MarkerSize=20);
    % plot(0,target_mat.Kinf,'squarer',MarkerSize=10,LineWidth=2);
    legend({"result"},'Location','northwest','FontSize',10); % 'southeast'
else
    % figure(1)
    % subplot(1,2,1);
    % title("G");
    % hold on;
    % ylim([0 max([target_mat.G;result(1:firstKindex-1)])+500]);
    % xlim([-1 20]);
    % % xlim([-1 max([target_mat.tauG;tau(1:firstKindex-1)])+5]);
    % plot(target_mat.tauG,target_mat.G,"squarer",MarkerSize=10,LineWidth=2);
    % plot(tau(1:firstKindex-1),result(1:firstKindex-1),'.b',MarkerSize=20);
    % plot(0,target_mat.Ginf,'xr',MarkerSize=20,LineWidth=2);
    % legend({"target","result"},'Location','southeast','FontSize',10);
    %
    % subplot(1,2,2);
    % title("K");
    % hold on;
    % ylim([0 max([target_mat.K;result(firstKindex:end)])+500]);
    % xlim([-1 20]);
    % % xlim([-1 max([target_mat.tauK;tau(firstKindex:end)])+5]);
    % plot(target_mat.tauK,target_mat.K,"squarer",MarkerSize=10,LineWidth=2);
    % plot(tau(firstKindex:end),result(firstKindex:end),'.b',MarkerSize=20);
    % plot(0,target_mat.Kinf,'xr',MarkerSize=20,LineWidth=2);
    % legend({"target","result"},'Location','southeast','FontSize',10);

    figure(2)
    subplot(1,2,1);
    title("G");
    hold on;
    ylim([0 max([target_mat.G;clusteredGs])+1000]);
    xlim([-1 20]);
    % xlim([-1 max([target_mat.tauG;tau(1:firstKindex-1)])+5]);
    plot(target_mat.tauG,target_mat.G,"squarer",MarkerSize=10,LineWidth=2);
    plot(tau(1:firstKindex-1),result(1:firstKindex-1),'.b',MarkerSize=5);
    plot(clusteredTauG,clusteredGs,'.b',MarkerSize=20);
    plot(0,target_mat.Ginf,'squarer',MarkerSize=10,LineWidth=2);
    legend({"target","result","clusteredGs"},'Location','northwest','FontSize',10);

    subplot(1,2,2);
    title("K");
    hold on;
    ylim([0 max([target_mat.K;result(firstKindex:end)])+1000]);
    xlim([-1 20]);
    % xlim([-1 max([target_mat.tauK;tau(firstKindex:end)])+5]);
    plot(target_mat.tauK,target_mat.K,"squarer",MarkerSize=10,LineWidth=2);
    plot(tau(firstKindex:end),result(firstKindex:end),'.b',MarkerSize=5);
    plot(clusteredTauK,clusteredKs,'.b',MarkerSize=20);
    plot(0,target_mat.Kinf,'squarer',MarkerSize=10,LineWidth=2);
    legend({"target","result","clusteredKs"},'Location','northwest','FontSize',10); % 'southeast'
end

%% save data
if exportData
    filename=sprintf('Inv_TD_q%1.0f_EL%2.0f_2ME_dt0.500_Tsim300.000.mat',F0,Nel);

    saveData.Ginf = clusteredGs(1);
    saveData.G = clusteredGs(2:end);
    saveData.tauG = clusteredTauG(2:end);
    saveData.Kinf = clusteredKs(1);
    saveData.K = clusteredKs(2:end);
    saveData.tauK = clusteredTauK(2:end);

    saveData.Ginf_noCluster = result(1);
    saveData.G_noCluster = result(2:firstKindex-1);
    saveData.tauG_noCluster = tau(2:firstKindex-1);
    saveData.Kinf_noCluster = result(firstKindex);
    saveData.K_noCluster = result(firstKindex+1:end);
    saveData.tauK_noCluster = tau(firstKindex+1:end);

    save(filename, 'saveData', '-v7.3');
end

%% sigma
disp('(9/8) - sigma'); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
G = zeros(timesteps,1);
K = zeros(timesteps,1);
G(:) = result(1);                 % G_inf
K(:) = result(firstKindex);       % K_inf
for t = 1:timesteps
    for gamma = 2:firstKindex-1
        G(t) = G(t) + result(gamma)*(1-(dt(t)/(2*tau(gamma)+dt(t))));
    end
    for gamma = firstKindex+1:length(result)
        K(t) = K(t) + result(gamma)*(1-(dt(t)/(2*tau(gamma)+dt(t))));
    end
end
for t = 1:timesteps
    for e = 1:nElements
        s{e}(:,t) = G(t)*Dmu*Idev*Be{e}*Ue{e}(:,t);
        pm{e}(:,t) = K(t)*m*b{e}*Ue{e}(:,t);

        for gamma = 2:firstKindex-1
            s{e}(:,t) = s{e}(:,t) - result(gamma)*Dmu*betaElementResult{gamma,e}(:,t);
        end
        for gamma = firstKindex+1:length(result)
            pm{e}(:,t) = pm{e}(:,t) - result(gamma)*m*betaElementResult{gamma,e}(t);
        end

        sigma{e}(:,t) = s{e}(:,t) + pm{e}(:,t);
    end
end

%% plot
plotTimeVisualisation = [3 5 length(time)];
nPlot = 100;
for sigmaCoordinate = 1:3
if planeStress
%     sigmaCoordinate = 2;     % sigma11 = 1 , sigma22 = 2 , sigma12 = 3
    sigmaCoordinateText = ["11" "22" "12"];
else
%     sigmaCoordinate = 2;     % sigma11 = 1 , sigma22 = 2 , sigma33 = 3 , sigma12 = 4
    sigmaCoordinateText = ["11" "22" "33" "12"];
end
elementStress = zeros(3,nElements,timesteps);

for t = 1:timesteps
    for e = 1:nElements
        elementStress(:,e,t) = sigma{e}(:,t);
    end
end

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
end

if ~realData
    lowLim = min(min(data.stress_time(sigmaCoordinate,:,plotTimeVisualisation) - elementStress(sigmaCoordinate,:,plotTimeVisualisation)))-10;
    highLim = max(max(data.stress_time(sigmaCoordinate,:,plotTimeVisualisation) - elementStress(sigmaCoordinate,:,plotTimeVisualisation)))+10;
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

            Z(1,1) = data.stress_time(sigmaCoordinate,j,plotTimeVisualisation(t)) - elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0001;
            Z(1,2) = data.stress_time(sigmaCoordinate,j,plotTimeVisualisation(t)) - elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0000;
            Z(2,1) = data.stress_time(sigmaCoordinate,j,plotTimeVisualisation(t)) - elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) + 0.0001;
            Z(2,2) = data.stress_time(sigmaCoordinate,j,plotTimeVisualisation(t)) - elementStress(sigmaCoordinate,j,plotTimeVisualisation(t)) - 0.0001;

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
end

%% edge forces, directions 11 22 12
% code with no real data is not working due to different variable size
% this is valid also for the next section plot
% if ~realData not updated
% edgeNodes = { {top edge} , {bottom edge} , {left edge} , {right edge} }
if realData
    edgeNodes = {nodes.data(nodes.data(:,4) == 1,1),...
                 nodes.data(nodes.data(:,4) == 2,1),...
                 nodes.data(nodes.data(:,5) == 4,1),...
                 nodes.data(nodes.data(:,5) == 3,1)};
else
    bottomNodes = data.dof_bottomY/2;
end

nodalElements = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
sigma11Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
force11Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
sigma22Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
force22Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
sigma12Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
force12Nodes = {cell(1,length(edgeNodes{1})),cell(1,length(edgeNodes{2})),...
                 cell(1,length(edgeNodes{3})),cell(1,length(edgeNodes{4}))};
area = {zeros(1,length(edgeNodes{1})),zeros(1,length(edgeNodes{2})),zeros(1,length(edgeNodes{3})),zeros(1,length(edgeNodes{4}))};
F_edge_11 = {zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps)};
F_edge_22 = {zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps)};
F_edge_12 = {zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps),zeros(1,timesteps)};

for edge = 1:4
    for i = 1:length(nodalElements{edge})
        nodalElements{edge}{i} = [];
        sigma11Nodes{edge}{i} = [];
        sigma22Nodes{edge}{i} = [];
        sigma12Nodes{edge}{i} = [];
    end

    for e = 1:nElements
        for i = 1:length(edgeNodes{edge})
            if ismember(edgeNodes{edge}(i),conne(e,:))
                nodalElements{edge}{i}(end+1) = e;
            end
        end
    end

    for t = 1:timesteps
        for i = 1:length(nodalElements{edge})
            sigma11Nodes{edge}{i}(end+1) = 0;
            sigma22Nodes{edge}{i}(end+1) = 0;
            sigma12Nodes{edge}{i}(end+1) = 0;
            for j = 1:length(nodalElements{edge}{i})
                sigma11Nodes{edge}{i}(end) = sigma11Nodes{edge}{i}(end) + ...
                                             elementStress(1,nodalElements{edge}{i}(j),t); % 1 --> sigma11
                sigma22Nodes{edge}{i}(end) = sigma22Nodes{edge}{i}(end) + ...
                                             elementStress(2,nodalElements{edge}{i}(j),t); % 2 --> sigma22
                sigma12Nodes{edge}{i}(end) = sigma12Nodes{edge}{i}(end) + ...
                                             elementStress(3,nodalElements{edge}{i}(j),t); % 3 --> sigma12
            end
            sigma11Nodes{edge}{i}(end) = sigma11Nodes{edge}{i}(end) / length(nodalElements{edge}{i});
            sigma22Nodes{edge}{i}(end) = sigma22Nodes{edge}{i}(end) / length(nodalElements{edge}{i});
            sigma12Nodes{edge}{i}(end) = sigma12Nodes{edge}{i}(end) / length(nodalElements{edge}{i});
        end
    end
    
    if edge == 1 || edge == 2
        for i = 1:length(edgeNodes{edge})
            pointsXposition = [];
            for j = 1:length(nodalElements{edge}{i})
                for k = 1:nNodesElement
                    pointsXposition(end+1) = coord(conne(nodalElements{edge}{i}(j),k),2);
                end
            end
            pointsXposition = pointsXposition(ismember(pointsXposition,coord(edgeNodes{edge}(:),2)));
            pointsXposition = sort(unique(pointsXposition));
            area{edge}(i) = thickness*abs(pointsXposition(end) - pointsXposition(1))/2;
        end
    else
        for i = 1:length(edgeNodes{edge})
            pointsYposition = [];
            for j = 1:length(nodalElements{edge}{i})
                for k = 1:nNodesElement
                    pointsYposition(end+1) = coord(conne(nodalElements{edge}{i}(j),k),3);
                end
            end
            pointsYposition = pointsYposition(ismember(pointsYposition,coord(edgeNodes{edge}(:),3)));
            pointsYposition = sort(unique(pointsYposition));
            area{edge}(i) = thickness*abs(pointsYposition(end) - pointsYposition(1))/2;
        end
    end

    for i = 1:length(edgeNodes{edge})
        force11Nodes{edge}{i} = sigma11Nodes{edge}{i}*area{edge}(i);
        force22Nodes{edge}{i} = sigma22Nodes{edge}{i}*area{edge}(i);
        force12Nodes{edge}{i} = sigma12Nodes{edge}{i}*area{edge}(i);
    end

    for t = 1:timesteps
        for i = 1:length(edgeNodes{edge})
            F_edge_11{edge}(t) = F_edge_11{edge}(t) + force11Nodes{edge}{i}(t);
            F_edge_22{edge}(t) = F_edge_22{edge}(t) + force22Nodes{edge}{i}(t);
            F_edge_12{edge}(t) = F_edge_12{edge}(t) + force12Nodes{edge}{i}(t);
        end
        if ~realData
            F_edge_22(t) = F_edge_22(t)*2;
        end
    end
end

%% plot bottom force over time
F0_ = zeros(timesteps);
F0_(:) = F0;
figure(nPlot)
hold on;
plot(time,F_edge_22{2},'r');
xlabel("time [s]");
ylabel("force [N]");
if realData
    ylim([0 10+max(max(-F(2,:),max(F_edge_22{2})))]);
    plot(time,-F(2,:),'b');
    legend({"EUCLID bottom force","experimental force"},'Location','southeast','FontSize',10);
else
    ylim([0 10+max(max(F0,FME))]);
    plot(time,F0_,'b');
    legend({"EUCLID bottom force","F0"},'Location','southeast','FontSize',10);
end
nPlot = nPlot +1;

%% export to txt
% copy variables F_topEdge, F_bottomEdge F_leftEdge, F_rightEdge
F_topEdge    = [F_edge_11{1}; F_edge_22{1}; F_edge_12{1}]';
F_bottomEdge = [F_edge_11{2}; F_edge_22{2}; F_edge_12{2}]';
F_leftEdge   = [F_edge_11{3}; F_edge_22{3}; F_edge_12{3}]';
F_rightEdge  = [F_edge_11{4}; F_edge_22{4}; F_edge_12{4}]';