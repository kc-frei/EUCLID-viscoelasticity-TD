function [betGnm1,betKnm1] = inverse_problem_input(filename,F0,Nel,NMeG,NMeK,tauG,tauK)
% addpath(genpath('routines'));

%% VISCOELASTIC MATERIAL PROPERTIES IDENTIFICATION FROM TIME DOMAIN SIMULATIONS
load(filename);

% Simulation set up 

dt        = data.dt;         
Tsim      = data.Tsim;
time      = data.time;
NN=size(coord,1);


m = [1 1 1 0]';
PrD =  eye(4) - 1/3*(m*m');
PrV =1/3*(m*m');
Dmu = [2 0 0 0;
    0 2 0 0;
    0 0 2 0;
    0 0 0 1];

% Initialisation of the Beta parameters, they come out from the time integration of the evolutionary equations using the trapezoidal rule
Stest(1).betGnm1  =   zeros(4,length(tauG),Nel);    % Beta parameters for the Shear moduli 
Stest(1).betKnm1  =   zeros(1,length(tauK),Nel);    % Beta parameters for the Bulk moduli 

%% Computing the viscous strain at the first time istant - INITIAL GUESS
nt=1;

% node displacements !!! not deformations
deformation_x=data.deformation_time(1:NN,nt);
deformation_y=data.deformation_time(NN+1:end,nt);

Z=zeros(2*NN,1);
Z(1:2:2*NN-1)=deformation_x;
Z(2:2:2*NN)=deformation_y;

epsG_test = zeros(4,Nel,length(time));
tht_test = zeros(1,Nel,length(time));


    for ie=1:Nel
        %------------------------ FEM formulation -------------------------
        vrtx = conne(ie,1:3); % vertexes array
        X = [coord(vrtx(1),2),coord(vrtx(2),2),coord(vrtx(3),2)]';
        Y = [coord(vrtx(1),3),coord(vrtx(2),3),coord(vrtx(3),3)]';
        N_xi  = [-1 1 0]';
        N_eta = [-1 0 1]';
        Jmat(1,1) = N_xi'*X; Jmat(1,2) = N_xi'*Y;
        Jmat(2,1) = N_eta'*X; Jmat(2,2) = N_eta'*Y;
        dofglob = reshape([2*vrtx-1; 2*vrtx],6,1);
        
        
        detJ = det(Jmat);
        Jm1 = 1/detJ*[Jmat(2,2), -Jmat(1,2); -Jmat(2,1), Jmat(1,1)];
        N1_x = Jm1(1,:)*[N_xi(1) N_eta(1)]';
        N1_y = Jm1(2,:)*[N_xi(1) N_eta(1)]';
        N2_x = Jm1(1,:)*[N_xi(2) N_eta(2)]';
        N2_y = Jm1(2,:)*[N_xi(2) N_eta(2)]';
        N3_x = Jm1(1,:)*[N_xi(3) N_eta(3)]';
        N3_y = Jm1(2,:)*[N_xi(3) N_eta(3)]';
                
        B = [N1_x 0    N2_x 0    N3_x 0;
            0     N1_y 0    N2_y 0    N3_y;
            0     0    0    0    0    0;
            N1_y  N1_x N2_y N2_x N3_y N3_x];
        BD = PrD*B;
        b = m'*B;
        vloce = Z(dofglob);

        epsG_test(:,ie,1)        =       BD*vloce;
        tht_test(1,ie,1)         =       b*vloce;
        eps_tot_test(:,ie)       =       epsG_test(:,ie) + tht_test(:,ie)*m;

        %------------ Update Beta parameters: zero nt=0 -------------------
        
        betGnm1(1:4,1:NMeG,ie) =  Stest(nt).betGnm1(1:4,1:NMeG,ie);
        betKnm1(1,1:NMeK,ie)   =  Stest(nt).betKnm1(1,1:NMeK,ie);

        %------------ Update viscous deformations -------------------------

        epsG(:,ie)                  =        epsG_test(:,ie,1); % St(nt).epsG(:,ie); % 
        epsvG(1:4,1:NMeG,ie,nt)     =        repmat(dt./(2*tauG + dt),1,4)'.*epsG(:,ie) +       betGnm1(1:4,1:NMeG,ie);
        
        tht(:,ie)                   =        tht_test(:,ie,1); % St(nt).tht(:,ie); % 
        thtv(1,1:NMeK,ie,nt)        =        (dt./(2*tauK + dt))'.*tht(:,ie)            +       betKnm1(1,1:NMeK,ie) ;

        %------------ Store variables -------------------------------------
        
        Stest(nt).epsvG(1:4,1:NMeG,ie)      =       epsvG(1:4,1:NMeG,ie,nt);
        Stest(nt).thtv(1,1:NMeK,ie)         =       thtv(1,1:NMeK,ie,nt);
        
        Stest(nt).tht(:,ie)                 =       tht(:,ie);
        Stest(nt).epsG(:,ie)                =       epsG(:,ie);
        Stest(nt).eps_tot(:,ie)             =       epsG(:,ie) + tht(:,ie)*m;
        %------------ Check of the results --------------------------------
        
        Stest(nt).check_epsG(ie)=norm(epsG(:,ie) - epsG_test(:,ie));
        Stest(nt).check_tht(ie)=norm(tht(:,ie) - tht_test(:,ie));
        Stest(nt).check_eps_tot(ie)=norm(eps_tot_test(:,ie)  - Stest(nt).eps_tot(:,ie));

    end

%% Computing the viscous strain FOR nt=1,2,3.... - INITIAL GUESS
for nt=2:length(time)
    deformation_x=data.deformation_time(1:NN,nt);
    deformation_y=data.deformation_time(NN+1:end,nt);

    Z=zeros(2*NN,1);
    Z(1:2:2*NN-1)=deformation_x;
    Z(2:2:2*NN)=deformation_y;

    for ie=1:Nel
        %------------------------ FEM formulation -------------------------
        vrtx = conne(ie,1:3); % vertexes array
        X = [coord(vrtx(1),2),coord(vrtx(2),2),coord(vrtx(3),2)]';
        Y = [coord(vrtx(1),3),coord(vrtx(2),3),coord(vrtx(3),3)]';
        N_xi  = [-1 1 0]';
        N_eta = [-1 0 1]';
        Jmat(1,1) = N_xi'*X; Jmat(1,2) = N_xi'*Y;
        Jmat(2,1) = N_eta'*X; Jmat(2,2) = N_eta'*Y;
        dofglob = reshape([2*vrtx-1; 2*vrtx],6,1);


        detJ = det(Jmat);
        Jm1 = 1/detJ*[Jmat(2,2), -Jmat(1,2); -Jmat(2,1), Jmat(1,1)];
        N1_x = Jm1(1,:)*[N_xi(1) N_eta(1)]';
        N1_y = Jm1(2,:)*[N_xi(1) N_eta(1)]';
        N2_x = Jm1(1,:)*[N_xi(2) N_eta(2)]';
        N2_y = Jm1(2,:)*[N_xi(2) N_eta(2)]';
        N3_x = Jm1(1,:)*[N_xi(3) N_eta(3)]';
        N3_y = Jm1(2,:)*[N_xi(3) N_eta(3)]';

        B = [N1_x 0    N2_x 0    N3_x 0;
            0     N1_y 0    N2_y 0    N3_y;
            0     0    0    0    0    0;
            N1_y  N1_x N2_y N2_x N3_y N3_x];
        BD = PrD*B;
        b = m'*B;
        vloce = Z(dofglob);

        epsG_test(:,ie,nt)        =       BD*vloce;
        tht_test(1,ie,nt)         =       b*vloce;
        eps_tot_test(:,ie)     =       epsG_test(:,ie) + tht_test(:,ie)*m;
        
        %------------- Update Beta parameters -----------------------------

        epsGnm1                            =   epsG_test(1:4,ie,nt-1); % St(nt-1).epsG(1:4,ie); % 
        epsvGnm1                           =   Stest(nt-1).epsvG(1:4,1:NMeG,ie);
        betGnm1(1:4,1:NMeG,ie)             =   repmat(dt./(2*tauG + dt),1,4)'.*(epsGnm1 - epsvGnm1) + repmat(2*tauG./(2*tauG + dt),1,4)'.*epsvGnm1;

        thtnm1                             =   tht_test(1,ie,nt-1); % St(nt-1).tht(1,ie); % 
        thtvnm1                            =   Stest(nt-1).thtv(1,1:NMeK,ie);
        betKnm1(1,1:NMeK,ie)               =   (dt./(2*tauK + dt))'.*(thtnm1 - thtvnm1)             + (2*tauK./(2*tauK + dt))'.*thtvnm1;

        %------------- Update viscous deformations ------------------------

        epsG(:,ie)                     =       epsG_test(1:4,ie,nt); % St(nt).epsG(:,ie);
        epsvG(1:4,1:NMeG,ie,nt)        =       repmat(dt./(2*tauG + dt),1,4)'.*epsG(:,ie)          +       betGnm1(1:4,1:NMeG,ie);

        tht(:,ie)                      =       tht_test(1,ie,nt); % St(nt).tht(:,ie);
        thtv(1,1:NMeK,ie,nt)           =       (dt./(2*tauK + dt))'.*tht(:,ie)                     +       betKnm1(1,1:NMeK,ie) ;

        %------------- Store variables ------------------------------------

        Stest(nt).epsvG(1:4,1:NMeG,ie)     =   epsvG(1:4,1:NMeG,ie,nt);
        Stest(nt).thtv(1,1:NMeK,ie)        =   thtv(1,1:NMeK,ie,nt);

        Stest(nt).tht(:,ie)                =   tht(:,ie);
        Stest(nt).epsG(:,ie)               =   epsG(:,ie);

        Stest(nt).betGnm1(1:4,1:NMeG,ie)   =   betGnm1(1:4,1:NMeG,ie);
        Stest(nt).betKnm1(1,1:NMeK,ie)     =   betKnm1(1,1:NMeK,ie);
        
        Stest(nt).eps_tot(:,ie)             =       epsG(:,ie) + tht(:,ie)*m;
        
        %------------ Check of the results --------------------------------
        
        Stest(nt).check_epsG(ie)=norm(epsG(:,ie) - epsG_test(:,ie));
        Stest(nt).check_tht(ie)=norm(tht(:,ie) - tht_test(:,ie));
        Stest(nt).check_eps_tot(ie)=norm(eps_tot_test(:,ie)  - Stest(nt).eps_tot(:,ie));

    end
    nt
end
%%
betGnm1 = cell(length(time),1);
betKnm1 = cell(length(time),1);
for nt=1:length(time)
    betGnm1{nt} = Stest(nt).betGnm1;
    betKnm1{nt} = Stest(nt).betKnm1;
end

% filename = sprintf('betGnm1_q%2.0f_EL%2.0f_2ME_v2.mat',F0,Nel);
% save(filename, 'betGnm1', '-v7.3');
% filename = sprintf('betKnm1_q%2.0f_EL%2.0f_2ME_v2.mat',F0,Nel);
% save(filename, 'betKnm1', '-v7.3');

end

