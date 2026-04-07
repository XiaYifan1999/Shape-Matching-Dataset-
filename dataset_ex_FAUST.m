clc;clear;close all;

addpath(genpath('../'));

%% load data
addpath(genpath('../'))

Spath = ('../Dataset/FAUST_PRE/');

num_pairs = 300; 

finf = dir([Spath,'*.mat']);
% pairs = gen_random_pairs(length(finf),num_pairs);
% save('../Results/FAUST_pairs.mat','pairs');
load('../Results/FAUST_pairs.mat');

methods = {'LPO','MWP','GCPD','ICP','BCICP','ZoomOut','PMF','DiscreteOp','KernalMatch','MWP_SinkHorn','Rcpd','SmoothShells','ICP_LPO','LPO-DF','PC-GAU','Shot'};

for methodIndex = 15
    
    Time = zeros(1,num_pairs);
    Errors = [];
    
    for i=1:size(pairs,1)
        
        %% data preparation
        % file1 = finf(pairs(i,1)).name;
        % file2 = finf(pairs(i,2)).name;
        
        disp([methods{methodIndex},'-',num2str(i)]);
        
        k = 500;
        load([Spath,num2str(pairs(i,1)),'.mat']); S1 = S; % MESH.preprocess(file1, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfFindEdge',true,'IfComputeGeoDist',true,'IfComputeNormals',true);
        load([Spath,num2str(pairs(i,2)),'.mat']); S2 = S; % MESH.preprocess(file2, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfFindEdge',true,'IfComputeGeoDist',true,'IfComputeNormals',true);
        
        [S1, S2] = surfaceNorm(S1, S2);
        
        S1.area = sum(calc_tri_areas(S1.surface));
        S2.area = sum(calc_tri_areas(S2.surface));
        
        D = S2.Gamma; % MESH.compute_geodesic_dist_matrix(S2);
        
        %% SHOT matching
        opts.shot_num_bins = 10;
        opts.shot_radius = 5;
        shot1 = calc_shot(S1.surface.VERT', S1.surface.TRIV', 1:S1.nv, opts.shot_num_bins, opts.shot_radius*sqrt(S1.area)/100, 3)';
        shot2 = calc_shot(S2.surface.VERT', S2.surface.TRIV', 1:S2.nv, opts.shot_num_bins, opts.shot_radius*sqrt(S2.area)/100, 3)';
        T0 = knnsearch(shot2, shot1, 'NSMethod','kdtree');
        T0_21 = knnsearch(shot1, shot2, 'NSMethod','kdtree');
        
        %% Shape matching methods
        tic;
        switch methodIndex
            case 1
                %% LPO - ours
                Nf = 5; max_iters = 5;
                [T, C] = MWP_L_v6(S1, S2, k, T0, Nf, max_iters,0.2);
                
            case 2
                %% MWP
                f = 5;
                num_iters = 4;
                [T, C] = MWP2(S1, S2, 500, T0, Nf, num_iters);
                
            case 3
                %% GCPD
                % MWP initialized
                Nf = 5;
                num_iters = 3;
                [~, C] = MWP2(S1, S2, 200, T0, Nf, num_iters);
                k0 = size(C,1);
                T12_0 = knnsearch(S2.evecs(:,1:k0)*C',S1.evecs(:,1:k0));
                
                X0 = S1.surface.VERT;
                Y0 = S2.surface.VERT;
                % GCPD intial
                conf = GCPD_initial_settings([]);
                U = S1.evecs(1:size(S1.surface.VERT,1), 1:k);
                e = S1.evals(1:k);
                [~, C, ~] = GCPD_initial(Y0(T12_0,:), X0, U, e, conf);
                % GCPD refine
                conf = GCPD_settings([]);
                U0 = S1.evecs(:, 1:k);
                [T0, ~, ~] = GCPD(Y0, X0, U0, e, C, conf);
                T = knnsearch(Y0, T0);
                
            case 4
                %% ICP
                B1 = S1.evecs(:,1:k); Ev1 = S1.evals(1:k);
                B2 = S2.evecs(:,1:k); Ev2 = S2.evals(1:k); 
                num_iters = 10;
                C_fmap=S1.evecs\S2.evecs(T0,:);
                C = fMAP.icp_refine(B2,B1,C_fmap,num_iters);
                T = fMAP.fMap2pMap(B2,B1,C);
                
            case 5
                %% BCICP
                B1 = S1.evecs(:,1:k); Ev1 = S1.evals(1:k);
                B2 = S2.evecs(:,1:k); Ev2 = S2.evals(1:k);
                para.beta = 1e-1; num_iters = 10; skipSize = 10;
                C_fmap=S1.evecs\S2.evecs(T0,:);
                T12_direct = fMAP.fMap2pMap(B2,B1,C_fmap);
                T21_direct = fMAP.fMap2pMap(B1,B2,C_fmap');
                [~, T] =  bcicp_refine(S1,S2,B1,B2,T21_direct, T12_direct, num_iters);
                
            case 6
                %% Zoom Out
                B1 = S1.evecs(:,1:4);
                B2 = S2.evecs(:,1:4);
                C21_ini = diag([1,1,1,1]);
                T12_ini = fMAP.fMap2pMap(B2,B1,C21_ini);
                para.k_init = 300;
                para.k_step = 5;
                para.k_final = 500;
                [T, C21, all_T12, all_C21] = zoomOut_refine(S1.evecs, S2.evecs, T0, para);
                % para.num_samples = 200;
                % T12_fast = zoomOut_refine_fast(S1, S2, T12_ini, para,1);
                
            case 7
                %% PMF
                niters = 3; X = S1; Y = S2;
                X.n = length(X.surface.VERT); X.m = length(X.surface.TRIV); [ X.fps , X.sradius ] = metricfps_new( X.n, X.Gamma);
                Y.n = length(Y.surface.VERT); Y.m = length(Y.surface.TRIV); [ Y.fps , Y.sradius ] = metricfps_new( Y.n, Y.Gamma);
                X.tri_areas = facearea(X.surface.VERT,X.surface.TRIV); Y.tri_areas = facearea(Y.surface.VERT,Y.surface.TRIV);
                varx = 2*X.area; vary = 2*Y.area; 
                xin = 1:X.n; yin = T0;
                for i=1:niters
                    F = kernel_density_estimation(X,Y,xin,yin,varx,vary);
                    xin = 1:X.n;
                    yin = sparseAssignmentProblemAuctionAlgorithm(F,[],[],0);
                end
                T = yin;
                
            case 8
                %% Discrete Optimization
                T21 = T0_21; B1 = S1.evecs(:,1:500);  B2 = S2.evecs(:, 1:500);
                C12 = B2\B1(T21,:);
                T21 = knnsearch(B1*C12', B2);
                for k = 400:20:500
                    B11 = S1.evecs(:,1:k); B22 = S2.evecs(:, 1:k);
                    Ev11 = S1.evals(1:k);  Ev22 = S2.evals(1:k);
                    Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); % normalize the delta to enforce the isometry
                    for iter = 1:5
                        C12 = B22\B11(T21,:);
                        T21 = knnsearch(B11*diag(Ev11)*C12', B22*diag(Ev22));
                        
                        C12 = B22\B11(T21,:);
                        T21 = knnsearch(B11*C12', B22);
                    end
                end
                T = knnsearch(B2,B1*C12');
                
            case 9
                %% KernalMatching
                addpath(genpath('../Methods/KernalMatching-master/'));
                opts = struct;
                opts.shot_num_bins = 10; opts.shot_radius = 5;
                opts.lambda = 1e7; opts.mu = 1; opts.maxIter = 10; 
                opts.k = 3; opts.problemSize = 1000; opts.problemSizeInit = 3000; 
                opts.use_par = false; opts.tX = logspace(log10(500),log10(10),opts.maxIter); opts.tY = logspace(log10(500),log10(10),opts.maxIter); 
                opts.t_it = @(t, k, i) (t(i) ./ 1); opts.partial = false; opts.vis = false; 
                opts.n_evecs = 500; opts.oneIteration = false; opts.convexify = false; opts.rho = 0.5; opts.filter = 'hk';
                X = S1; Y = S2; X.TRIV = X.surface.TRIV; Y.TRIV = Y.surface.TRIV; 
                X.VERT = X.surface.VERT; Y.VERT = Y.surface.VERT;
                X.n = size(X.VERT, 1); X.m = size(X.TRIV, 1);
                Y.n = size(Y.VERT, 1); Y.m = size(Y.TRIV, 1);
                X.desc = shot1; Y.desc = shot2;
                X.lores = X; Y.lores = Y;
                [X.evecs,X.evals] = HeatKernels(X,1:X.n,opts.n_evecs);
                [Y.evecs,Y.evals] = HeatKernels(Y,1:Y.n,opts.n_evecs);
                matches = run_PMF(X, Y, opts);
                T = matches(:,2);
                
            case 10
                %% MWP_SinkHorn
                f = 5;
                num_iters = 4;
                [T,C]=MWP_SinkHorn(S1,S2,T0,f,num_iters);
                
            case 11
                %% CPD - regularized recovery
                f = 5;
                num_iters = 4;
                [~, C] = MWP2(S1, S2, 500, T0, Nf, num_iters);
                
                M.VERT = S1.surface.VERT; M.TRIV = S1.surface.TRIV; M.n = S1.nv; M.m = S1.nf; M.S = S1.Gamma;
                N.VERT = S2.surface.VERT; N.TRIV = S2.surface.TRIV; N.n = S2.nv; N.m = S2.nf; N.S = S2.Gamma;
                T = run_cpd(M, N, C, 10);
                
            case 12
                %% Smooth Shells
                param = struct;
                param = standardparams(param);
                
                param.noPlot = true; %turn on/off for intermediate plots
                param.GPUcorrespondences = false; %recommended option = true. Turn off, if no GPU is available
                
                %method parameters with recommended settings
                param.facFeat = 0.25;
                param.kArrayLength = 50;
                param.kMax = 490;
                param.lambdaArapInit = 0.02;
                param.lambdaFeat = 25;
                param.normalDamping = 0.04;
                param.numMCMC = 100;
                
                %execute main script
                [P,tau,C,X,Y] = smoothshells2(S1,S2,param);
                k0 = size(C,1);
                T = knnsearch(S2.evecs(:,1:k0)*C',S1.evecs(:,1:k0));
                
            case 13
                %% ICP - LPO
                num_iters = 10; lambda = 0.02;
                C_fmap=S1.evecs\S2.evecs(T0,:);
                C = ICP_LPO(S2,S1,C_fmap,num_iters,lambda);
                T = fMAP.fMap2pMap(S2,S1,C);
                
            case 14
                %% LPO -deform 
                Nf = 5;
                max_iters = 5;
                
                [T, ~] = MWP_LOPR_DeformField(S1, S2, 500, T0, Nf, max_iters,0.2);
                
            case 15
                %% PC-GAU
                T = PC_Gau_STAG(S1,S2,T0);
                
            otherwise
                T = T0;
                
        end
        Time(i)=toc;
        
        err =  calc_geo_err_sparse(1:S1.nv,T,1:S2.nv,D);
        Errors = [Errors,err];
        printSeparator('=');
    end
    
    % print error curves
    Error = mean(Errors);
    % handle(methodIndex)=plot(thresholds,curve);hold on;
    % L{methodIndex} = sprintf([methods{methodIndex},'  %.4f%%'],mean(Error));
    % exhibit average time
    fprintf('%s: ave-time: %.4f,\t ave-errors: %.4f\n',methods{methodIndex}, mean(Time),mean(Error))
    save(['../Results/FAUST_',methods{methodIndex},'.mat'],'Error','Time');
end
% legend(handle,L,'Location','southeast');

