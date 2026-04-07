clc;clear;close all;

addpath(genpath('../'));

%% load data

% need process shape - more runtime

Spath = ('../Dataset/TOPKIDS/low resolution/');

finf = dir([Spath,'*.off']);

num_pairs = 200;

pairs = gen_random_pairs(length(finf),num_pairs);

save('../Results/TOPKIDS/ls/pairs.mat','pairs');

methods = {'LPO','MWP','GCPD','ICP','BCICP','ZoomOut','PMF','DiscreteOp','KernalMatch','MWP_SinkHorn','Rcpd','SmoothShells'};

for methodIndex = [1,2,3,4,5,6,10,12]
    
    Time = zeros(1,num_pairs);
    Errors = [];
    
    for i=1:num_pairs
        
        %% data preparation
        
        file1 = [Spath,finf(pairs(i,1)).name]; fname1 = strsplit(file1,'.');
        file2 = [Spath,finf(pairs(i,2)).name]; fname2 = strsplit(file2,'.');
        
        k = 500;
        S1 = MESH.preprocess(file1, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true);
        S2 = MESH.preprocess(file2, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfComputeGeoDist',true);
        
        [S1, S2] = surfaceNorm(S1, S2);
        
        S1.area = sum(calc_tri_areas(S1.surface));
        S2.area = sum(calc_tri_areas(S2.surface));
        
        if strcmp(fname1{2},'/Dataset/TOPKIDS/low resolution/kid00')
            GT1 = 1:S1.nv;
        else
            GT1 = load(['..',fname1{2},'_ref.txt']);GT1 = GT1(:,2);
        end
        if strcmp(fname2{2},'/Dataset/TOPKIDS/low resolution/kid00')
            GT2 = 1:S2.nv;
        else
            GT2 = load(['..',fname2{2},'_ref.txt']);GT2 = GT2(:,2);
        end
        GT = zeros(length(GT1),1);
        for ij=1:length(GT1)
            [~,GT(ij)] = min(abs(GT2-GT1(ij)));
        end
        
        opts.shot_num_bins = 10;
        opts.shot_radius = 5;
        shot1 = calc_shot(S1.surface.VERT', S1.surface.TRIV', 1:S1.nv, opts.shot_num_bins, opts.shot_radius*sqrt(S1.area)/100, 3)';
        shot2 = calc_shot(S2.surface.VERT', S2.surface.TRIV', 1:S2.nv, opts.shot_num_bins, opts.shot_radius*sqrt(S2.area)/100, 3)';
        T0 = knnsearch(shot2, shot1, 'NSMethod','kdtree');
        T0_21 = knnsearch(shot1, shot2, 'NSMethod','kdtree');
        
        disp([methods{methodIndex},'-',num2str(i)]);
        
        D = S2.Gamma;
        
        %% Shape matching methods
        tic;
        switch methodIndex
            case 1
                %% LPO - ours
                Nf = 5; max_iters = 8;
                [T, C] = MWP_L_v6(S1, S2, k, T0, Nf, max_iters);
                
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
                para.beta = 1e-1; num_iters = 3;
                C_fmap=S1.evecs\S2.evecs(T0,:);
                C = fMAP.icp_refine(B2,B1,C_fmap,num_iters);
                T = fMAP.fMap2pMap(B2,B1,C);
                
            case 5
                %% BCICP
                B1 = S1.evecs(:,1:k); Ev1 = S1.evals(1:k);
                B2 = S2.evecs(:,1:k); Ev2 = S2.evals(1:k);
                para.beta = 1e-1; num_iters = 3; skipSize = 10;
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
                para.k_init = 3;
                para.k_step = 1;
                para.k_final = 50;
                [T, C21, all_T12, all_C21] = zoomOut_refine(S1.evecs, S2.evecs, T12_ini, para);
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
                T21 = T0_21; B1 = S1.evecs(:,1:50);  B2 = S2.evecs(:, 1:50);
                C12 = B2\B1(T21,:);
                T21 = knnsearch(B1*C12', B2);
                for k = 2:50
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
                opts = struct;
                opts.shot_num_bins = 10; opts.shot_radius = 5;
                opts.lambda = 1e7; opts.mu = 1; opts.maxIter = 10; 
                opts.k = 3; opts.problemSize = 1000; opts.problemSizeInit = 3000; 
                opts.use_par = false; opts.tX = logspace(log10(500),log10(10),opts.maxIter); opts.tY = logspace(log10(500),log10(10),opts.maxIter); 
                opts.t_it = @(t, k, i) (t(i) ./ 1); opts.partial = false; opts.vis = false; 
                opts.n_evecs = 500; opts.oneIteration = false; opts.convexify = false; opts.rho = 0.5; opts.filter = 'hk';
                X = S1.surface; Y = S2.surface;
                X.n = size(X.VERT, 1);
                X.m = size(X.TRIV, 1);
                Y.n = size(Y.VERT, 1);
                Y.m = size(Y.TRIV, 1);
                X.desc = shot1; Y.desc = shot2;
                matches = run_PMF(X, Y, opts);
                T = matches(:,2);
                
            case 10
                %% MWP_SinkHorn
                f = 5;
                num_iters = 4;
                [T,C]=MWP_SinkHorn(S1,S2,T0,f,num_iters);
                
        end
        Time(i)=toc;
        
        err =  calc_geo_err_sparse(1:S1.nv,T,GT,D);
        Errors = [Errors,mean(err)];
        printSeparator('=');
    end
    
    % print error curves
    Error = Errors;
    thresholds = 0:.001:.1;
    for j=1:length(thresholds)
        curve(j) = 100*sum(Error <= thresholds(j)) / length(Error);
    end
    % handle(methodIndex)=plot(thresholds,curve);hold on;
    % L{methodIndex} = sprintf([methods{methodIndex},'  %.4f%%'],mean(Error));
    % exhibit average time
    fprintf('%s: ave-time: %.4f,\t ave-errors: %.4f\n',methods{methodIndex}, mean(Time),mean(Error))
    save(['../Results/TOPKIDS_',methods{methodIndex},'.mat'],'Error','Time');
end
% legend(handle,L,'Location','southeast');

