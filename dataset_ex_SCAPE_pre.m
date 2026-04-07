%% SCAPE - around 5000 vertices for a shape
% quantative experiments

clc;clear;close all;

%% load data

addpath(genpath('../'))

Spath = ('../Dataset/data_learning/SCAPE_r/off/');

finf = dir([Spath,'*.off']);

methods = {'MWP','GCPR','LOPR','LORD'};

for methodIndex = 1:3
    
    Time = zeros(1,length(finf)-1);
    Errors = [];
    
    k=500;
    
    for i = 2:length(finf)
        
        disp([methods{methodIndex},num2str(i-1)]);
        
        file1 = [Spath,finf(1).name];
        file2 = [Spath,finf(i).name];
        
        S1 = MESH.preprocess(file1, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfFindEdge',true,'IfComputeGeoDist',true,'IfComputeNormals',true);
        S2 = MESH.preprocess(file2, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfFindEdge',true,'IfComputeGeoDist',true,'IfComputeNormals',true);
        
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
        
        vx = S1.surface.VERT; fx = S1.surface.TRIV;
        vy = S2.surface.VERT; fy = S2.surface.TRIV;
        save(['../Results/SCAPE_r/data',num2str(i-1),'.mat'],'vx','vy','fx','fy','shot1','shot2');
        
         %% Shape matching methods
        tic;
        switch methodIndex
            case 1
                %% LPO - ours
                Nf = 5; max_iters = 8;
                [T, C] = MWP_L_v6(S1, S2, k, T0, Nf, max_iters,0.2);
                
            case 2
                %% MWP
                f = 5;
                num_iters = 8;
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
                %% LPO -deform
                Nf = 5;
                max_iters = 10;
                
                [T, ~] = MWP_LOPR_DeformField(S1, S2, 500, T0, Nf, max_iters,0.4);
        end
        
        Time=toc;
        save(['../Results/SCAPE_r/',methods{methodIndex},num2str(i-1),'.mat'],'T','Time');
    end
end
