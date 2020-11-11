%These scripts are Written by
% Srikanth Ryali, PhD & Vinod Menon, PhD
%Department of Psychiatry & Behavioral Sciences
%Stanford Cognitive and Systems Neuroscience Laboratory
%Stanford School of Medicine
%Stanford, USA


Vm(:,1) = stim_design';
Vm(:,2) = 1-stim_design';
L = round(32/TR);          % Embedded Dimension
method = 'L1_woi';
Nsubjects = 1;
Ntimepoints = size(timeseries,2);
for subj = 1:Nsubjects
    %subj
    Y = timeseries;   
    J = size(Vm,2);        % Experimental conditions
    [M N] = size(Y); % M: region number N: time length
    v = zeros(N,1);  % External stimulus
    %v = 1-stim_design;
    cnt = 1;
    S = 1;
    for s = 1:S
        data(s).Y = Y;
        data(s).Vm = Vm;
    end
    %%%%%%%%%%% Initialization %%%%%%%%%%%%%%
    Vmc = [];
    vc = [];
    Xc = [];
    Yc = [];
%   [Phi, Bv] = get_model_hrf_rev(L,TR,M,1/1000);
 [Phi, Bv] =  get_model_hrf_3basis(L,TR,M,1/1000);
    for s = 1:S
        Vmc = [Vmc;data(s).Vm];
        vc = [vc;v];
        Y = data(s).Y(:,1:N);
        Yc = [Yc Y];
    end
    for m = 1:M
        y = Yc(m,:);
        R(m,m) = var(y);
    end
    R = eye(M);
    % Weiner Deconvolution
    h = Bv(:,1);
    Xest = weiner_deconv(Yc,h,R);
    [Bm,d,Q] = estAR_wo_intrinsic(Xest',vc,Vmc,1);
    A = zeros(M);
    [B,R1] = initialize_B_R_WD(Yc,Xest,Bv',L);
    %%%%%%%%%%%%%%%% Data for KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = zeros(M,N,S);
    Um = zeros(N,J,S);
    Ue = zeros(N,S);
    for s = 1:S
        Y(:,:,s) = data(s).Y(:,1:N);
        Vm = data(s).Vm;
        Um(:,:,s) = Vm;
        Ue(:,s) = v;
    end
    BDS.Phi = Phi;
    BDS.Phit = Bv';
    BDS.A = A;
    BDS.Bm = Bm;   %Weights due to Modulatory Inputs
    BDS.Q = Q;  %State Covariance
    BDS.d = d;
    BDS.B = B;  %Beta
    BDS.R = R;  %Ouput Covariance
    BDS.L = L; % Embedded Dimension
    BDS.M = M; %# of regionssubj_model_parameters
    BDS.mo = zeros(L*M,1); %Initial state vector mean
    BDS.Vo =  10^-3*eye(L*M);    % Initial Covariance
     BDS.ao = 10^-10; BDS.bo = 10^-10;
     BDS.co = 10^-10; BDS.do = 10^-10;
%      BDS.ao = 10^-0; BDS.bo = (10^-4);
%      BDS.co = 10^-0; BDS.do = (10^-4);
   % BDS.ao = 0; BDS.bo = 0;
   % BDS.co = 0; BDS.do = 0;
    BDS.method = method;
    switch BDS.method
        case 'L1'
            if sum(v(:)) ~= 0
                BDS.Alphab = 0*ones(M,(J+1)*M+1); %L1
            else
                BDS.Alphab = 0*ones(M,(J+1)*M); %L1
            end
            BDS.Alphab_op = 10^-0*ones(M,size(Bv,2)); %L1
        case 'L2'
            BDS.Alphab = 0*ones(M,1); %L2
            BDS.Alphab_op = 0*ones(M,1); %L2
        case 'L1_woi'
            BDS.Alphab = 0*ones(M,(J)*M+1); %L1
            BDS.Alphab_op = 0*ones(M,size(Bv,2)); %L1
        case 'L2_woi'
            BDS.Alphab = 0*ones(M,1); %L2
            BDS.Alphab_op = 0*ones(M,1); %L2
    end
    BDS.flag = 0;
    BDS.tol = tol;
    BDS.maxIter = maxIter;  %Max allowable Iterations
    %[BDS,LL] = vb_em_iterations_all_subjs(BDS,Y,Um,Ue);
    [BDS,LL,KS] = vb_em_iterations_combined_par_convergence(BDS,Y,Um,Ue);
    length(LL)
    BDS.Var_A = ones(M);
    [Mapi_normal,Mapm_normal,Mapd_normal] = compute_normalized_stats_Multiple_inputs(BDS.A,BDS.Bm,BDS.Var_A,BDS.Var_Bm,BDS.d,BDS.Var_d);
    
    subj_model_parameters(subj).Theta_normal = Mapm_normal;
    subj_model_parameters(subj).Theta = BDS.Theta;
    subj_model_parameters(subj).Cov_mat = BDS.Cov_mat;
    subj_model_parameters(subj).Q = BDS.Q;
    
end
% subplot(2,1,1)
% plot(Y')
% legend('M1','Thalamus')
% hold on
% plot(stim_design, 'r','Linewidth',2)
%subplot(2,1,2)
%plot(LL,'o-')
% Results of MDS
%save(fname_result,'subj_model_parameters');
%%%
Mapm_normal(:,:,1);
Pvals(:,:,1) = 1-normcdf(abs(Mapm_normal(:,:,1)));
Pvals(:,:,2) = 1-normcdf(abs(Mapm_normal(:,:,2)));
Map = Pvals <= 0.05/(M*(M-1));%



% if M == 3
%     save(fname_result,'subj_model_parameters','Mapm_normal','Pvals','BDS','KS');
% end
