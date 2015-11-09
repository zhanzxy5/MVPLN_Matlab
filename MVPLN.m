%% MVPLN model
% Author: Xianyuan Zhan, Purdue University

%% Prepare data
Nselect=length(selection);
NumVars=Nselect;
n_obs=length(data(:,1));
X=zeros(n_obs,Nselect);
for i=1:Nselect
    if selection(i)==0
        X(:,1)=ones(n_obs,1);
    else
        X(:,i)=data(:,selection(i));
    end
end
% Epsilon are an n x 2 matrix
% x are an n x k matrix
% beta are an k x 2 matrix
epsilon=zeros(n_obs,NumSevrty);
beta=zeros(NumVars,NumSevrty);

%% Poisson Regression Getting Initial betas
beta_00=zeros(NumVars,NumSevrty);
for i=1:NumSevrty
    mdl =GeneralizedLinearModel.fit(X(:,2:NumVars),Y(:,i),'Distribution','poisson');
    beta_00(:,i)=table2array(mdl.Coefficients(:,1));
end


%% Hyper Parameters
beta_0=zeros(NumVars,NumSevrty);
betaSigma_0=100*eye(NumVars);
df_Sigma=NumSevrty;
% V_Sigma=inv([0.1 0.005;0.005 0.1]); % Follow Jonathan
V_Sigma=eye(NumSevrty);
invV_Sigma=inv(V_Sigma);
Sigma0=eye(NumSevrty); % Initial test
SE_df=NumSevrty;
beta_df=NumVars-1;

%% MCMC
result_beta=zeros(maxIter,NumVars,NumSevrty);
result_sigma=zeros(maxIter,NumSevrty,NumSevrty);
result_flag=zeros(maxIter,NumSevrty);
result_epsilon=zeros(maxIter,n_obs,NumSevrty);
result_LL=zeros(maxIter,1);%% LL evaluated at current iteration
fprintf('MCMC Running:\n');
% Initializing epsilon
SE_SIGMA=Sigma0;
S_epsilon=epsilon;
for i=1:n_obs
    S_epsilon(i,:)=SampleEpsiloni(X(i,:),Y(i,:),beta_00,S_epsilon(i,:),SE_SIGMA,SE_df);
end
% Initialize beta
% S_beta=beta_00; %start from Poisson
S_beta=beta_00+0.1*(rand(NumVars,NumSevrty)-0.5); %start with disturbance
% S_beta=beta_00+0.5*(rand(NumVars,NumSevrty)-0.5); %start with greater disturbance
% Main loop
tic;
for iter=1:maxIter
    fprintf('Iteration: %d\t',iter);
    S_SIGMA=SampleSIGMA(n_obs,S_epsilon,df_Sigma,invV_Sigma);
%     S_SIGMA
%     fprintf('Sample Sigma finished!\t');
    parfor i=1:n_obs
    %for i=1:n_obs
        S_epsilon(i,:)=SampleEpsiloni(X(i,:),Y(i,:),beta_00,S_epsilon(i,:),S_SIGMA,SE_df);
    end
%     fprintf('Sample Epsilon finished!\n');
    parfor j=1:NumSevrty
    %for j=1:NumSevrty
        [betaj,flag]=SampleBeta(n_obs,X,Y(:,j),S_beta(:,j),S_epsilon(:,j),betaSigma_0,beta_0(:,j),beta_df);
        S_beta(:,j)=betaj;
        S_flag(j)=flag;
%         fprintf('Beta %d sampling result: %d\n',j,flag);
    end
    result_beta(iter,:,:)=S_beta;
    result_sigma(iter,:,:)=S_SIGMA;
    result_flag(iter,:)=S_flag;
    result_epsilon(iter,:,:)= S_epsilon;
    
    % Evaluate LL at current iteration
    % Compute likelihood
    LL=0;
    for i=1:n_obs
        for j=1:NumSevrty
            Lambda_i=exp(X(i,:)*S_beta(:,j)+S_epsilon(i,j));
            LL=LL+log(Lambda_i^Y(i,j)*exp(-Lambda_i)/factorial(Y(i,j)));
        end 
        LL=LL+log(mvnpdf(S_epsilon(i,:),zeros(1,NumSevrty),S_SIGMA));
    end
    fprintf('LL at current iteration: %d\n',LL);
    result_LL(iter)=LL;
end
toc;
if save_intermediate_data == 1
    save('output_beta.mat','result_beta');
    save('output_sigma.mat','result_sigma');
    save('output_flag.mat','result_flag');
    save('output_epsilon.mat','result_epsilon');
    save('output_LL.mat','result_LL');
end

%% Report results
checkResult;
computeDIC;
computeCorr;