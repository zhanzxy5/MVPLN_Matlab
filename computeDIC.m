%% Compute DIC
% Author: Xianyuan Zhan, Purdue University

%% Likelihood at the expectation of paprameters
NN=N-burnin;
n_obs=length(result_epsilon(1,:,1));
m_epsilon=zeros(n_obs,Nsev);
for i=1:n_obs
    for j=1:Nsev
        m_epsilon(i,j)=mean(result_epsilon(:,i,j));
    end
end

% Compute likelihood
LL=0;
for i=1:n_obs
    for j=1:Nsev
        Lambda_i=exp(X(i,:)*m_beta(:,j)+m_epsilon(i,j));
        LL=LL+log(Lambda_i^Y(i,j)*exp(-Lambda_i)/factorial(Y(i,j)));
    end
    LL=LL+log(mvnpdf(m_epsilon(i,:),zeros(1,Nsev),m_sigma));
end
LLPM=LL;
fprintf('Log-likelihood evaluated at the posterior mean of the parameters:\n%d\n',LLPM);

% %% Posterior mean deviance - no input
% % Compute likelihood
% LL=0;
% for n=burnin:N
%     n
%     for i=1:n_obs
%         for j=1:Nsev
%             Lambda_i=exp(X(i,:)*result_beta(n,:,j)'+result_epsilon(n,i,j));
%             LL=LL+log(Lambda_i^Y(i,j)*exp(-Lambda_i)/factorial(Y(i,j)));
%         end 
%         LL=LL+log(mvnpdf(reshape(result_epsilon(n,i,:),[1,Nsev]),zeros(1,Nsev),reshape(result_sigma(n,:,:),[Nsev,Nsev])));
%     end
% end
% Dhat=-2*LL/(N-burnin+1);
% fprintf('Posterior mean deviance:\n%d\n',Dhat);
% 
% 
%% Posterior mean - with input
Dhat=-2*sum(result_LL(burnin:N,1))/(N-burnin+1);
fprintf('Posterior mean deviance:\n%d\n',Dhat);
PD=Dhat+2*LLPM;
fprintf('P_D=%d\n',PD);
%% Compute DIC
DIC=2*Dhat+2*LLPM;
fprintf('DIC:%d\n',DIC);

