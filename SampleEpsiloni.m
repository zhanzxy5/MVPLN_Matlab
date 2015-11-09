%% Sample Epsilon
% Author: Xianyuan Zhan, Purdue University
function S_epsiloni=SampleEpsiloni(SE_xi,SE_yi,SE_beta,SE_epsiloni,SE_SIGMA,SE_df)
Nsev=length(SE_epsiloni);
SE_invSIGMA=inv(SE_SIGMA);
[m_epsiloni, V_epsiloni]=optEpsilon(SE_xi,SE_yi,SE_beta,SE_epsiloni,SE_invSIGMA);
V_epsiloni=(V_epsiloni+V_epsiloni')/2;
N_eps=length(SE_epsiloni);
S=zeros(N_eps,N_eps);
diag_invV=S;
for i=1:N_eps
    S(i,i)=sqrt(V_epsiloni(i,i));
    diag_invV(i,i)=1/S(i,i);
end
C=diag_invV*V_epsiloni*diag_invV;
epsiloni_star_standard=mvtrnd(C,SE_df)';
epsilon_star=m_epsiloni'+S*epsiloni_star_standard;
% Compute the proposed density
prop_new=log(mvtpdf(epsiloni_star_standard',C,SE_df));
SE_epsilon_standard=diag_invV*(SE_epsiloni-m_epsiloni)';
prop_old=log(mvtpdf(SE_epsilon_standard',C,SE_df));
% Compute the posterior density
sum_LL_old=0;
sum_LL_new=0;
for i=1:N_eps
    sum_LL_old=sum_LL_old-exp(SE_xi*SE_beta(:,i)+SE_epsiloni(i))+SE_yi(i)*(SE_xi*SE_beta(:,i)+SE_epsiloni(i));
    sum_LL_new=sum_LL_new-exp(SE_xi*SE_beta(:,i)+epsilon_star(i))+SE_yi(i)*(SE_xi*SE_beta(:,i)+epsilon_star(i));
end
post_old=log(mvnpdf(SE_epsiloni,zeros(1,Nsev),SE_SIGMA))+sum_LL_old;
post_new=log(mvnpdf(epsilon_star',zeros(1,Nsev),SE_SIGMA))+sum_LL_new;
alpha=min((post_new+prop_old)-(post_old+prop_new),0);
if alpha>log(rand())
    S_epsiloni=epsilon_star';
else
    S_epsiloni=SE_epsiloni;
end