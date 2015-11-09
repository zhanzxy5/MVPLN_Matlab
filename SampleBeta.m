%% Sample Beta
function [S_betaj,flag]=SampleBeta(n_obs,SB_x,SB_yj,SB_betaj,SB_epsilonj,SB_SIGMA,betaj_0,beta_df)
invSB_SIGMA=inv(SB_SIGMA);
[m_betaj, V_betaj]=optBeta(n_obs,SB_x,SB_yj,SB_betaj,SB_epsilonj,betaj_0,invSB_SIGMA);
V_betaj=(V_betaj+V_betaj')/2;
N_b=length(SB_betaj);
S=zeros(N_b,N_b);
diag_invV=S;
for i=1:N_b
    S(i,i)=sqrt(V_betaj(i,i));
    diag_invV(i,i)=1/S(i,i);
end
C=diag_invV*V_betaj*diag_invV;
betaj_star_standard=mvtrnd(C,beta_df)';
betaj_star=m_betaj+S*betaj_star_standard;
% Compute the proposed density
prop_new=log(mvtpdf(betaj_star_standard',C,beta_df));
SB_betaj_standard=diag_invV*(SB_betaj-m_betaj);
prop_old=log(mvtpdf(SB_betaj_standard',C,beta_df));
% Compute the posterior density
sum_LL_old=0;
sum_LL_new=0;
for i=1:n_obs
    sum_LL_old=sum_LL_old-exp(SB_x(i,:)*SB_betaj+SB_epsilonj(i))+SB_yj(i)*(SB_x(i,:)*SB_betaj+SB_epsilonj(i));
    sum_LL_new=sum_LL_new-exp(SB_x(i,:)*betaj_star+SB_epsilonj(i))+SB_yj(i)*(SB_x(i,:)*betaj_star+SB_epsilonj(i));
end
post_old=log(mvnpdf(SB_betaj',betaj_0',SB_SIGMA))+sum_LL_old;
post_new=log(mvnpdf(betaj_star',betaj_0',SB_SIGMA))+sum_LL_new;
alpha=min((post_new+prop_old)-(post_old+prop_new),0);
%betaj_star'
if alpha>log(rand())
    S_betaj=betaj_star;
    flag=1;
else
    S_betaj=SB_betaj;
    flag=0;
end