%% Sample Sigma
% Author: Xianyuan Zhan, Purdue University
function Sigma=SampleSIGMA(n_obs,S_epsilon,df_Sigma,invV_Sigma)
df=n_obs+df_Sigma;
V=invV_Sigma;
for i=1:n_obs
    V=V+S_epsilon(i,:)'*S_epsilon(i,:);
end
% invSigma=wishrnd(V,df);
% Sigma=inv(invSigma);
% invV=inv(V);
Sigma=iwishrnd(V,df);
