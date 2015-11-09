%% Optimize posterior of beta distribution
% Author: Xianyuan Zhan, Purdue University
function [m_betaj, V_betaj]=optBeta(n_obs,opt_x,opt_yj,opt_betaj,opt_epsilonj,betaj_0,invbetaSigma_0)
% Newton Raphson Method
eps=0.001;
maxit=5;
betaj=opt_betaj;
for k=1:maxit
    temp_betaj=betaj;
    G=-invbetaSigma_0*(betaj-betaj_0);
    H=-invbetaSigma_0;
    for i=1:n_obs
        G=G+(-exp(opt_x(i,:)*betaj+opt_epsilonj(i))+opt_yj(i))*opt_x(i,:)';
        H=H-exp(opt_x(i,:)*betaj+opt_epsilonj(i))*opt_x(i,:)'*opt_x(i,:);
    end
    betaj=betaj-H\G;
    delta=sum(abs(betaj-temp_betaj));
    if delta<eps
        break;
    end
end
m_betaj=betaj;
V_betaj=inv(-H);
