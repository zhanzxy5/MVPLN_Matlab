%% Optimize posterior of epsilon distribution
% Author: Xianyuan Zhan, Purdue University
function [m_epsiloni, V_epsiloni]=optEpsilon(opt_xi,opt_yi,opt_beta,opt_epsiloni,opt_invSigma)
% Newton Raphson Method
eps=0.001;
maxit=5;
epsiloni=opt_epsiloni;
for k=1:maxit
    temp_epsiloni=epsiloni;
    G=-opt_invSigma*epsiloni'+opt_yi'-exp(opt_xi*opt_beta+epsiloni)';
    H=-opt_invSigma-diag(exp(opt_xi*opt_beta+epsiloni));
    epsiloni=epsiloni-(H\G)';
    delta=sum(abs(epsiloni-temp_epsiloni));
    if delta<eps
        break;
    end
end
m_epsiloni=epsiloni;
V_epsiloni=inv(-H);