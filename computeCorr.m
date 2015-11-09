%% Compute Correlation
% Author: Xianyuan Zhan, Purdue University

%% Extract data
[N,Nvar,Nsev]=size(result_beta);
% r_beta=result_beta(burnin:N,:,:);
r_sigma=result_sigma(burnin:N,:,:);

%% Compute Correlation
fprintf('Correlation:\n')
for cor1 = 1:Nsev-1
    for cor2 = cor1+1:Nsev
        fprintf('Correlation between dependent variable %d and %d\n', cor1, cor2);
        corr_sigma=r_sigma(:,cor1,cor2)./sqrt(r_sigma(:,cor1,cor1).*r_sigma(:,cor2,cor2));
        m_corr=mean(corr_sigma);
        sd_corr=std(corr_sigma);
        % 95% HDR
        [f,xi] = ksdensity(corr_sigma);
        %total Area
        bandwidth=xi(2)-xi(1);
        S=bandwidth*sum(f);
        S_25=S*0.025;
        tempSL=0;
        tempSU=0;
        for s=1:100
            tempSL=tempSL+bandwidth*f(s);
            if tempSL>S_25
                lower=xi(s);
                break;
            end
        end
        for s=100:-1:1
            tempSU=tempSU+bandwidth*f(s);
            if tempSU>S_25
                upper=xi(s);
                break;
            end
        end
        fprintf('Mean\tStd\t95_HDR\n');
        fprintf('%d\t%d\t%d\t%d\n',m_corr,sd_corr,lower,upper);
    end
end