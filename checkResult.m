%% Check Result
% Author: Xianyuan Zhan, Purdue University

% Extract data
[N,Nvar,Nsev]=size(result_beta);
r_sigma=result_sigma(burnin:N,:,:);

%% Plot
X_iter=[1:N]';
switch select
    case 1
        Y_iter=result_beta(:,Var_select,Sevrty_select);
        plot(X_iter,Y_iter);
    case 2
        Y_iter=result_sigma(:,Sevrty_select,Sevrty_select);
        plot(X_iter,Y_iter);
    case 3
        pcount=1;
        if selection(1) == 0
            vstart = 2;
        else
            vstart = 1;
        end
        for Var_select=vstart:Nvar
            subplot(subpM,subpN,pcount);
            [f1,xi1]=ksdensity(result_beta(:,Var_select,1));
            plot(xi1,f1,'-r','LineWidth',2);
            hold on;
            [f2,xi2]=ksdensity(result_beta(:,Var_select,2));
            plot(xi2,f2,'--b','LineWidth',2);
            hold on;
            [f3,xi3]=ksdensity(result_beta(:,Var_select,3));
            plot(xi3,f3,'-.c','LineWidth',2);
            xlabel('Coefficient Value','fontsize',10.5);
            ylabel('Probability Density','fontsize',10.5);
            title_str=strcat('Variable ',int2str(Var_select));
            title(title_str,'fontsize',12);
            pcount=pcount+1;
        end
        legend('show')
    case 4
        ksdensity(result_sigma(:,Var_select,Sevrty_select));        
end

%% Print
fprintf('Number of Variables:%d\tTotal Number of Iterations:%d\tBurnin:%d\n',Nvar,N,burnin);
fprintf('Sigma:\n')
fprintf('Dependent variable\tMean\tStd\t95_HDR\n');
m_sigma=zeros(Nsev,Nsev);
for i=1:Nsev
    for j=1:Nsev
        fprintf('Severity: %d-%d\t',i,j);
        m=mean(r_sigma(:,i,j));
        sd=std(r_sigma(:,i,j));
        m_sigma(i,j)=m;
        %95% HDR
        [f,xi] = ksdensity(r_sigma(:,i,j));
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
%         t=m/sd;
        fprintf('%d\t%d\t%d\t%d\n',m,sd,lower,upper);
    end
end
m_beta=zeros(Nvar,Nsev);
for i=1:Nsev
    n=0;
    dataSev=zeros(1,Nvar);
    for m=burnin:N
        if result_flag(m,i)>=0%result_flag(m,i)~=0
            n=n+1;
            dataSev(n,:)=result_beta(m,:,i);
        end
    end
    fprintf('Severity: %d:\n',i);
    fprintf('Variables\tMean\tStd\t95_HDR\tElasticity\n');
    for j=1:Nvar
        fprintf('Variable:%d\t',j);
        m=mean(dataSev(:,j));
        sd=std(dataSev(:,j));
        m_beta(j,i)=m;
        %95% HDR
        [f,xi] = ksdensity(dataSev(:,j));
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
%         t=m/sd;
        elasticity=mean(m*X(:,j));
        fprintf('%d\t%d\t%d\t%d\t%d\n',m,sd,lower,upper,elasticity);
    end
end
