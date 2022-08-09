% Author: Sobhan Goudarzi
% IMage Processing And Characterization of Tissue (IMPACT) Group
% Concordia University
% email address: sobhan.goudarzi@concordia.ca
% August 2022
function [u, v, lambd, objective] = ADMM_v0(Phi, y, u0, v0, lambd0, eps, mu, beta, maxiter)
global y v_k lambd_k beta Phi
u_k = u0;
v_k = v0;
lambd_k = lambd0;
objective = 1e6;
options = [];
options.display = 'none';
options.maxFunEvals = 100;
options.Method = 'lbfgs';
for i=1:maxiter
    %%%%%%%%%%%%%%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_u  1/2*||y-u||^2 + beta/2*||u - v + lambd/beta||^2
    u_k = minFunc(@funcu,zeros(484137,1),options);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update v %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_v  mu||v||_1 + beta/2*||u - v + lambd/beta||^2
    v_k = max(abs(u_k+lambd_k/beta)-mu/beta,0).*sign(u_k+lambd_k/beta);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update lambd %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambd_k = lambd_k+beta*(u_k-v_k);
    %%%%%%%%%%%%%%%%%%%% stopping criterion %%%%%%%%%%%%%%%%%%%%%%%
    objective(i+1) = 0.5*norm(y-full(Phi*sparse(u_k)))^2 + mu * sum(abs(v_k))...
        +(beta/2)*norm(u_k-v_k+(lambd_k/beta))^2;
    criterion(i) = abs(objective(i+1)-objective(i))/objective(i);
    disp(['iteration ',num2str(i),', criterion: ',num2str(criterion(i))])
    if ( criterion(i) < eps )
        u(:,i) = u_k;
        v(:,i) = v_k;
        lambd(:,i) = lambd_k;
        break
    end
    u(:,i) = u_k;
    v(:,i) = v_k;
    lambd(:,i) = lambd_k;
end
objective = objective(2:end);
end

