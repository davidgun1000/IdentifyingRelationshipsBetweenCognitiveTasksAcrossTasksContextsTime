function [rnorm_tau_kron]=reshape_tau_matchsearchstop(data_task_repmat,rnorm_tau0_kron,rnorm_tau1_kron,rnorm_tau2_kron)

    ind1=data_task_repmat==1;
    ind2=data_task_repmat==2;
    ind3=data_task_repmat==3;
    
    rnorm_tau_kron(ind1,:)=rnorm_tau0_kron(ind1,1);
    rnorm_tau_kron(ind2,:)=rnorm_tau1_kron(ind2,1);
    rnorm_tau_kron(ind3,:)=rnorm_tau2_kron(ind3,1);

end