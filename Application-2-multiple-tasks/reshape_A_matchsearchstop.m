function [rnorm_A_kron]=reshape_A_matchsearchstop(data_task_repmat,rnorm_A0_kron,rnorm_A1_kron,rnorm_A2_kron)

    ind1=data_task_repmat==1;
    ind2=data_task_repmat==2;
    ind3=data_task_repmat==3;
    
    rnorm_A_kron(ind1,:)=rnorm_A0_kron(ind1,1);
    rnorm_A_kron(ind2,:)=rnorm_A1_kron(ind2,1);
    rnorm_A_kron(ind3,:)=rnorm_A2_kron(ind3,1);


end