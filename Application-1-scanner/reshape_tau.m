function [rnorm_tau_kron]=reshape_tau(data_scan_repmat,rnorm_tau0_kron,rnorm_tau1_kron)

    ind0=data_scan_repmat==0;
    ind1=data_scan_repmat==1;
    
    rnorm_tau_kron(ind0,:)=rnorm_tau0_kron(ind0,1);
    rnorm_tau_kron(ind1,:)=rnorm_tau1_kron(ind1,1);



end