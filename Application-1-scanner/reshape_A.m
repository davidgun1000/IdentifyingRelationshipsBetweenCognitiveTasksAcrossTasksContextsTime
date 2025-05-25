function [rnorm_A_kron]=reshape_A(data_scan_repmat,rnorm_A0_kron,rnorm_A1_kron)

    ind0=data_scan_repmat==0;
    ind1=data_scan_repmat==1;
    
    rnorm_A_kron(ind0,:)=rnorm_A0_kron(ind0,1);
    rnorm_A_kron(ind1,:)=rnorm_A1_kron(ind1,1);



end