function [rnorm_b_kron]=reshape_b(data_cond_repmat,data_scan_repmat,rnorm_b10_kron,rnorm_b20_kron,rnorm_b30_kron,rnorm_b11_kron,rnorm_b21_kron,rnorm_b31_kron)

      ind10=data_cond_repmat==1 & data_scan_repmat==0;
      rnorm_b_kron(ind10,:)=rnorm_b10_kron(ind10,1);
      
      ind20=data_cond_repmat==2 & data_scan_repmat==0;
      rnorm_b_kron(ind20,:)=rnorm_b20_kron(ind20,1);
      
      ind30=data_cond_repmat==3 & data_scan_repmat==0;
      rnorm_b_kron(ind30,:)=rnorm_b30_kron(ind30,1);
      
      ind11=data_cond_repmat==1 & data_scan_repmat==1;
      rnorm_b_kron(ind11,:)=rnorm_b11_kron(ind11,1);
      
      ind21=data_cond_repmat==2 & data_scan_repmat==1;
      rnorm_b_kron(ind21,:)=rnorm_b21_kron(ind21,1);
      
      ind31=data_cond_repmat==3 & data_scan_repmat==1;
      rnorm_b_kron(ind31,:)=rnorm_b31_kron(ind31,1);

end