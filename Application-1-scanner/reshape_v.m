function [rnorm_v_kron]=reshape_v(data_response_repmat,data_scan_repmat,rnorm_v10_kron,rnorm_v20_kron,rnorm_v11_kron,rnorm_v21_kron)

     ind10=data_response_repmat==1 & data_scan_repmat==0;
     ind20=data_response_repmat==2 & data_scan_repmat==0;
      
     rnorm_v_kron(ind10,:)=[rnorm_v10_kron(ind10,1),rnorm_v20_kron(ind10,1)];
     rnorm_v_kron(ind20,:)=[rnorm_v20_kron(ind20,1),rnorm_v10_kron(ind20,1)];

     ind11=data_response_repmat==1 & data_scan_repmat==1;
     ind21=data_response_repmat==2 & data_scan_repmat==1;
     
     rnorm_v_kron(ind11,:)=[rnorm_v11_kron(ind11,1),rnorm_v21_kron(ind11,1)];
     rnorm_v_kron(ind21,:)=[rnorm_v21_kron(ind21,1),rnorm_v11_kron(ind21,1)];
end