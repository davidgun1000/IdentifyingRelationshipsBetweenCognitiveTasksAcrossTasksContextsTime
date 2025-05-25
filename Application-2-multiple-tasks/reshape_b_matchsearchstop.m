function [rnorm_b_kron]=reshape_b_matchsearchstop(data_cond_repmat,data_task_repmat,data_popOut_repmat,rnorm_b10_kron,rnorm_b20_kron,rnorm_b30_kron,rnorm_b11_kron,rnorm_b21_kron,rnorm_b31_kron,...
    rnorm_b12_kron,rnorm_b22_kron,rnorm_b32_kron)

      %match task-----------------------------------
        
      ind10=data_cond_repmat==1 & data_task_repmat==1;
      rnorm_b_kron(ind10,:)=rnorm_b10_kron(ind10,1);
      
      ind20=data_cond_repmat==2 & data_task_repmat==1;
      rnorm_b_kron(ind20,:)=rnorm_b20_kron(ind20,1);
      
      ind30=data_cond_repmat==3 & data_task_repmat==1;
      rnorm_b_kron(ind30,:)=rnorm_b30_kron(ind30,1);
      
      %search task-----------------------------------
      
      ind11=                      data_task_repmat==2 & data_popOut_repmat==1; 
      rnorm_b_kron(ind11,:)=rnorm_b11_kron(ind11,1);
      
      ind21=data_cond_repmat==2 & data_task_repmat==2 & data_popOut_repmat==2;
      rnorm_b_kron(ind21,:)=rnorm_b21_kron(ind21,1);
      
      ind31=data_cond_repmat==3 & data_task_repmat==2 & data_popOut_repmat==2;
      rnorm_b_kron(ind31,:)=rnorm_b31_kron(ind31,1);
      
      %stop task-------------------------------------
      
      ind12=                      data_task_repmat==3 & data_popOut_repmat==1;
      rnorm_b_kron(ind12,:)=rnorm_b12_kron(ind12,1);

      ind22=data_cond_repmat==2 & data_task_repmat==3 & data_popOut_repmat==2;
      rnorm_b_kron(ind22,:)=rnorm_b22_kron(ind22,1);
      
      ind32=data_cond_repmat==3 & data_task_repmat==3 & data_popOut_repmat==2;
      rnorm_b_kron(ind32,:)=rnorm_b32_kron(ind32,1);
end