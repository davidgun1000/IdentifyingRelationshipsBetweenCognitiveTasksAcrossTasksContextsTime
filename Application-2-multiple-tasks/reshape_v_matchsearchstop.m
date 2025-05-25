function [rnorm_v_kron]=reshape_v_matchsearchstop(data_response_repmat,data_cond_repmat,data_task_repmat,data_popOut_repmat,rnorm_v10_kron,rnorm_v210_kron,rnorm_v220_kron,rnorm_v230_kron,...
    rnorm_v11_kron,rnorm_v211_kron,rnorm_v221_kron,rnorm_v231_kron,...
    rnorm_v12_kron,rnorm_v212_kron,rnorm_v222_kron,rnorm_v232_kron)

       %---------------------------match task
       ind011=data_response_repmat==0  & data_cond_repmat==1 & data_task_repmat==1;
       rnorm_v_kron(ind011,:)=[rnorm_v10_kron(ind011,1),rnorm_v210_kron(ind011,1)];

       ind021=data_response_repmat==0 & data_cond_repmat==2 & data_task_repmat==1;
       rnorm_v_kron(ind021,:)=[rnorm_v10_kron(ind021,1),rnorm_v220_kron(ind021,1)];
       
       ind031=data_response_repmat==0 & data_cond_repmat==3 & data_task_repmat==1;
       rnorm_v_kron(ind031,:)=[rnorm_v10_kron(ind031,1),rnorm_v230_kron(ind031,1)];
       
       ind111=data_response_repmat==1  & data_cond_repmat==1 & data_task_repmat==1;
       rnorm_v_kron(ind111,:)=[rnorm_v210_kron(ind111,1),rnorm_v10_kron(ind111,1)];

       ind121=data_response_repmat==1  & data_cond_repmat==2 & data_task_repmat==1;
       rnorm_v_kron(ind121,:)=[rnorm_v220_kron(ind121,1),rnorm_v10_kron(ind121,1)];

       ind131=data_response_repmat==1  & data_cond_repmat==3 & data_task_repmat==1;
       rnorm_v_kron(ind131,:)=[rnorm_v230_kron(ind131,1),rnorm_v10_kron(ind131,1)];

       %-----------------------------search task
       
       ind012=data_response_repmat==0  &                       data_task_repmat==2 & data_popOut_repmat==1;
       rnorm_v_kron(ind012,:)=[rnorm_v11_kron(ind012,1),rnorm_v211_kron(ind012,1)];  

       ind022=data_response_repmat==0  & data_cond_repmat==2 & data_task_repmat==2 & data_popOut_repmat==2;
       rnorm_v_kron(ind022,:)=[rnorm_v11_kron(ind022,1),rnorm_v221_kron(ind022,1)];  

       ind032=data_response_repmat==0  & data_cond_repmat==3 & data_task_repmat==2 & data_popOut_repmat==2;
       rnorm_v_kron(ind032,:)=[rnorm_v11_kron(ind032,1),rnorm_v231_kron(ind032,1)];
       
       ind112=data_response_repmat==1  &                       data_task_repmat==2 & data_popOut_repmat==1;
       rnorm_v_kron(ind112,:)=[rnorm_v211_kron(ind112,1),rnorm_v11_kron(ind112,1)];

       ind122=data_response_repmat==1  & data_cond_repmat==2 & data_task_repmat==2 & data_popOut_repmat==2;
       rnorm_v_kron(ind122,:)=[rnorm_v221_kron(ind122,1),rnorm_v11_kron(ind122,1)];
       
       ind132=data_response_repmat==1  & data_cond_repmat==3 & data_task_repmat==2 & data_popOut_repmat==2;
       rnorm_v_kron(ind132,:)=[rnorm_v231_kron(ind132,1),rnorm_v11_kron(ind132,1)];

       %-------------------------stop task
       ind013=data_response_repmat==0  &                       data_task_repmat==3 &  data_popOut_repmat==1;
       rnorm_v_kron(ind013,:)=[rnorm_v12_kron(ind013,1),rnorm_v212_kron(ind013,1)];
       
       ind023=data_response_repmat==0  & data_cond_repmat==2 & data_task_repmat==3 &  data_popOut_repmat==2;
       rnorm_v_kron(ind023,:)=[rnorm_v12_kron(ind023,1),rnorm_v222_kron(ind023,1)];
       
       ind033=data_response_repmat==0  & data_cond_repmat==3 & data_task_repmat==3 &  data_popOut_repmat==2;
       rnorm_v_kron(ind033,:)=[rnorm_v12_kron(ind033,1),rnorm_v232_kron(ind033,1)];
       
       ind113=data_response_repmat==1  &                       data_task_repmat==3 &  data_popOut_repmat==1;
       rnorm_v_kron(ind113,:)=[rnorm_v212_kron(ind113,1),rnorm_v12_kron(ind113,1)];
       
       ind123=data_response_repmat==1  & data_cond_repmat==2 & data_task_repmat==3 &  data_popOut_repmat==2;
       rnorm_v_kron(ind123,:)=[rnorm_v222_kron(ind123,1),rnorm_v12_kron(ind123,1)];
       
       ind133=data_response_repmat==1  & data_cond_repmat==3 & data_task_repmat==3 &  data_popOut_repmat==2;
       rnorm_v_kron(ind133,:)=[rnorm_v232_kron(ind133,1),rnorm_v12_kron(ind133,1)];
       
end