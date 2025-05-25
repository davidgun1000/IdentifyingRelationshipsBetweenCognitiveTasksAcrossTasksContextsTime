function [theta_latent]=LBA_MC_v1(data,param,num_subjects,num_trials,num_particles)
% This is the Monte Carlo algorithm for generating initial random effects
% values, this is an optional code. We can get the initial values randomly,
% or some other methods.

parfor j=1:num_subjects
    % generating the particles for the random effects and adjust the size
    % of the vector of each random effect particles 
    rnorm=mvnrnd(param.theta_mu',param.theta_sig2,num_particles);
    rnorm_theta_b10=rnorm(:,1); 
    rnorm_theta_b20=rnorm(:,2);
    rnorm_theta_b30=rnorm(:,3);
    rnorm_theta_A0=rnorm(:,4);
    rnorm_theta_v10=rnorm(:,5);
    rnorm_theta_v20=rnorm(:,6);
    rnorm_theta_tau0=rnorm(:,7);
    rnorm_theta_b11=rnorm(:,8); 
    rnorm_theta_b21=rnorm(:,9);
    rnorm_theta_b31=rnorm(:,10);
    rnorm_theta_A1=rnorm(:,11);
    rnorm_theta_v11=rnorm(:,12);
    rnorm_theta_v21=rnorm(:,13);
    rnorm_theta_tau1=rnorm(:,14);
    
    rnorm_theta_b10_kron=kron(rnorm_theta_b10,ones(num_trials(j,1),1));
    rnorm_theta_b20_kron=kron(rnorm_theta_b20,ones(num_trials(j,1),1));
    rnorm_theta_b30_kron=kron(rnorm_theta_b30,ones(num_trials(j,1),1));
    rnorm_theta_A0_kron=kron(rnorm_theta_A0,ones(num_trials(j,1),1));
    rnorm_theta_v10_kron=kron(rnorm_theta_v10,ones(num_trials(j,1),1));
    rnorm_theta_v20_kron=kron(rnorm_theta_v20,ones(num_trials(j,1),1));
    rnorm_theta_tau0_kron=kron(rnorm_theta_tau0,ones(num_trials(j,1),1));
    
    rnorm_theta_b11_kron=kron(rnorm_theta_b11,ones(num_trials(j,1),1));
    rnorm_theta_b21_kron=kron(rnorm_theta_b21,ones(num_trials(j,1),1));
    rnorm_theta_b31_kron=kron(rnorm_theta_b31,ones(num_trials(j,1),1));
    rnorm_theta_A1_kron=kron(rnorm_theta_A1,ones(num_trials(j,1),1));
    rnorm_theta_v11_kron=kron(rnorm_theta_v11,ones(num_trials(j,1),1));
    rnorm_theta_v21_kron=kron(rnorm_theta_v21,ones(num_trials(j,1),1));
    rnorm_theta_tau1_kron=kron(rnorm_theta_tau1,ones(num_trials(j,1),1));
    
    %adjust the size of the dataset 
    
    data_response_repmat=repmat(data.response{j,1}(:,1),num_particles,1);
    data_rt_repmat=repmat(data.rt{j,1}(:,1),num_particles,1);
    data_cond_repmat=repmat(data.cond{j,1}(:,1),num_particles,1);  
    data_scan_repmat=repmat(data.scan{j,1}(:,1),num_particles,1);
    
    [rnorm_theta_b_kron]=reshape_b(data_cond_repmat,data_scan_repmat,rnorm_theta_b10_kron,rnorm_theta_b20_kron,rnorm_theta_b30_kron,...
        rnorm_theta_b11_kron,rnorm_theta_b21_kron,rnorm_theta_b31_kron); %choose the threshold particles to match with the conditions of the experiment and whether it is inside or outside scanner
    [rnorm_theta_v_kron]=reshape_v(data_response_repmat,data_scan_repmat,rnorm_theta_v10_kron,rnorm_theta_v20_kron,rnorm_theta_v11_kron,rnorm_theta_v21_kron); %set the drift rate particles to match with the response
    [rnorm_theta_A_kron]=reshape_A(data_scan_repmat,rnorm_theta_A0_kron,rnorm_theta_A1_kron);
     [rnorm_theta_tau_kron]=reshape_tau(data_scan_repmat,rnorm_theta_tau0_kron,rnorm_theta_tau1_kron);
    %computing the log density of the LBA given the particles of random
    %effects
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    lw_reshape=reshape(lw,num_trials(j,1),num_particles);
    logw_first=sum(lw_reshape);
    logw=logw_first';
    
    %check if there is imaginary number of logw
    id=imag(logw)~=0;
    id=1-id;
    id=logical(id);
    logw=logw(id,1); 
    logw=real(logw);

    if sum(isinf(logw))>0 | sum(isnan(logw))>0
     id=isinf(logw) | isnan(logw);
     id=1-id;
     id=logical(id);
     logw=logw(id,1);
    end
    
    max_logw=max(real(logw));
    weight=real(exp(logw-max_logw));
    weight=weight./sum(weight);
    
    if sum(weight<0)>0
        id=weight<0;
        id=1-id;
        id=logical(id);
        weight=weight(id,1);
    end
    Nw=length(weight);
    rnorm_theta=[rnorm_theta_b10,rnorm_theta_b20,rnorm_theta_b30,rnorm_theta_A0,rnorm_theta_v10,rnorm_theta_v20,rnorm_theta_tau0,...
        rnorm_theta_b11,rnorm_theta_b21,rnorm_theta_b31,rnorm_theta_A1,rnorm_theta_v11,rnorm_theta_v21,rnorm_theta_tau1];
    %choose the random effects according to its weights
    if Nw>0 
        ind=randsample(Nw,1,true,weight);
        theta_latent(j,:)=rnorm_theta(ind,:);
       
    end
        
%----------------------------------------------------------------------------------------------------------------------------------    
    
end

end

 %theta_latent_b1(j,1)=rnorm_theta_b1(ind,1);
        %theta_latent_b2(j,1)=rnorm_theta_b2(ind,1);
        %theta_latent_b3(j,1)=rnorm_theta_b3(ind,1);
        %theta_latent_A(j,1)=rnorm_theta_A(ind,1);
        %theta_latent_v1(j,1)=rnorm_theta_v1(ind,1);
        %theta_latent_v2(j,1)=rnorm_theta_v2(ind,1);
        %theta_latent_tau(j,1)=rnorm_theta_tau(ind,1);
%       rnorm_b=normt_rnd(theta.b_mu,theta.b_sig2,theta.left_truncation,theta.right_truncation,num_particles);
%       rnorm_A=normt_rnd(theta.A_mu,theta.A_sig2,theta.left_truncation,theta.right_truncation,num_particles);
%       rnorm_v1=normt_rnd(theta.v1_mu,theta.v1_sig2,theta.left_truncation,theta.right_truncation,num_particles);
%       rnorm_v2=normt_rnd(theta.v2_mu,theta.v2_sig2,theta.left_truncation,theta.right_truncation,num_particles);
%       
%       rnorm_b(1,1)=latent_b(j,1);
%       rnorm_A(1,1)=latent_A(j,1);
%       rnorm_v1(1,1)=latent_v1(j,1);
%       rnorm_v2(1,1)=latent_v2(j,1);
%       
%       data_response_j=data.response{j,1};
%       data_rt_j=data.rt{j,1};
%       for k=1:num_particles
%       for i=1:num_trials
%           if data_response_j(i,1)==1
%              rnorm_v=[rnorm_v1(k,1),rnorm_v2(k,1)];
%           else
%              rnorm_v=[rnorm_v2(k,1),rnorm_v1(k,1)];
%           end
%           log_lik_ind(i,1)=real(log(LBA_n1PDF(data_rt_j(i,1),rnorm_A(k,1),rnorm_b(k,1),rnorm_v,1)));
%       end
%       lw(k,1)=sum(log_lik_ind);
%       end
%       max_lw=max(lw);
%       weight=real(exp(lw-max_lw));
%       weight=weight./sum(weight);
%       ind=randsample(num_particles,1,true,weight);
%       latent_b(j,1)=rnorm_b(ind,1);
%       latent_A(j,1)=rnorm_A(ind,1);
%       latent_v1(j,1)=rnorm_v1(ind,1);
%       latent_v2(j,1)=rnorm_v2(ind,1);