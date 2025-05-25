function [theta_latent,llh]=LBA_CMC_v1(data,param,...
theta_latent,...
num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,count,switch_num)
%this is the Conditional Monte Carlo algorithm to sample the random effects
%for each subjects.
parfor j=1:num_subjects
        
    % generating the particles for the random effects and adjust the size
    % of the vector of each random effect particles 
    
    epsilon=0.4; %the scaling parameter for the proposal during burn in and initial adaptation.
    reference_par=theta_latent(j,:); % the set of random effects from previous iteration of MCMC for conditioning.
    num_randeffect=param.num_randeffect;
    % If count>(switch_num=20), then we use the more
    % efficient proposal, otherwise we use the initial proposal. Both proposals are outlined clearly in Gunawan et al (2018).     

    if sum(count>switch_num)==num_subjects
    w_mix=0.9; % setting the weight of the mixture for the proposal in the sampling stage. 
    % generating the proposals from the mixture distribution in the sampling stage
    %-----------------------
    u=rand(num_particles,1); 
    id1=(u<w_mix);
    id2=1-id1;
    n1=sum(id1);
    n2=num_particles-n1;
    chol_theta_sig2=chol(param.theta_sig2,'lower');
    chol_theta_sig2_1=log(chol_theta_sig2(1,1));
    chol_theta_sig2_2=[chol_theta_sig2(2,1),log(chol_theta_sig2(2,2))];
    chol_theta_sig2_3=[chol_theta_sig2(3,1:2),log(chol_theta_sig2(3,3))];
    chol_theta_sig2_4=[chol_theta_sig2(4,1:3),log(chol_theta_sig2(4,4))];
    chol_theta_sig2_5=[chol_theta_sig2(5,1:4),log(chol_theta_sig2(5,5))];
    chol_theta_sig2_6=[chol_theta_sig2(6,1:5),log(chol_theta_sig2(6,6))];
    chol_theta_sig2_7=[chol_theta_sig2(7,1:6),log(chol_theta_sig2(7,7))];
    chol_theta_sig2_8=[chol_theta_sig2(8,1:7),log(chol_theta_sig2(8,8))];
    chol_theta_sig2_9=[chol_theta_sig2(9,1:8),log(chol_theta_sig2(9,9))];
    chol_theta_sig2_10=[chol_theta_sig2(10,1:9),log(chol_theta_sig2(10,10))];
    chol_theta_sig2_11=[chol_theta_sig2(11,1:10),log(chol_theta_sig2(11,11))];
    chol_theta_sig2_12=[chol_theta_sig2(12,1:11),log(chol_theta_sig2(12,12))];
    chol_theta_sig2_13=[chol_theta_sig2(13,1:12),log(chol_theta_sig2(13,13))];
    chol_theta_sig2_14=[chol_theta_sig2(14,1:13),log(chol_theta_sig2(14,14))];
    
        
    xx=[param.theta_mu';chol_theta_sig2_1';chol_theta_sig2_2';chol_theta_sig2_3';...
        chol_theta_sig2_4';chol_theta_sig2_5';chol_theta_sig2_6';chol_theta_sig2_7';chol_theta_sig2_8';
        chol_theta_sig2_9';chol_theta_sig2_10';chol_theta_sig2_11';chol_theta_sig2_12';
        chol_theta_sig2_13';chol_theta_sig2_14'];% we need this to compute the mean of the proposal in the sampling stage
    cond_mean=mean_theta(j,1:num_randeffect)'+covmat_theta(1:num_randeffect,num_randeffect+1:end,j)*((covmat_theta(num_randeffect+1:end,num_randeffect+1:end,j))\(xx-mean_theta(j,num_randeffect+1:end)')); % computing the mean of the proposal in the sampling stage
    cond_var=covmat_theta(1:num_randeffect,1:num_randeffect,j)-covmat_theta(1:num_randeffect,num_randeffect+1:end,j)*(covmat_theta(num_randeffect+1:end,num_randeffect+1:end,j)\covmat_theta(num_randeffect+1:end,1:num_randeffect,j)); % computing the variance of the proposal in the sampling stage
    chol_cond_var=chol(cond_var,'lower');
    rnorm1=cond_mean+chol_cond_var*randn(param.num_randeffect,n1);
    chol_covmat=chol(param.theta_sig2,'lower');
    rnorm2=param.theta_mu'+chol_covmat*randn(param.num_randeffect,n2);
    rnorm=[rnorm1,rnorm2];
    rnorm=rnorm';
    %-----------------------
    else
    % generating the proposals from the mixture distribution in the burn in
    % and initial sampling stage
    %-----------------------
    w_mix=0.95; % setting the weights of the mixture in the burn in and initial sampling stage.
    u=rand(num_particles,1);
    id1=(u<w_mix);
    n1=sum(id1);
    n2=num_particles-n1;
    chol_covmat=chol(param.theta_sig2,'lower');
    rnorm1=reference_par'+epsilon.*chol_covmat*randn(num_randeffect,n1);
    rnorm2=param.theta_mu'+chol_covmat*randn(num_randeffect,n2);
    
    rnorm=[rnorm1,rnorm2];
    rnorm=rnorm';
    
    %------------------------
    end   
    
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
    
    
    % set the first particles to the values of the random effects from the
    % previous iterations of MCMC for conditioning
    rnorm_theta_b10(1,1)=reference_par(1,1);
    rnorm_theta_b20(1,1)=reference_par(1,2);
    rnorm_theta_b30(1,1)=reference_par(1,3);
    rnorm_theta_A0(1,1)=reference_par(1,4);
    rnorm_theta_v10(1,1)=reference_par(1,5);
    rnorm_theta_v20(1,1)=reference_par(1,6);
    rnorm_theta_tau0(1,1)=reference_par(1,7);
    rnorm_theta_b11(1,1)=reference_par(1,8);
    rnorm_theta_b21(1,1)=reference_par(1,9);
    rnorm_theta_b31(1,1)=reference_par(1,10);
    rnorm_theta_A1(1,1)=reference_par(1,11);
    rnorm_theta_v11(1,1)=reference_par(1,12);
    rnorm_theta_v21(1,1)=reference_par(1,13);
    rnorm_theta_tau1(1,1)=reference_par(1,14);

    
    rnorm_theta=[rnorm_theta_b10,rnorm_theta_b20,rnorm_theta_b30,rnorm_theta_A0,rnorm_theta_v10,rnorm_theta_v20,rnorm_theta_tau0,...
        rnorm_theta_b11,rnorm_theta_b21,rnorm_theta_b31,rnorm_theta_A1,rnorm_theta_v11,rnorm_theta_v21,rnorm_theta_tau1];
    
    %adjust the size of the vectors of the random effects
    
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
        rnorm_theta_b11_kron,rnorm_theta_b21_kron,rnorm_theta_b31_kron); %choose the threshold particles to match with the conditions of the experiment
    [rnorm_theta_v_kron]=reshape_v(data_response_repmat,data_scan_repmat,rnorm_theta_v10_kron,rnorm_theta_v20_kron,rnorm_theta_v11_kron,rnorm_theta_v21_kron); %set the drift rate particles to match with the response
    [rnorm_theta_A_kron]=reshape_A(data_scan_repmat,rnorm_theta_A0_kron,rnorm_theta_A1_kron);
    [rnorm_theta_tau_kron]=reshape_tau(data_scan_repmat,rnorm_theta_tau0_kron,rnorm_theta_tau1_kron);
    
    %computing the log density of the LBA given the particles of random
    %effects
    
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    lw_reshape=reshape(lw,num_trials(j,1),num_particles);
    logw_first=sum(lw_reshape);
    
    %computing the log of p(\alpha|\theta) and density of the proposal for
    %burn in and initial sampling stage (count<=switch_num) and sampling
    %stage (count>switch_num)
    
    if  sum(count>switch_num)==num_subjects
        logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw_third=log(w_mix.*mvnpdf(rnorm_theta,cond_mean',chol_cond_var*chol_cond_var')+...
            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw=logw_first'+logw_second'-logw_third;
    else
        logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw_third=log(w_mix.*mvnpdf(rnorm_theta,reference_par,(epsilon^2).*(chol_covmat*chol_covmat'))+...
            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw=logw_first'+logw_second'-logw_third;
    end
    
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
    llh_i(j) = max_logw+log(mean(weight)); 
    llh_i(j) = real(llh_i(j)); 	
    weight=weight./sum(weight);
    if sum(weight<0)>0
        id=weight<0;
        id=1-id;
        id=logical(id);
        weight=weight(id,1);
    end
    Nw=length(weight);
    
    if Nw>0 
        ind=randsample(Nw,1,true,weight);
        theta_latent(j,:)=rnorm_theta(ind,:);
    end        
%----------------------------------------------------------------------------------------------------------------------------------    
    
end
llh=sum(llh_i);
end


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