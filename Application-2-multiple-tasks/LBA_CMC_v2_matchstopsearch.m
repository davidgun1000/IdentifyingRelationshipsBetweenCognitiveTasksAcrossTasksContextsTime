function [theta_latent,llh]=LBA_CMC_v2_matchstopsearch(data,param,...
theta_latent,...
num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,burn,adapt)
%this is the Conditional Monte Carlo algorithm to sample the random effects
%for each subjects.
parfor j=1:num_subjects
        
    % generating the particles for the random effects and adjust the size
    % of the vector of each random effect particles 
    
    epsilon=0.3; %the scaling parameter for the proposal during burn in and initial adaptation.
    reference_par=theta_latent(j,:); % the set of random effects from previous iteration of MCMC for conditioning.
    num_randeffect=param.num_randeffect;
    % If count>(switch_num=20), then we use the more
    % efficient proposal, otherwise we use the initial proposal. Both proposals are outlined clearly in Gunawan et al (2018).     

    %if sum(count>switch_num)==num_subjects
    if i>=burn+adapt+1
    w1_mix=0.75;
    w2_mix=0.05;
    w3_mix=1-w1_mix-w2_mix;% setting the weight of the mixture for the proposal in the sampling stage. 
    % generating the proposals from the mixture distribution in the sampling stage
    %-----------------------
    u=rand(num_particles,1); 
    id1=(u<w1_mix);
    id2=(u>w1_mix) & (u<=(w1_mix+w2_mix));
    id3=(u>(w1_mix+w2_mix)) & (u<=(w1_mix+w2_mix+w3_mix));
    n1=sum(id1);
    n2=sum(id2);
    n3=sum(id3);
    id1=logical(id1);
    id2=logical(id2);
    id3=logical(id3);
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
    chol_theta_sig2_15=[chol_theta_sig2(15,1:14),log(chol_theta_sig2(15,15))];
    chol_theta_sig2_16=[chol_theta_sig2(16,1:15),log(chol_theta_sig2(16,16))];
    chol_theta_sig2_17=[chol_theta_sig2(17,1:16),log(chol_theta_sig2(17,17))];
    chol_theta_sig2_18=[chol_theta_sig2(18,1:17),log(chol_theta_sig2(18,18))];
    chol_theta_sig2_19=[chol_theta_sig2(19,1:18),log(chol_theta_sig2(19,19))];
    chol_theta_sig2_20=[chol_theta_sig2(20,1:19),log(chol_theta_sig2(20,20))];
    chol_theta_sig2_21=[chol_theta_sig2(21,1:20),log(chol_theta_sig2(21,21))];
    chol_theta_sig2_22=[chol_theta_sig2(22,1:21),log(chol_theta_sig2(22,22))];
    chol_theta_sig2_23=[chol_theta_sig2(23,1:22),log(chol_theta_sig2(23,23))];
    chol_theta_sig2_24=[chol_theta_sig2(24,1:23),log(chol_theta_sig2(24,24))];
    chol_theta_sig2_25=[chol_theta_sig2(25,1:24),log(chol_theta_sig2(25,25))];
    chol_theta_sig2_26=[chol_theta_sig2(26,1:25),log(chol_theta_sig2(26,26))];
    chol_theta_sig2_27=[chol_theta_sig2(27,1:26),log(chol_theta_sig2(27,27))];

        
    xx=[param.theta_mu';chol_theta_sig2_1';chol_theta_sig2_2';chol_theta_sig2_3';...
        chol_theta_sig2_4';chol_theta_sig2_5';chol_theta_sig2_6';chol_theta_sig2_7';chol_theta_sig2_8';
        chol_theta_sig2_9';chol_theta_sig2_10';chol_theta_sig2_11';chol_theta_sig2_12';
        chol_theta_sig2_13';chol_theta_sig2_14';chol_theta_sig2_15';chol_theta_sig2_16';
        chol_theta_sig2_17';chol_theta_sig2_18';chol_theta_sig2_19';chol_theta_sig2_20';
        chol_theta_sig2_21';chol_theta_sig2_22';chol_theta_sig2_23';chol_theta_sig2_24';
        chol_theta_sig2_25';chol_theta_sig2_26';chol_theta_sig2_27'];% we need this to compute the mean of the proposal in the sampling stage
    cond_mean=mean_theta(j,1:num_randeffect)'+covmat_theta(1:num_randeffect,num_randeffect+1:end,j)*((covmat_theta(num_randeffect+1:end,num_randeffect+1:end,j))\(xx-mean_theta(j,num_randeffect+1:end)')); % computing the mean of the proposal in the sampling stage
    %cond_mean_ref=reference_par'+covmat_theta(1:num_randeffect,num_randeffect+1:end,j)*((covmat_theta(num_randeffect+1:end,num_randeffect+1:end,j))\(xx-mean_theta(j,num_randeffect+1:end)'));
    cond_mean_ref=reference_par';
    cond_var=covmat_theta(1:num_randeffect,1:num_randeffect,j)-covmat_theta(1:num_randeffect,num_randeffect+1:end,j)*(covmat_theta(num_randeffect+1:end,num_randeffect+1:end,j)\covmat_theta(num_randeffect+1:end,1:num_randeffect,j)); % computing the variance of the proposal in the sampling stage
    chol_cond_var=chol(cond_var,'lower');
    rnorm1=cond_mean+chol_cond_var*randn(param.num_randeffect,n1);
    chol_covmat=chol(param.theta_sig2,'lower');
    rnorm2=param.theta_mu'+chol_covmat*randn(param.num_randeffect,n2);
    rnorm3=cond_mean_ref+chol_cond_var*randn(param.num_randeffect,n3);
    rnorm=[rnorm1,rnorm2,rnorm3];
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
    rnorm_theta_v210=rnorm(:,6);
    rnorm_theta_v220=rnorm(:,7);
    rnorm_theta_v230=rnorm(:,8);
    rnorm_theta_tau0=rnorm(:,9);
    rnorm_theta_b11=rnorm(:,10);
    rnorm_theta_b21=rnorm(:,11);
    rnorm_theta_b31=rnorm(:,12);
    rnorm_theta_A1=rnorm(:,13);
    rnorm_theta_v11=rnorm(:,14);
    rnorm_theta_v211=rnorm(:,15);
    rnorm_theta_v221=rnorm(:,16);
    rnorm_theta_v231=rnorm(:,17);
    rnorm_theta_tau1=rnorm(:,18);
    rnorm_theta_b12=rnorm(:,19);
    rnorm_theta_b22=rnorm(:,20);
    rnorm_theta_b32=rnorm(:,21);
    rnorm_theta_A2=rnorm(:,22);
    rnorm_theta_v12=rnorm(:,23);
    rnorm_theta_v212=rnorm(:,24);
    rnorm_theta_v222=rnorm(:,25);
    rnorm_theta_v232=rnorm(:,26);
    rnorm_theta_tau2=rnorm(:,27);
    
    % set the first particles to the values of the random effects from the
    % previous iterations of MCMC for conditioning
    rnorm_theta_b10(1,1)=reference_par(1,1);
    rnorm_theta_b20(1,1)=reference_par(1,2);
    rnorm_theta_b30(1,1)=reference_par(1,3);
    rnorm_theta_A0(1,1)=reference_par(1,4);
    rnorm_theta_v10(1,1)=reference_par(1,5);
    rnorm_theta_v210(1,1)=reference_par(1,6);
    rnorm_theta_v220(1,1)=reference_par(1,7);
    rnorm_theta_v230(1,1)=reference_par(1,8);    
    rnorm_theta_tau0(1,1)=reference_par(1,9);
    rnorm_theta_b11(1,1)=reference_par(1,10);
    rnorm_theta_b21(1,1)=reference_par(1,11);
    rnorm_theta_b31(1,1)=reference_par(1,12);
    rnorm_theta_A1(1,1)=reference_par(1,13);
    rnorm_theta_v11(1,1)=reference_par(1,14);
    rnorm_theta_v211(1,1)=reference_par(1,15);
    rnorm_theta_v221(1,1)=reference_par(1,16);
    rnorm_theta_v231(1,1)=reference_par(1,17);
    rnorm_theta_tau1(1,1)=reference_par(1,18);
    rnorm_theta_b12(1,1)=reference_par(1,19);
    rnorm_theta_b22(1,1)=reference_par(1,20);
    rnorm_theta_b32(1,1)=reference_par(1,21);
    rnorm_theta_A2(1,1)=reference_par(1,22);
    rnorm_theta_v12(1,1)=reference_par(1,23);
    rnorm_theta_v212(1,1)=reference_par(1,24);
    rnorm_theta_v222(1,1)=reference_par(1,25);
    rnorm_theta_v232(1,1)=reference_par(1,26);
    rnorm_theta_tau2(1,1)=reference_par(1,27);

    
    rnorm_theta=[rnorm_theta_b10,rnorm_theta_b20,rnorm_theta_b30,rnorm_theta_A0,rnorm_theta_v10,rnorm_theta_v210,rnorm_theta_v220,rnorm_theta_v230,rnorm_theta_tau0,...
        rnorm_theta_b11,rnorm_theta_b21,rnorm_theta_b31,rnorm_theta_A1,rnorm_theta_v11,rnorm_theta_v211,rnorm_theta_v221,rnorm_theta_v231,rnorm_theta_tau1,...
        rnorm_theta_b12,rnorm_theta_b22,rnorm_theta_b32,rnorm_theta_A2,rnorm_theta_v12,rnorm_theta_v212,rnorm_theta_v222,rnorm_theta_v232,rnorm_theta_tau2,...
        ];
    
    %adjust the size of the vectors of the random effects
    
    rnorm_theta_b10_kron=kron(rnorm_theta_b10,ones(num_trials(j,1),1));
    rnorm_theta_b20_kron=kron(rnorm_theta_b20,ones(num_trials(j,1),1));
    rnorm_theta_b30_kron=kron(rnorm_theta_b30,ones(num_trials(j,1),1));
    rnorm_theta_A0_kron=kron(rnorm_theta_A0,ones(num_trials(j,1),1));
    rnorm_theta_v10_kron=kron(rnorm_theta_v10,ones(num_trials(j,1),1));
    rnorm_theta_v210_kron=kron(rnorm_theta_v210,ones(num_trials(j,1),1));
    rnorm_theta_v220_kron=kron(rnorm_theta_v220,ones(num_trials(j,1),1));
    rnorm_theta_v230_kron=kron(rnorm_theta_v230,ones(num_trials(j,1),1));
    rnorm_theta_tau0_kron=kron(rnorm_theta_tau0,ones(num_trials(j,1),1));
    
    rnorm_theta_b11_kron=kron(rnorm_theta_b11,ones(num_trials(j,1),1));
    rnorm_theta_b21_kron=kron(rnorm_theta_b21,ones(num_trials(j,1),1));
    rnorm_theta_b31_kron=kron(rnorm_theta_b31,ones(num_trials(j,1),1));
    rnorm_theta_A1_kron=kron(rnorm_theta_A1,ones(num_trials(j,1),1));
    rnorm_theta_v11_kron=kron(rnorm_theta_v11,ones(num_trials(j,1),1));
    rnorm_theta_v211_kron=kron(rnorm_theta_v211,ones(num_trials(j,1),1));
    rnorm_theta_v221_kron=kron(rnorm_theta_v221,ones(num_trials(j,1),1));
    rnorm_theta_v231_kron=kron(rnorm_theta_v231,ones(num_trials(j,1),1));
    rnorm_theta_tau1_kron=kron(rnorm_theta_tau1,ones(num_trials(j,1),1));
    
    rnorm_theta_b12_kron=kron(rnorm_theta_b12,ones(num_trials(j,1),1));
    rnorm_theta_b22_kron=kron(rnorm_theta_b22,ones(num_trials(j,1),1));
    rnorm_theta_b32_kron=kron(rnorm_theta_b32,ones(num_trials(j,1),1));
    rnorm_theta_A2_kron=kron(rnorm_theta_A2,ones(num_trials(j,1),1));
    rnorm_theta_v12_kron=kron(rnorm_theta_v12,ones(num_trials(j,1),1));
    rnorm_theta_v212_kron=kron(rnorm_theta_v212,ones(num_trials(j,1),1));
    rnorm_theta_v222_kron=kron(rnorm_theta_v222,ones(num_trials(j,1),1));
    rnorm_theta_v232_kron=kron(rnorm_theta_v232,ones(num_trials(j,1),1));
    rnorm_theta_tau2_kron=kron(rnorm_theta_tau2,ones(num_trials(j,1),1));

    
    %adjust the size of the dataset
    
    data_response_repmat=repmat(data.response{j,1}(:,1),num_particles,1);
    data_rt_repmat=repmat(data.rt{j,1}(:,1),num_particles,1);
    data_cond_repmat=repmat(data.cond{j,1}(:,1),num_particles,1);  
    data_task_repmat=repmat(data.task{j,1}(:,1),num_particles,1);
    data_popOut_repmat=repmat(data.popOut{j,1}(:,1),num_particles,1);
    
    [rnorm_theta_b_kron]=reshape_b_matchsearchstop(data_cond_repmat,data_task_repmat,data_popOut_repmat,rnorm_theta_b10_kron,rnorm_theta_b20_kron,rnorm_theta_b30_kron,...
        rnorm_theta_b11_kron,rnorm_theta_b21_kron,rnorm_theta_b31_kron,...
        rnorm_theta_b12_kron,rnorm_theta_b22_kron,rnorm_theta_b32_kron); %choose the threshold particles to match with the conditions of the experiment and whether it is inside or outside scanner
    [rnorm_theta_v_kron]=reshape_v_matchsearchstop(data_response_repmat,data_cond_repmat,data_task_repmat,data_popOut_repmat,rnorm_theta_v10_kron,rnorm_theta_v210_kron,rnorm_theta_v220_kron,rnorm_theta_v230_kron,...
        rnorm_theta_v11_kron,rnorm_theta_v211_kron,rnorm_theta_v221_kron,rnorm_theta_v231_kron,...
        rnorm_theta_v12_kron,rnorm_theta_v212_kron,rnorm_theta_v222_kron,rnorm_theta_v232_kron); %set the drift rate particles to match with the response
    [rnorm_theta_A_kron]=reshape_A_matchsearchstop(data_task_repmat,rnorm_theta_A0_kron,rnorm_theta_A1_kron,rnorm_theta_A2_kron);
    [rnorm_theta_tau_kron]=reshape_tau_matchsearchstop(data_task_repmat,rnorm_theta_tau0_kron,rnorm_theta_tau1_kron,rnorm_theta_tau2_kron);
    
    %computing the log density of the LBA given the particles of random
    %effects
    
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    lw_reshape=reshape(lw,num_trials(j,1),num_particles);
    logw_first=sum(lw_reshape);
    
    %computing the log of p(\alpha|\theta) and density of the proposal for
    %burn in and initial sampling stage (count<=switch_num) and sampling
    %stage (count>switch_num)
    
    %if  sum(count>switch_num)==num_subjects
    if i>=burn+adapt+1
        logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw_third=log(w1_mix.*mvnpdf(rnorm_theta,cond_mean',chol_cond_var*chol_cond_var')+...
            (w2_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat')+...
            w3_mix.*mvnpdf(rnorm_theta,cond_mean_ref',chol_cond_var*chol_cond_var'));
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