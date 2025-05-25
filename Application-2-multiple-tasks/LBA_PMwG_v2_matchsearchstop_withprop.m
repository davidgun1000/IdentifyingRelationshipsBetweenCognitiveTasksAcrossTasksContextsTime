% estimating the hierarchical LBA model using PMwG method for the Forstmann (2008) in and out scanner dataset
% The LBA specification can be found in the paper
% The data is stored in the matlab file 'LBA_Forstmann_scanner.mat', it has three
% components: 
% data.cond: the conditions of the experiments, we have three conditions in
% the Forstmann data
% data.rt: the response time
% data.response: response = 1 for incorrect and response = 2 for correct.
% data.scan: scan=0 for out of scanner and scan=1 for in the scanner

load('datamatchsearchstop.mat'); %load the dataset, see an example in the 'LBA_Forstmann_scanner.mat' 
load('prop_matchsearchstop.mat');
num_subjects=length(data.rt); %number of subjects in the experiments
for j=1:num_subjects
    num_trials(j,1)=length(data.rt{j,1}); %computing the number of trials per subject
end
num_particles=500; %number of particles in the conditional Monte Carlo algorithm
parpool(28) %number of processors available to be used.

%initial values of the hyperparameters for lower level parameters
num_randeffect=27; %the total number of random effects in the LBA model. For Forstmann dataset, we have 7 random effects.
prior_mu_mean=zeros(num_randeffect,1); %the prior for \mu_{\alpha}
prior_mu_sig2=eye(num_randeffect);%the prior for \Sigma_{\alpha}
v_half=2; %the hyperparameters of the prior of \Sigma_{\alpha}
A_half=1; %the hyperparameters of the prior of \Sigma_{\alpha}

param.theta_mu=[0.201;0.249;-0.069;-0.373;0.091;1.058;1.182;1.258;-1.761;
    0.412;0.507;0.120;-0.171;0.371;1.150;1.250;1.350;-1.670;
    0.43;0.51;0.13;-0.123;0.34;1.12;1.23;1.31;-1.65]; %the initial values for parameter \mu
param.theta_sig2=iwishrnd(eye(num_randeffect),30); % the initial values for \Sigma
param.sv=1; %we assume that the random effect for the variance of the drift rate is set to 1.
param.num_randeffect=27; %the total number of random effects in the LBA model. For Forstmann dataset, we have 7 random effects.

num_choice=2;% the number of choice
burn=1000;% the burn in iterations
adapt=5000;
nit=15000; % we take the last 10000 draws out of 12000 draws
s=burn+adapt+nit;% the maximum total number of iterations

[theta_latent]=LBA_MC_v1_matchstopsearch(data,param,num_subjects,num_trials,num_particles); %obtain initial values of the random effects.




tot_param=432; %total number of parameters for each subjects
a_half=1./random('gam',1/2,1,num_randeffect,1); %initial values for a_{1},...,a_{7}
count=zeros(1,num_subjects);
%mean_theta=zeros(num_subjects,tot_param); %allocation for proposal mean for the random effects
%covmat_theta=zeros(tot_param,tot_param,num_subjects); %allocation for proposal covariance matrix for the random effects 
switch_num=300;
temp=1;
i=1;

while i<=s
    i
    tic
    theta_latent(:,1)
    %sample parameter \mu_{\alpha} in Gibbs step, look at the paper for the
    %full conditional distribution of the \mu_{\alpha} given rest
    
    var_mu=inv(num_subjects*inv(param.theta_sig2)+inv(prior_mu_sig2));
    mean_mu=var_mu*(inv(param.theta_sig2)*(sum(theta_latent)'));
    chol_var_mu=chol(var_mu,'lower');
    param.theta_mu=mvnrnd(mean_mu,chol_var_mu*chol_var_mu');
    
    %sample parameter \Sigma_{\alpha} in Gibbs step, look at the paper for
    %the full conditional distribution of the \Sigma_{\alpha} given rest
    k_half=v_half+num_randeffect-1+num_subjects;
    cov_temp=zeros(num_randeffect,num_randeffect);
    for j=1:num_subjects
        theta_j=theta_latent(j,:)';
        cov_temp=cov_temp+(theta_j-param.theta_mu')*(theta_j-param.theta_mu')';
    end
    B_half=2*v_half*diag([1./a_half])+cov_temp;
    param.theta_sig2=iwishrnd(B_half,k_half);    
    theta_sig2_inv=inv(param.theta_sig2);

    %sample a_{1},...,a_{7}, look at the paper for the full conditional
    %distribution
    for j=1:num_randeffect
        temp_v_half=(v_half+num_randeffect)/2;
        temp_s_half=(v_half*theta_sig2_inv(j,j)+A_half);
        a_half(j,1)=1./random('gam',temp_v_half,1/temp_s_half);
    end
 
   % training the proposals for conditional Monte
   % Carlo algorithm
   %if sum(count>switch_num)==num_subjects
    if i>=burn+adapt & mod(i,100)==0
        length_draws=length(theta_mu_store);
        for j=1:num_subjects 
             theta=[theta_latent_b10_store(length_draws-(adapt-1):length_draws,j),theta_latent_b20_store(length_draws-(adapt-1):length_draws,j),theta_latent_b30_store(length_draws-(adapt-1):length_draws,j),theta_latent_A0_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_v10_store(length_draws-(adapt-1):length_draws,j),theta_latent_v210_store(length_draws-(adapt-1):length_draws,j),theta_latent_v220_store(length_draws-(adapt-1):length_draws,j),theta_latent_v230_store(length_draws-(adapt-1):length_draws,j),theta_latent_tau0_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_b11_store(length_draws-(adapt-1):length_draws,j),theta_latent_b21_store(length_draws-(adapt-1):length_draws,j),theta_latent_b31_store(length_draws-(adapt-1):length_draws,j),theta_latent_A1_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_v11_store(length_draws-(adapt-1):length_draws,j),theta_latent_v211_store(length_draws-(adapt-1):length_draws,j),theta_latent_v221_store(length_draws-(adapt-1):length_draws,j),theta_latent_v231_store(length_draws-(adapt-1):length_draws,j),theta_latent_tau1_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_b12_store(length_draws-(adapt-1):length_draws,j),theta_latent_b22_store(length_draws-(adapt-1):length_draws,j),theta_latent_b32_store(length_draws-(adapt-1):length_draws,j),theta_latent_A2_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_v12_store(length_draws-(adapt-1):length_draws,j),theta_latent_v212_store(length_draws-(adapt-1):length_draws,j),theta_latent_v222_store(length_draws-(adapt-1):length_draws,j),...
                 theta_latent_v232_store(length_draws-(adapt-1):length_draws,j),theta_latent_tau2_store(length_draws-(adapt-1):length_draws,j),...
                 theta_mu_store(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store1(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store2(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store3(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store4(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store5(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store6(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store7(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store8(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store9(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store10(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store11(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store12(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store13(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store14(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store15(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store16(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store17(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store18(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store19(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store20(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store21(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store22(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store23(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store24(length_draws-(adapt-1):length_draws,:),...
                 chol_theta_sig2_store25(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store26(length_draws-(adapt-1):length_draws,:),chol_theta_sig2_store27(length_draws-(adapt-1):length_draws,:)];
             % in the matrix called theta above, you have to list 
             % (1) all your random effects in the LBA model,(2) followed by the parameters \mu_{\alpha}, and
             % cholesky factor (lower triangular matrix) of the covariance
             % matrix \Sigma_{\alpha}
             covmat_theta(:,:,j)=cov(theta); %computing sample covariance matrix for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
             covmat_theta(:,:,j)=topdm(covmat_theta(:,:,j)); %little correction if the covariance matrix for the proposal is not positive definite matrix.
             mean_theta(j,:)=mean(theta); %computing the sample mean for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
         end           
   end
   %if there are at least 250 unique values of the random effects for each
   %subject, we then switch to better proposal defined in the paper, and
   %reduce the number of particles.
   %if sum(count>switch_num)==num_subjects
   if i>=burn+adapt
      num_particles=500; 
      %t(temp,1)=i;
      %s=t(1,1)+nit;
      %temp=temp+1;
   end
   
   %conditional Monte Carlo algorithm to update the random effects   
   [theta_latent]=LBA_CMC_v2_matchstopsearch_withprop(data,param,theta_latent,num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,burn,adapt);
   
    %storing the MCMC draws
    
    % storing the cholesky factor of the covariance matrix \sigma_{\alpha}
    chol_theta_sig2=chol(param.theta_sig2,'lower');
    chol_theta_sig2_store1(i,:)=log(chol_theta_sig2(1,1));
    chol_theta_sig2_store2(i,:)=[chol_theta_sig2(2,1),log(chol_theta_sig2(2,2))];
    chol_theta_sig2_store3(i,:)=[chol_theta_sig2(3,1:2),log(chol_theta_sig2(3,3))];
    chol_theta_sig2_store4(i,:)=[chol_theta_sig2(4,1:3),log(chol_theta_sig2(4,4))];
    chol_theta_sig2_store5(i,:)=[chol_theta_sig2(5,1:4),log(chol_theta_sig2(5,5))];
    chol_theta_sig2_store6(i,:)=[chol_theta_sig2(6,1:5),log(chol_theta_sig2(6,6))];
    chol_theta_sig2_store7(i,:)=[chol_theta_sig2(7,1:6),log(chol_theta_sig2(7,7))];
    chol_theta_sig2_store8(i,:)=[chol_theta_sig2(8,1:7),log(chol_theta_sig2(8,8))];
    chol_theta_sig2_store9(i,:)=[chol_theta_sig2(9,1:8),log(chol_theta_sig2(9,9))];
    chol_theta_sig2_store10(i,:)=[chol_theta_sig2(10,1:9),log(chol_theta_sig2(10,10))];
    chol_theta_sig2_store11(i,:)=[chol_theta_sig2(11,1:10),log(chol_theta_sig2(11,11))];
    chol_theta_sig2_store12(i,:)=[chol_theta_sig2(12,1:11),log(chol_theta_sig2(12,12))];
    chol_theta_sig2_store13(i,:)=[chol_theta_sig2(13,1:12),log(chol_theta_sig2(13,13))];
    chol_theta_sig2_store14(i,:)=[chol_theta_sig2(14,1:13),log(chol_theta_sig2(14,14))];
    chol_theta_sig2_store15(i,:)=[chol_theta_sig2(15,1:14),log(chol_theta_sig2(15,15))];
    chol_theta_sig2_store16(i,:)=[chol_theta_sig2(16,1:15),log(chol_theta_sig2(16,16))];
    chol_theta_sig2_store17(i,:)=[chol_theta_sig2(17,1:16),log(chol_theta_sig2(17,17))];
    chol_theta_sig2_store18(i,:)=[chol_theta_sig2(18,1:17),log(chol_theta_sig2(18,18))];
    chol_theta_sig2_store19(i,:)=[chol_theta_sig2(19,1:18),log(chol_theta_sig2(19,19))];
    chol_theta_sig2_store20(i,:)=[chol_theta_sig2(20,1:19),log(chol_theta_sig2(20,20))];
    chol_theta_sig2_store21(i,:)=[chol_theta_sig2(21,1:20),log(chol_theta_sig2(21,21))];
    chol_theta_sig2_store22(i,:)=[chol_theta_sig2(22,1:21),log(chol_theta_sig2(22,22))];
    chol_theta_sig2_store23(i,:)=[chol_theta_sig2(23,1:22),log(chol_theta_sig2(23,23))];
    chol_theta_sig2_store24(i,:)=[chol_theta_sig2(24,1:23),log(chol_theta_sig2(24,24))];
    chol_theta_sig2_store25(i,:)=[chol_theta_sig2(25,1:24),log(chol_theta_sig2(25,25))];
    chol_theta_sig2_store26(i,:)=[chol_theta_sig2(26,1:25),log(chol_theta_sig2(26,26))];
    chol_theta_sig2_store27(i,:)=[chol_theta_sig2(27,1:26),log(chol_theta_sig2(27,27))];
    
    
    
    theta_mu_store(i,:)=param.theta_mu';
    theta_sig2_store1(i,:)=param.theta_sig2(1,:);
    theta_sig2_store2(i,:)=param.theta_sig2(2,2:end);
    theta_sig2_store3(i,:)=param.theta_sig2(3,3:end);
    theta_sig2_store4(i,:)=param.theta_sig2(4,4:end);
    theta_sig2_store5(i,:)=param.theta_sig2(5,5:end);
    theta_sig2_store6(i,:)=param.theta_sig2(6,6:end);
    theta_sig2_store7(i,:)=param.theta_sig2(7,7:end);
    theta_sig2_store8(i,:)=param.theta_sig2(8,8:end);
    theta_sig2_store9(i,:)=param.theta_sig2(9,9:end);
    theta_sig2_store10(i,:)=param.theta_sig2(10,10:end);
    theta_sig2_store11(i,:)=param.theta_sig2(11,11:end);
    theta_sig2_store12(i,:)=param.theta_sig2(12,12:end);
    theta_sig2_store13(i,:)=param.theta_sig2(13,13:end);
    theta_sig2_store14(i,:)=param.theta_sig2(14,14:end);
    theta_sig2_store15(i,:)=param.theta_sig2(15,15:end);
    theta_sig2_store16(i,:)=param.theta_sig2(16,16:end);
    theta_sig2_store17(i,:)=param.theta_sig2(17,17:end);
    theta_sig2_store18(i,:)=param.theta_sig2(18,18:end);
    theta_sig2_store19(i,:)=param.theta_sig2(19,19:end);
    theta_sig2_store20(i,:)=param.theta_sig2(20,20:end);
    theta_sig2_store21(i,:)=param.theta_sig2(21,21:end);
    theta_sig2_store22(i,:)=param.theta_sig2(22,22:end);
    theta_sig2_store23(i,:)=param.theta_sig2(23,23:end);
    theta_sig2_store24(i,:)=param.theta_sig2(24,24:end);
    theta_sig2_store25(i,:)=param.theta_sig2(25,25:end);
    theta_sig2_store26(i,:)=param.theta_sig2(26,26:end);
    theta_sig2_store27(i,:)=param.theta_sig2(27,27:end);
    
    
    theta_latent_b10_store(i,:)=theta_latent(:,1)';
    theta_latent_b20_store(i,:)=theta_latent(:,2)';
    theta_latent_b30_store(i,:)=theta_latent(:,3)';
    theta_latent_A0_store(i,:)=theta_latent(:,4)';
    theta_latent_v10_store(i,:)=theta_latent(:,5)';
    theta_latent_v210_store(i,:)=theta_latent(:,6)';
    theta_latent_v220_store(i,:)=theta_latent(:,7)';
    theta_latent_v230_store(i,:)=theta_latent(:,8)';
    theta_latent_tau0_store(i,:)=theta_latent(:,9)';
    
    theta_latent_b11_store(i,:)=theta_latent(:,10)';
    theta_latent_b21_store(i,:)=theta_latent(:,11)';
    theta_latent_b31_store(i,:)=theta_latent(:,12)';
    theta_latent_A1_store(i,:)=theta_latent(:,13)';
    theta_latent_v11_store(i,:)=theta_latent(:,14)';
    theta_latent_v211_store(i,:)=theta_latent(:,15)';
    theta_latent_v221_store(i,:)=theta_latent(:,16)';
    theta_latent_v231_store(i,:)=theta_latent(:,17)';
    theta_latent_tau1_store(i,:)=theta_latent(:,18)';
    
    theta_latent_b12_store(i,:)=theta_latent(:,19)';
    theta_latent_b22_store(i,:)=theta_latent(:,20)';
    theta_latent_b32_store(i,:)=theta_latent(:,21)';
    theta_latent_A2_store(i,:)=theta_latent(:,22)';
    theta_latent_v12_store(i,:)=theta_latent(:,23)';
    theta_latent_v212_store(i,:)=theta_latent(:,24)';
    theta_latent_v222_store(i,:)=theta_latent(:,25)';
    theta_latent_v232_store(i,:)=theta_latent(:,26)';
    theta_latent_tau2_store(i,:)=theta_latent(:,27)';

    
    a_half_store(i,:)=a_half';
    
    %after burn in, we count the number of unique values in the random
    %effects for each subject. 
    if i>burn
      for j=1:num_subjects
          count(1,j)=size(unique(theta_latent_A0_store(burn:end,j)),1);
      end  
    end
    
    %save the output to your directory
     %if i==500 | i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==8000 | i==9000 | i==10000 | i==20000 | i==30000 | i==40000 | i==50000 | i==60000 | i==70000 | i==80000 | i==90000 | i==100000
     if mod(i,1000)==0
        save('/short/jz21/dg2271/LBA_matchsearchstop_covariance_popOut_withprop.mat','theta_mu_store','theta_sig2_store1','theta_sig2_store2','theta_sig2_store3','theta_sig2_store4',...
             'theta_sig2_store5','theta_sig2_store6','theta_sig2_store7','theta_sig2_store8',...
             'theta_sig2_store9','theta_sig2_store10','theta_sig2_store11','theta_sig2_store12',...
             'theta_sig2_store13','theta_sig2_store14','theta_sig2_store15','theta_sig2_store16','theta_sig2_store17','theta_sig2_store18',...
             'theta_sig2_store19','theta_sig2_store20','theta_sig2_store21','theta_sig2_store22','theta_sig2_store23','theta_sig2_store24',...
             'theta_sig2_store25','theta_sig2_store26','theta_sig2_store27',...
             'chol_theta_sig2_store1','chol_theta_sig2_store2','chol_theta_sig2_store3',...
             'chol_theta_sig2_store4','chol_theta_sig2_store5','chol_theta_sig2_store6',...
             'chol_theta_sig2_store7','chol_theta_sig2_store8',...
             'chol_theta_sig2_store9','chol_theta_sig2_store10','chol_theta_sig2_store11',...
             'chol_theta_sig2_store12','chol_theta_sig2_store13','chol_theta_sig2_store14',...
             'chol_theta_sig2_store15','chol_theta_sig2_store16','chol_theta_sig2_store17','chol_theta_sig2_store18',...
             'chol_theta_sig2_store19','chol_theta_sig2_store20','chol_theta_sig2_store21','chol_theta_sig2_store22',...
             'chol_theta_sig2_store23','chol_theta_sig2_store24','chol_theta_sig2_store25','chol_theta_sig2_store26',...
             'chol_theta_sig2_store27',...
             'theta_latent_b10_store','theta_latent_b20_store','theta_latent_b30_store','theta_latent_A0_store','theta_latent_v10_store',...
             'theta_latent_v210_store','theta_latent_v220_store','theta_latent_v230_store','theta_latent_tau0_store',...
             'theta_latent_b11_store','theta_latent_b21_store','theta_latent_b31_store','theta_latent_A1_store','theta_latent_v11_store',...
             'theta_latent_v211_store','theta_latent_v221_store','theta_latent_v231_store',...
             'theta_latent_tau1_store',...
             'theta_latent_b12_store','theta_latent_b22_store','theta_latent_b32_store','theta_latent_A2_store','theta_latent_v12_store',...
             'theta_latent_v212_store','theta_latent_v222_store','theta_latent_v232_store',...
             'theta_latent_tau2_store','a_half_store'); 
     end
     i=i+1; 
     toc
end
% save the last 10000 draws for further analysis
length_draws=length(theta_latent_A0_store);
theta_latent_b10_store=theta_latent_b10_store(length_draws-9999:end,:);
theta_latent_b20_store=theta_latent_b20_store(length_draws-9999:end,:);
theta_latent_b30_store=theta_latent_b30_store(length_draws-9999:end,:);
theta_latent_A0_store=theta_latent_A0_store(length_draws-9999:end,:);
theta_latent_v10_store=theta_latent_v10_store(length_draws-9999:end,:);
theta_latent_v210_store=theta_latent_v210_store(length_draws-9999:end,:);
theta_latent_v220_store=theta_latent_v220_store(length_draws-9999:end,:);
theta_latent_v230_store=theta_latent_v230_store(length_draws-9999:end,:);
theta_latent_tau0_store=theta_latent_tau0_store(length_draws-9999:end,:);

theta_latent_b11_store=theta_latent_b11_store(length_draws-9999:end,:);
theta_latent_b21_store=theta_latent_b21_store(length_draws-9999:end,:);
theta_latent_b31_store=theta_latent_b31_store(length_draws-9999:end,:);
theta_latent_A1_store=theta_latent_A1_store(length_draws-9999:end,:);
theta_latent_v11_store=theta_latent_v11_store(length_draws-9999:end,:);
theta_latent_v211_store=theta_latent_v211_store(length_draws-9999:end,:);
theta_latent_v221_store=theta_latent_v221_store(length_draws-9999:end,:);
theta_latent_v231_store=theta_latent_v231_store(length_draws-9999:end,:);
theta_latent_tau1_store=theta_latent_tau1_store(length_draws-9999:end,:);

theta_latent_b12_store=theta_latent_b12_store(length_draws-9999:end,:);
theta_latent_b22_store=theta_latent_b22_store(length_draws-9999:end,:);
theta_latent_b32_store=theta_latent_b32_store(length_draws-9999:end,:);
theta_latent_A2_store=theta_latent_A2_store(length_draws-9999:end,:);
theta_latent_v12_store=theta_latent_v12_store(length_draws-9999:end,:);
theta_latent_v212_store=theta_latent_v212_store(length_draws-9999:end,:);
theta_latent_v222_store=theta_latent_v222_store(length_draws-9999:end,:);
theta_latent_v232_store=theta_latent_v232_store(length_draws-9999:end,:);
theta_latent_tau2_store=theta_latent_tau2_store(length_draws-9999:end,:);

theta_mu_store=theta_mu_store(length_draws-9999:end,:);
theta_sig2_store1=theta_sig2_store1(length_draws-9999:end,:);
theta_sig2_store2=theta_sig2_store2(length_draws-9999:end,:);
theta_sig2_store3=theta_sig2_store3(length_draws-9999:end,:);
theta_sig2_store4=theta_sig2_store4(length_draws-9999:end,:);
theta_sig2_store5=theta_sig2_store5(length_draws-9999:end,:);
theta_sig2_store6=theta_sig2_store6(length_draws-9999:end,:);
theta_sig2_store7=theta_sig2_store7(length_draws-9999:end,:);
theta_sig2_store8=theta_sig2_store8(length_draws-9999:end,:);
theta_sig2_store9=theta_sig2_store9(length_draws-9999:end,:);
theta_sig2_store10=theta_sig2_store10(length_draws-9999:end,:);
theta_sig2_store11=theta_sig2_store11(length_draws-9999:end,:);
theta_sig2_store12=theta_sig2_store12(length_draws-9999:end,:);
theta_sig2_store13=theta_sig2_store13(length_draws-9999:end,:);
theta_sig2_store14=theta_sig2_store14(length_draws-9999:end,:);
theta_sig2_store15=theta_sig2_store15(length_draws-9999:end,:);
theta_sig2_store16=theta_sig2_store16(length_draws-9999:end,:);
theta_sig2_store17=theta_sig2_store17(length_draws-9999:end,:);
theta_sig2_store18=theta_sig2_store18(length_draws-9999:end,:);
theta_sig2_store19=theta_sig2_store19(length_draws-9999:end,:);
theta_sig2_store20=theta_sig2_store20(length_draws-9999:end,:);
theta_sig2_store21=theta_sig2_store21(length_draws-9999:end,:);
theta_sig2_store22=theta_sig2_store22(length_draws-9999:end,:);
theta_sig2_store23=theta_sig2_store23(length_draws-9999:end,:);
theta_sig2_store24=theta_sig2_store24(length_draws-9999:end,:);
theta_sig2_store25=theta_sig2_store25(length_draws-9999:end,:);
theta_sig2_store26=theta_sig2_store26(length_draws-9999:end,:);
theta_sig2_store27=theta_sig2_store27(length_draws-9999:end,:);

chol_theta_sig2_store1=chol_theta_sig2_store1(length_draws-9999:end,:);
chol_theta_sig2_store2=chol_theta_sig2_store2(length_draws-9999:end,:);
chol_theta_sig2_store3=chol_theta_sig2_store3(length_draws-9999:end,:);
chol_theta_sig2_store4=chol_theta_sig2_store4(length_draws-9999:end,:);
chol_theta_sig2_store5=chol_theta_sig2_store5(length_draws-9999:end,:);
chol_theta_sig2_store6=chol_theta_sig2_store6(length_draws-9999:end,:);
chol_theta_sig2_store7=chol_theta_sig2_store7(length_draws-9999:end,:);
chol_theta_sig2_store8=chol_theta_sig2_store8(length_draws-9999:end,:);
chol_theta_sig2_store9=chol_theta_sig2_store9(length_draws-9999:end,:);
chol_theta_sig2_store10=chol_theta_sig2_store10(length_draws-9999:end,:);
chol_theta_sig2_store11=chol_theta_sig2_store11(length_draws-9999:end,:);
chol_theta_sig2_store12=chol_theta_sig2_store12(length_draws-9999:end,:);
chol_theta_sig2_store13=chol_theta_sig2_store13(length_draws-9999:end,:);
chol_theta_sig2_store14=chol_theta_sig2_store14(length_draws-9999:end,:);
chol_theta_sig2_store15=chol_theta_sig2_store15(length_draws-9999:end,:);
chol_theta_sig2_store16=chol_theta_sig2_store16(length_draws-9999:end,:);
chol_theta_sig2_store17=chol_theta_sig2_store17(length_draws-9999:end,:);
chol_theta_sig2_store18=chol_theta_sig2_store18(length_draws-9999:end,:);
chol_theta_sig2_store19=chol_theta_sig2_store19(length_draws-9999:end,:);
chol_theta_sig2_store20=chol_theta_sig2_store20(length_draws-9999:end,:);
chol_theta_sig2_store21=chol_theta_sig2_store21(length_draws-9999:end,:);
chol_theta_sig2_store22=chol_theta_sig2_store22(length_draws-9999:end,:);
chol_theta_sig2_store23=chol_theta_sig2_store23(length_draws-9999:end,:);
chol_theta_sig2_store24=chol_theta_sig2_store24(length_draws-9999:end,:);
chol_theta_sig2_store25=chol_theta_sig2_store25(length_draws-9999:end,:);
chol_theta_sig2_store26=chol_theta_sig2_store26(length_draws-9999:end,:);
chol_theta_sig2_store27=chol_theta_sig2_store27(length_draws-9999:end,:);

a_half_store=a_half_store(length_draws-9999:end,:);
%save the output to your directory
save('LBA_matchsearchstop_covariance_popOut_withprop.mat','theta_mu_store','theta_sig2_store1','theta_sig2_store2','theta_sig2_store3','theta_sig2_store4',...
             'theta_sig2_store5','theta_sig2_store6','theta_sig2_store7','theta_sig2_store8',...
             'theta_sig2_store9','theta_sig2_store10','theta_sig2_store11','theta_sig2_store12',...
             'theta_sig2_store13','theta_sig2_store14','theta_sig2_store15','theta_sig2_store16','theta_sig2_store17','theta_sig2_store18',...
             'theta_sig2_store19','theta_sig2_store20','theta_sig2_store21','theta_sig2_store22','theta_sig2_store23','theta_sig2_store24',...
             'theta_sig2_store25','theta_sig2_store26','theta_sig2_store27',...
             'chol_theta_sig2_store1','chol_theta_sig2_store2','chol_theta_sig2_store3',...
             'chol_theta_sig2_store4','chol_theta_sig2_store5','chol_theta_sig2_store6',...
             'chol_theta_sig2_store7','chol_theta_sig2_store8',...
             'chol_theta_sig2_store9','chol_theta_sig2_store10','chol_theta_sig2_store11',...
             'chol_theta_sig2_store12','chol_theta_sig2_store13','chol_theta_sig2_store14',...
             'chol_theta_sig2_store15','chol_theta_sig2_store16','chol_theta_sig2_store17','chol_theta_sig2_store18',...
             'chol_theta_sig2_store19','chol_theta_sig2_store20','chol_theta_sig2_store21','chol_theta_sig2_store22',...
             'chol_theta_sig2_store23','chol_theta_sig2_store24','chol_theta_sig2_store25','chol_theta_sig2_store26',...
             'chol_theta_sig2_store27',...
             'theta_latent_b10_store','theta_latent_b20_store','theta_latent_b30_store','theta_latent_A0_store','theta_latent_v10_store',...
             'theta_latent_v210_store','theta_latent_v220_store','theta_latent_v230_store','theta_latent_tau0_store',...
             'theta_latent_b11_store','theta_latent_b21_store','theta_latent_b31_store','theta_latent_A1_store','theta_latent_v11_store',...
             'theta_latent_v211_store','theta_latent_v221_store','theta_latent_v231_store',...
             'theta_latent_tau1_store',...
             'theta_latent_b12_store','theta_latent_b22_store','theta_latent_b32_store','theta_latent_A2_store','theta_latent_v12_store',...
             'theta_latent_v212_store','theta_latent_v222_store','theta_latent_v232_store',...
             'theta_latent_tau2_store','a_half_store'); 
         
% save('/short/jz21/dg2271/LBA_matchsearchstop_covariance_popOut_withprop.mat','theta_mu_store','theta_sig2_store1','theta_sig2_store2','theta_sig2_store3','theta_sig2_store4',...
%              'theta_sig2_store5','theta_sig2_store6','theta_sig2_store7','theta_sig2_store8',...
%              'theta_sig2_store9','theta_sig2_store10','theta_sig2_store11','theta_sig2_store12',...
%              'theta_sig2_store13','theta_sig2_store14','theta_sig2_store15','theta_sig2_store16','theta_sig2_store17','theta_sig2_store18',...
%              'theta_sig2_store19','theta_sig2_store20','theta_sig2_store21','theta_sig2_store22','theta_sig2_store23','theta_sig2_store24',...
%              'theta_sig2_store25','theta_sig2_store26','theta_sig2_store27',...
%              'chol_theta_sig2_store1','chol_theta_sig2_store2','chol_theta_sig2_store3',...
%              'chol_theta_sig2_store4','chol_theta_sig2_store5','chol_theta_sig2_store6',...
%              'chol_theta_sig2_store7','chol_theta_sig2_store8',...
%              'chol_theta_sig2_store9','chol_theta_sig2_store10','chol_theta_sig2_store11',...
%              'chol_theta_sig2_store12','chol_theta_sig2_store13','chol_theta_sig2_store14',...
%              'chol_theta_sig2_store15','chol_theta_sig2_store16','chol_theta_sig2_store17','chol_theta_sig2_store18',...
%              'chol_theta_sig2_store19','chol_theta_sig2_store20','chol_theta_sig2_store21','chol_theta_sig2_store22',...
%              'chol_theta_sig2_store23','chol_theta_sig2_store24','chol_theta_sig2_store25','chol_theta_sig2_store26',...
%              'chol_theta_sig2_store27',...
%              'theta_latent_b10_store','theta_latent_b20_store','theta_latent_b30_store','theta_latent_A0_store','theta_latent_v10_store',...
%              'theta_latent_v210_store','theta_latent_v220_store','theta_latent_v230_store','theta_latent_tau0_store',...
%              'theta_latent_b11_store','theta_latent_b21_store','theta_latent_b31_store','theta_latent_A1_store','theta_latent_v11_store',...
%              'theta_latent_v211_store','theta_latent_v221_store','theta_latent_v231_store',...
%              'theta_latent_tau1_store',...
%              'theta_latent_b12_store','theta_latent_b22_store','theta_latent_b32_store','theta_latent_A2_store','theta_latent_v12_store',...
%              'theta_latent_v212_store','theta_latent_v222_store','theta_latent_v232_store',...
%              'theta_latent_tau2_store','a_half_store'); 