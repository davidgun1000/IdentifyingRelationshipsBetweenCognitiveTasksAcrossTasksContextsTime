%compute correlation matrix

load('LBA_Forstmann_covariance.mat');
length_draws=length(theta_mu_store(:,1));

for i=1:length_draws
    i
    %mean lognormal
    mu_b10(i,1)=exp(theta_mu_store(i,1)+0.5*theta_sig2_store1(i,1));
    mu_b20(i,1)=exp(theta_mu_store(i,2)+0.5*theta_sig2_store2(i,1));
    mu_b30(i,1)=exp(theta_mu_store(i,3)+0.5*theta_sig2_store3(i,1));
    mu_A0(i,1)=exp(theta_mu_store(i,4)+0.5*theta_sig2_store4(i,1));
    mu_v10(i,1)=exp(theta_mu_store(i,5)+0.5*theta_sig2_store5(i,1));
    mu_v20(i,1)=exp(theta_mu_store(i,6)+0.5*theta_sig2_store6(i,1));
    mu_tau0(i,1)=exp(theta_mu_store(i,7)+0.5*theta_sig2_store7(i,1));
    
    mu_b11(i,1)=exp(theta_mu_store(i,8)+0.5*theta_sig2_store8(i,1));
    mu_b21(i,1)=exp(theta_mu_store(i,9)+0.5*theta_sig2_store9(i,1));
    mu_b31(i,1)=exp(theta_mu_store(i,10)+0.5*theta_sig2_store10(i,1));
    mu_A1(i,1)=exp(theta_mu_store(i,11)+0.5*theta_sig2_store11(i,1));
    mu_v11(i,1)=exp(theta_mu_store(i,12)+0.5*theta_sig2_store12(i,1));
    mu_v21(i,1)=exp(theta_mu_store(i,13)+0.5*theta_sig2_store13(i,1));
    mu_tau1(i,1)=exp(theta_mu_store(i,14)+0.5*theta_sig2_store14(i,1));
        

    %covmat normal
    chol_sig2(1,1)=exp(chol_theta_sig2_store1(i,1));
    chol_sig2(2,1:2)=[chol_theta_sig2_store2(i,1),exp(chol_theta_sig2_store2(i,2))];
    chol_sig2(3,1:3)=[chol_theta_sig2_store3(i,1:2),exp(chol_theta_sig2_store3(i,3))];
    chol_sig2(4,1:4)=[chol_theta_sig2_store4(i,1:3),exp(chol_theta_sig2_store4(i,4))];
    chol_sig2(5,1:5)=[chol_theta_sig2_store5(i,1:4),exp(chol_theta_sig2_store5(i,5))];
    chol_sig2(6,1:6)=[chol_theta_sig2_store6(i,1:5),exp(chol_theta_sig2_store6(i,6))];
    chol_sig2(7,1:7)=[chol_theta_sig2_store7(i,1:6),exp(chol_theta_sig2_store7(i,7))];
    chol_sig2(8,1:8)=[chol_theta_sig2_store8(i,1:7),exp(chol_theta_sig2_store8(i,8))];
    chol_sig2(9,1:9)=[chol_theta_sig2_store9(i,1:8),exp(chol_theta_sig2_store9(i,9))];
    chol_sig2(10,1:10)=[chol_theta_sig2_store10(i,1:9),exp(chol_theta_sig2_store10(i,10))];
    chol_sig2(11,1:11)=[chol_theta_sig2_store11(i,1:10),exp(chol_theta_sig2_store11(i,11))];
    chol_sig2(12,1:12)=[chol_theta_sig2_store12(i,1:11),exp(chol_theta_sig2_store12(i,12))];
    chol_sig2(13,1:13)=[chol_theta_sig2_store13(i,1:12),exp(chol_theta_sig2_store13(i,13))];
    chol_sig2(14,1:14)=[chol_theta_sig2_store14(i,1:13),exp(chol_theta_sig2_store14(i,14))];
    covmat_normal(:,:,i)=chol_sig2*chol_sig2';
    
    %covmat lognormal
    
    covmat_LN(1,:,i)=exp(theta_mu_store(i,1)+theta_mu_store(i,:)+0.5.*(covmat_normal(1,1,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(1,:,i))-1);
    covmat_LN(2,:,i)=exp(theta_mu_store(i,2)+theta_mu_store(i,:)+0.5.*(covmat_normal(2,2,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(2,:,i))-1);
    covmat_LN(3,:,i)=exp(theta_mu_store(i,3)+theta_mu_store(i,:)+0.5.*(covmat_normal(3,3,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(3,:,i))-1);
    covmat_LN(4,:,i)=exp(theta_mu_store(i,4)+theta_mu_store(i,:)+0.5.*(covmat_normal(4,4,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(4,:,i))-1);
    covmat_LN(5,:,i)=exp(theta_mu_store(i,5)+theta_mu_store(i,:)+0.5.*(covmat_normal(5,5,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(5,:,i))-1);
    covmat_LN(6,:,i)=exp(theta_mu_store(i,6)+theta_mu_store(i,:)+0.5.*(covmat_normal(6,6,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(6,:,i))-1);
    covmat_LN(7,:,i)=exp(theta_mu_store(i,7)+theta_mu_store(i,:)+0.5.*(covmat_normal(7,7,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(7,:,i))-1);
    covmat_LN(8,:,i)=exp(theta_mu_store(i,8)+theta_mu_store(i,:)+0.5.*(covmat_normal(8,8,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(8,:,i))-1);
    covmat_LN(9,:,i)=exp(theta_mu_store(i,9)+theta_mu_store(i,:)+0.5.*(covmat_normal(9,9,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(9,:,i))-1);
    covmat_LN(10,:,i)=exp(theta_mu_store(i,10)+theta_mu_store(i,:)+0.5.*(covmat_normal(10,10,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(10,:,i))-1);
    covmat_LN(11,:,i)=exp(theta_mu_store(i,11)+theta_mu_store(i,:)+0.5.*(covmat_normal(11,11,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(11,:,i))-1);
    covmat_LN(12,:,i)=exp(theta_mu_store(i,12)+theta_mu_store(i,:)+0.5.*(covmat_normal(12,12,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(12,:,i))-1);
    covmat_LN(13,:,i)=exp(theta_mu_store(i,13)+theta_mu_store(i,:)+0.5.*(covmat_normal(13,13,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(13,:,i))-1);
    covmat_LN(14,:,i)=exp(theta_mu_store(i,14)+theta_mu_store(i,:)+0.5.*(covmat_normal(14,14,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(14,:,i))-1);
    
    
    [~,corr_LN(:,:,i)]=cov2corr(covmat_LN(:,:,i));
end
mu_LN=[mu_b10,mu_b20,mu_b30,mu_A0,mu_v10,mu_v20,mu_tau0,mu_b11,mu_b21,mu_b31,mu_A1,mu_v11,mu_v21,mu_tau1];
save('result_Forstmann_fullcovariance.mat','mu_LN','covmat_LN','corr_LN');


