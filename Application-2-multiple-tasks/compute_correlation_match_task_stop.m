%compute correlation matrix

load('LBA_matchsearchstop_covariance_popOut_withprop.mat');
length_draws=length(theta_mu_store(:,1));

for i=1:length_draws
    i
    %mean lognormal
    mu_b10(i,1)=exp(theta_mu_store(i,1)+0.5*theta_sig2_store1(i,1));
    mu_b20(i,1)=exp(theta_mu_store(i,2)+0.5*theta_sig2_store2(i,1));
    mu_b30(i,1)=exp(theta_mu_store(i,3)+0.5*theta_sig2_store3(i,1));
    mu_A0(i,1)=exp(theta_mu_store(i,4)+0.5*theta_sig2_store4(i,1));
    mu_v10(i,1)=exp(theta_mu_store(i,5)+0.5*theta_sig2_store5(i,1));
    mu_v210(i,1)=exp(theta_mu_store(i,6)+0.5*theta_sig2_store6(i,1));
    mu_v220(i,1)=exp(theta_mu_store(i,7)+0.5*theta_sig2_store7(i,1));
    mu_v230(i,1)=exp(theta_mu_store(i,8)+0.5*theta_sig2_store8(i,1));
    mu_tau0(i,1)=exp(theta_mu_store(i,9)+0.5*theta_sig2_store9(i,1));
    
    mu_b11(i,1)=exp(theta_mu_store(i,10)+0.5*theta_sig2_store10(i,1));
    mu_b21(i,1)=exp(theta_mu_store(i,11)+0.5*theta_sig2_store11(i,1));
    mu_b31(i,1)=exp(theta_mu_store(i,12)+0.5*theta_sig2_store12(i,1));
    mu_A1(i,1)=exp(theta_mu_store(i,13)+0.5*theta_sig2_store13(i,1));
    mu_v11(i,1)=exp(theta_mu_store(i,14)+0.5*theta_sig2_store14(i,1));
    mu_v211(i,1)=exp(theta_mu_store(i,15)+0.5*theta_sig2_store15(i,1));
    mu_v221(i,1)=exp(theta_mu_store(i,16)+0.5*theta_sig2_store16(i,1));
    mu_v231(i,1)=exp(theta_mu_store(i,17)+0.5*theta_sig2_store17(i,1));    
    mu_tau1(i,1)=exp(theta_mu_store(i,18)+0.5*theta_sig2_store18(i,1));
        
    mu_b12(i,1)=exp(theta_mu_store(i,19)+0.5*theta_sig2_store19(i,1));
    mu_b22(i,1)=exp(theta_mu_store(i,20)+0.5*theta_sig2_store20(i,1));
    mu_b32(i,1)=exp(theta_mu_store(i,21)+0.5*theta_sig2_store21(i,1));
    mu_A2(i,1)=exp(theta_mu_store(i,22)+0.5*theta_sig2_store22(i,1));
    mu_v12(i,1)=exp(theta_mu_store(i,23)+0.5*theta_sig2_store23(i,1));
    mu_v212(i,1)=exp(theta_mu_store(i,24)+0.5*theta_sig2_store24(i,1));
    mu_v222(i,1)=exp(theta_mu_store(i,25)+0.5*theta_sig2_store25(i,1));
    mu_v232(i,1)=exp(theta_mu_store(i,26)+0.5*theta_sig2_store26(i,1));    
    mu_tau2(i,1)=exp(theta_mu_store(i,27)+0.5*theta_sig2_store27(i,1));
    
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
    chol_sig2(15,1:15)=[chol_theta_sig2_store15(i,1:14),exp(chol_theta_sig2_store15(i,15))];
    chol_sig2(16,1:16)=[chol_theta_sig2_store16(i,1:15),exp(chol_theta_sig2_store16(i,16))];
    chol_sig2(17,1:17)=[chol_theta_sig2_store17(i,1:16),exp(chol_theta_sig2_store17(i,17))];
    chol_sig2(18,1:18)=[chol_theta_sig2_store18(i,1:17),exp(chol_theta_sig2_store18(i,18))];
    chol_sig2(19,1:19)=[chol_theta_sig2_store19(i,1:18),exp(chol_theta_sig2_store19(i,19))];
    chol_sig2(20,1:20)=[chol_theta_sig2_store20(i,1:19),exp(chol_theta_sig2_store20(i,20))];
    chol_sig2(21,1:21)=[chol_theta_sig2_store21(i,1:20),exp(chol_theta_sig2_store21(i,21))];
    chol_sig2(22,1:22)=[chol_theta_sig2_store22(i,1:21),exp(chol_theta_sig2_store22(i,22))];
    chol_sig2(23,1:23)=[chol_theta_sig2_store23(i,1:22),exp(chol_theta_sig2_store23(i,23))];
    chol_sig2(24,1:24)=[chol_theta_sig2_store24(i,1:23),exp(chol_theta_sig2_store24(i,24))];
    chol_sig2(25,1:25)=[chol_theta_sig2_store25(i,1:24),exp(chol_theta_sig2_store25(i,25))];
    chol_sig2(26,1:26)=[chol_theta_sig2_store26(i,1:25),exp(chol_theta_sig2_store26(i,26))];
    chol_sig2(27,1:27)=[chol_theta_sig2_store27(i,1:26),exp(chol_theta_sig2_store27(i,27))];   
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
    covmat_LN(15,:,i)=exp(theta_mu_store(i,15)+theta_mu_store(i,:)+0.5.*(covmat_normal(15,15,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(15,:,i))-1);
    covmat_LN(16,:,i)=exp(theta_mu_store(i,16)+theta_mu_store(i,:)+0.5.*(covmat_normal(16,16,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(16,:,i))-1);
    covmat_LN(17,:,i)=exp(theta_mu_store(i,17)+theta_mu_store(i,:)+0.5.*(covmat_normal(17,17,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(17,:,i))-1);
    covmat_LN(18,:,i)=exp(theta_mu_store(i,18)+theta_mu_store(i,:)+0.5.*(covmat_normal(18,18,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(18,:,i))-1);
    covmat_LN(19,:,i)=exp(theta_mu_store(i,19)+theta_mu_store(i,:)+0.5.*(covmat_normal(19,19,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(19,:,i))-1);
    covmat_LN(20,:,i)=exp(theta_mu_store(i,20)+theta_mu_store(i,:)+0.5.*(covmat_normal(20,20,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(20,:,i))-1);
    covmat_LN(21,:,i)=exp(theta_mu_store(i,21)+theta_mu_store(i,:)+0.5.*(covmat_normal(21,21,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(21,:,i))-1);
    covmat_LN(22,:,i)=exp(theta_mu_store(i,22)+theta_mu_store(i,:)+0.5.*(covmat_normal(22,22,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(22,:,i))-1);
    covmat_LN(23,:,i)=exp(theta_mu_store(i,23)+theta_mu_store(i,:)+0.5.*(covmat_normal(23,23,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(23,:,i))-1);
    covmat_LN(24,:,i)=exp(theta_mu_store(i,24)+theta_mu_store(i,:)+0.5.*(covmat_normal(24,24,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(24,:,i))-1);
    covmat_LN(25,:,i)=exp(theta_mu_store(i,25)+theta_mu_store(i,:)+0.5.*(covmat_normal(25,25,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(25,:,i))-1);
    covmat_LN(26,:,i)=exp(theta_mu_store(i,26)+theta_mu_store(i,:)+0.5.*(covmat_normal(26,26,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(26,:,i))-1);
    covmat_LN(27,:,i)=exp(theta_mu_store(i,27)+theta_mu_store(i,:)+0.5.*(covmat_normal(27,27,i)+diag(covmat_normal(:,:,i))')).*(exp(covmat_normal(27,:,i))-1);

    
    [~,corr_LN(:,:,i)]=cov2corr(covmat_LN(:,:,i));
end
mu_LN=[mu_b10,mu_b20,mu_b30,mu_A0,mu_v10,mu_v210,mu_v220,mu_v230,mu_tau0,mu_b11,mu_b21,mu_b31,mu_A1,mu_v11,mu_v211,mu_v221,mu_v231,mu_tau1,...
       mu_b12,mu_b22,mu_b32,mu_A2,mu_v12,mu_v212,mu_v222,mu_v232,mu_tau2];
save('result_matchsearchstop_fullcovariance_popOut.mat','mu_LN','covmat_LN','corr_LN');


