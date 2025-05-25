[fv_mu_b1_1,xv_mu_b1_1]=ksdensity(mu_LN(:,1));
[fv_mu_b2_1,xv_mu_b2_1]=ksdensity(mu_LN(:,2));
[fv_mu_b3_1,xv_mu_b3_1]=ksdensity(mu_LN(:,3));
[fv_mu_A_1,xv_mu_A_1]=ksdensity(mu_LN(:,4));
[fv_mu_v1_1,xv_mu_v1_1]=ksdensity(mu_LN(:,5));
[fv_mu_v2_1,xv_mu_v2_1]=ksdensity(mu_LN(:,6));
[fv_mu_tau_1,xv_mu_tau_1]=ksdensity(mu_LN(:,7));

[fv_mu_b1_2,xv_mu_b1_2]=ksdensity(mu_LN(:,8));
[fv_mu_b2_2,xv_mu_b2_2]=ksdensity(mu_LN(:,9));
[fv_mu_b3_2,xv_mu_b3_2]=ksdensity(mu_LN(:,10));
[fv_mu_A_2,xv_mu_A_2]=ksdensity(mu_LN(:,11));
[fv_mu_v1_2,xv_mu_v1_2]=ksdensity(mu_LN(:,12));
[fv_mu_v2_2,xv_mu_v2_2]=ksdensity(mu_LN(:,13));
[fv_mu_tau_2,xv_mu_tau_2]=ksdensity(mu_LN(:,14));

subplot(2,4,1);plot(xv_mu_b1_1,fv_mu_b1_1,'linewidth',3,'color','b');hold on;plot(xv_mu_b1_2,fv_mu_b1_2,'linewidth',3,'color','r');title('\mu_{\alpha_{b^{1}}}');legend('OutScanner','InScanner');
subplot(2,4,2);plot(xv_mu_b2_1,fv_mu_b2_1,'linewidth',3,'color','b');hold on;plot(xv_mu_b2_2,fv_mu_b2_2,'linewidth',3,'color','r');title('\mu_{\alpha_{b^{2}}}');legend('OutScanner','InScanner');
subplot(2,4,3);plot(xv_mu_b3_1,fv_mu_b3_1,'linewidth',3,'color','b');hold on;plot(xv_mu_b3_2,fv_mu_b3_2,'linewidth',3,'color','r');title('\mu_{\alpha_{b^{3}}}');legend('OutScanner','InScanner');
subplot(2,4,4);plot(xv_mu_A_1,fv_mu_A_1,'linewidth',3,'color','b');hold on;plot(xv_mu_A_2,fv_mu_A_2,'linewidth',3,'color','r');title('\mu_{\alpha_{A}}');legend('OutScanner','InScanner');
subplot(2,4,5);plot(xv_mu_v1_1,fv_mu_v1_1,'linewidth',3,'color','b');hold on;plot(xv_mu_v1_2,fv_mu_v1_2,'linewidth',3,'color','r');title('\mu_{\alpha_{v^{1}}}');legend('OutScanner','InScanner');
subplot(2,4,6);plot(xv_mu_v2_1,fv_mu_v2_1,'linewidth',3,'color','b');hold on;plot(xv_mu_v2_2,fv_mu_v2_2,'linewidth',3,'color','r');title('\mu_{\alpha_{v^{2}}}');legend('OutScanner','InScanner');
subplot(2,4,7);plot(xv_mu_tau_1,fv_mu_tau_1,'linewidth',3,'color','b');hold on;plot(xv_mu_tau_2,fv_mu_tau_2,'linewidth',3,'color','r');title('\mu_{\alpha_{\tau}}');legend('OutScanner','InScanner');
