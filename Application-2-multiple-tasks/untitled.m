for i=1:27
    mean_temp(i,1)=mean(covmat_LN(i,i,:));
    std_temp(i,1)=std(covmat_LN(i,i,:));
end
   