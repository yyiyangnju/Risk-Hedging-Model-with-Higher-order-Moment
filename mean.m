function [Es,Ef] =mean(ARMAparams_rs,ARMAparams_rf,rs,rf,error_rs,error_rf,t)
%利用ARMA模型得到的额参数，收益率，ARMA模型得到的残差进行迭代求得t时刻的预测值
%在带入之前，这些rf,rs,error_rs,error_rf就已经在t+1时刻阶段，也就是说，
%t是rf,rs的end,t-1是error_rs,error_rf的end
%说是t时刻，rs,rf的长度为t,error_rs,error_rf的长度为t-1
Es=ARMAparams_rs(1,1)+error_rs(t-1)+ARMAparams_rs(1,2)*rs(t-1)+ARMAparams_rs(1,3)*error_rs(t-2);
Ef=ARMAparams_rf(1,1)+error_rf(t-1)+ARMAparams_rf(1,2)*rf(t-1)+ARMAparams_rf(1,3)*error_rf(t-2);

end

