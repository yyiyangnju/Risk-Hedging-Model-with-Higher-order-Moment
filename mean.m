function [Es,Ef] =mean(ARMAparams_rs,ARMAparams_rf,rs,rf,error_rs,error_rf,t)
%����ARMAģ�͵õ��Ķ�����������ʣ�ARMAģ�͵õ��Ĳв���е������tʱ�̵�Ԥ��ֵ
%�ڴ���֮ǰ����Щrf,rs,error_rs,error_rf���Ѿ���t+1ʱ�̽׶Σ�Ҳ����˵��
%t��rf,rs��end,t-1��error_rs,error_rf��end
%˵��tʱ�̣�rs,rf�ĳ���Ϊt,error_rs,error_rf�ĳ���Ϊt-1
Es=ARMAparams_rs(1,1)+error_rs(t-1)+ARMAparams_rs(1,2)*rs(t-1)+ARMAparams_rs(1,3)*error_rs(t-2);
Ef=ARMAparams_rf(1,1)+error_rf(t-1)+ARMAparams_rf(1,2)*rf(t-1)+ARMAparams_rf(1,3)*error_rf(t-2);

end

