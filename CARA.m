function [cara] =CARA(var,delta,rf )
%��������ϵ��
%input
%varΪt+ʱ�̼�֮ǰ������3�������գ�rf����������,������
%deltaΪ�����������յ�t+1ʱ�̵�Ͷ��������ָ����������
%rfΪ�����������յ�t+1ʱ�������ʣ�������
%output
%tʱ�̵ķ������ϵ��
X=[delta var var.*delta];
B=regress(rf,X);
cara=B(3)+B(4)*delta(end);
end

