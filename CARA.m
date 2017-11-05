function [cara] =CARA(var,delta,rf )
%求风险厌恶系数
%input
%var为t+时刻及之前（到第3个交易日）rf的条件方差,列向量
%delta为第三个交易日到t+1时刻的投资者情绪指数，列向量
%rf为第三个交易日到t+1时刻收益率，列向量
%output
%t时刻的风险厌恶系数
X=[delta var var.*delta];
B=regress(rf,X);
cara=B(3)+B(4)*delta(end);
end

