
%Just copy and paste this into your Matlab window for greater ease. The GARCH_code.m found on the homepage will look better thanks to proper spacing. This is not meant to be run as command line.

%Garth Mortensen

%%
%%DESCRIPTION
%Bivariate GARCH model

%REQUIREMENTS
%This code requires the James P. LeSage Econometrics Toolbox and UCSD
%GARCH toolboxes to run properly. Verify you have them installed using command 'ver'
%Install/uninstall toolboxes using command 'pathtool.' This code doesnt use
%the adftest that comes with the Econometrics toolbox. Instead, the original 
% Matlab adftest is chosen due to its ease of use.


%% WIPE
%wipe the memory

clear all
close all
clc


%%
%This code is used to read data from various excel forms.

%///CHANGE EXCEL DATE FORMAT TO GENERAL, NOT STRING///



%This is for importing from Datastream. note this reads .xls, not .xlsx

rs0     = xlsread('rate.xlsx','A1:A732');
rf0     = xlsread('rate.xlsx','B1:B732');
rp=zeros(732,1);
%���������
corrboth = corr(rs0, rf0);
i=1;
%ע��������tʱ�̶��ǵ�t+1��������
for t=633:732
    rs=rs0(1:t);
    rf=rf0(1:t);
    % ARMA filter
    %��һ�������������ֵ��ʹ��������ǰ���������ֵΪ�㣩�������Ļ���Garchģ�͵�����������ܸ��õĹ���
    %Pull out the conditional mean with ARMA. 
    [ARMAparams_rs,ARMAerrors_rs,~,ARMAsterrors_rs,~,~,~,~]           = armaxfilter(rs,1,1,1);
    [ARMAparams_rf,ARMAerrors_rf,~,ARMAsterrors_rf,~,~,~,~]           = armaxfilter(rf,1,1,1);


    %ARMA parameters %����ARMAģ�͵Ĳ���
    %ARMAconstant_rs          = ARMAparams_rs(1,1);
    %ARMA_AR_rs                = ARMAparams_rs(1,2);
    %ARMA_MA_rs                = ARMAparams_rs(1,3);

    %ARMAconstant_rf          = ARMAparams_rf(1,1);
    %ARMA_AR_rf                = ARMAparams_rf(1,2);
    %ARMA_MA_rf                = ARMAparams_rf(1,3);

    % WRITE ARMA PARAMS
    %rs
    %xlswrite('record.xlsx',ARMAconstant_rs,'ARMA','c4');
    %xlswrite('record.xlsx',ARMA_AR_rs,'ARMA','C5');
    %xlswrite('record.xlsx',ARMA_MA_rs,'ARMA','C6');

    %rf
    %xlswrite('record.xlsx',ARMAconstant_rf,'ARMA','E4');
    %xlswrite('record.xlsx',ARMA_AR_rf,'ARMA','E5');
    %xlswrite('record.xlsx',ARMA_MA_rf,'ARMA','E6');


   % GARCH
   %Using the residuals from the ARMA model, estimate GARCH parameters.
   %����ARMAģ����ò��������� GARCHģ�Ͳ���
   %Pull out the conditional variance with GARCH.
   [GARCHpqparameters_rs,GARCHpqmaxliklihood_rs,GARCHpgvariances_rs,GARCHpgstderror_rs,GARCHpgscores_rs,~]    = garchpq(ARMAerrors_rs,1,1);
   [GARCHpqparameters_rf,GARCHpqmaxliklihood_rf,GARCHpgvariances_rf,GARCHpgstderror_rf,GARCHpgscores_rf,~]    = garchpq(ARMAerrors_rf,1,1);
   
   % MVGARCH COMBINE VECTORS
   %Prepare an MVGARCH matrix from the ARMA errors for CC-GARCH
   both = [ARMAerrors_rs,ARMAerrors_rf];
   %DCC-GARCH
   [DCCparameters,DCClogliklihood,Ht,Qt,~,~,DCCstderror,~,~,~]  = dcc_mvgarch(both,1,1,1,1);

   cara=CARA(GARCHpgvariances_rf(1:end),delta,rf(2:end));
   [Es,Ef]=mean(ARMAparams_rs,ARMAparams_rf,rs,rf,ARMAerrors_rs,ARMAerrors_rf,t);%������ֵ
   Vs=GARCHpgvariances_rs(t-1);
   Vf=GARCHpgvariances_rf(t-1);
   Vsf=Ht(1,2,t-1);
   h(i)=XiaoYong_solve1(rs,rf,cara,Es,Ef,Vs,Vf,Vsf );  %��Գ����
   rp(t)=rs0(t)-rf0(t)*h(i);  %�Գ�������ǰ632����0��
   i=i+1;
end

%֮���бȽ϶Գ�Ч���Ĵ��롭����ʯ��һ��











