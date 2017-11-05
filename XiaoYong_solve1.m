function [ h_opt ] = XiaoYong_solve1( x,y,A,Ex,Ey,Vx,Vy,COV )
%给现货期货各自的数据x y,以及投资者情绪A，求解效用函数得到h_opt
%  mx（i）是x的i阶矩，my（i）是y的i阶矩
   
    Mx(1)=Ex;
    My(1)=Ey;
    Mx(2)=Vx;
    My(2)=Vy;
        for i=3:4
            Mx(i)=mean((x-Mx(1)).^i);
            My(i)=mean((y-My(1)).^i);
        end
    
   
function [ CoM ] = LCoM(nx,ny )
%假定线性关系求协矩(不超过4阶)x取nx阶 y取ny阶
    function[ alpha,beta ]=HG(x,y)
      n=length(x);
      beta=(dot(x,y)-n*Mx(1)*My(1))/(dot(x,x)-n*Mx(1)*My(1));
      alpha=My(1)-beta*Mx(1);
    end
[alpha,beta]=HG(x,y);%最小二乘法估计的系数
k=nx+ny;
    if nx==1
        CoM=beta*My(k);
    elseif nx==2
        CoM=2*alpha*beta*My(k-1)+beta^2*My(k);
    elseif ny==1
        CoM=beta*Mx(k);
    elseif ny==2
        CoM=2*(-alpha/beta)*(1/beta)*Mx(k-1)+(1/beta)^2*Mx(k);
    end
end

%要用的均值
Ef=My(1);
%求方差
Drs= Mx(2);
Drf= My(2);
sigma_sf= COV;
%求三阶矩
Srs=Mx(3);
Srf=My(3);
Ssff=LCoM( 1,2 );
Sssf=LCoM( 2,1 );
%求四阶矩
Krs=Mx(4);
Krf=My(4);
Ksfff=LCoM( 1,3 );
Kssff=LCoM( 2,2 );
Ksssf=LCoM( 3,1 );

        Y4=A^4*Krf/24;
        X3=A^4*Krf/6;
        Y3=A^3*Srf/6-A^4*Ksfff/6;
        X2=(A^3*Srf-A^4*Ksfff)/2;
        Y2=A^2*Drf/2-A^3*Ssff/2+Kssff/4;
        X1=A^2*Drf-A^3*Ssff+Kssff/2;
        Y1=A^3*Sssf/2-A^2*sigma_sf-A^4*Ksssf/6;
        X0=A^3*Sssf/2-A^2*sigma_sf-A^4*Ksssf/6;
        Y0=1+A^2*Drs/2-A^3*Srs/6+A^4*Krs/24;
    function [ y ]=XiaoYong_partial( h )
        y=X3.*h.^3+X2.*h.^2+X1.*h+X0+A.*Ef.*(Y4.*h.^4+Y3.*h.^3+Y2.*h.^2+Y1.*h+Y0);
    end
h_opt=fsolve(@(h) XiaoYong_partial(h),1);

end

