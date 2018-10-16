clear
clc

%%%%基于最佳中继的多用户调度%%%%
%%%%十个主用户,八个次用户,八个次要中继%%%%

%%%%各类参数（pu为主用户,su为次用户,sr为次要中继,r为主接收端,d为次接收端）%%%%
%%%%I=20dBm%%%%
%主用户参数
M=10;                  %主用户数（大于或等于2）
I=0.1;                 %主用户的干扰可容值  
Pp=10;                 %主用户功率  
Sigma_hir=0.6;         %pu->r,E(|hir|^2)=Sigma_ir
Sigma_hid=0.2;         %pu->d,E(|hid|^2)=Sigma_id
Sigma_hib=0.2;         %pu->sr,E(|hib|^2)=Sigma_ib
Theta_ir=1;
Theta_ie=0.7;

%次用户参数
N=8;                   %次用户数（大于或等于2）
P=10^0.2;              %次用户的分配功率P
Sigma_hjb=1;           %su->sr,E(|hjb|^2)=Sigma_jb
Sigma_hjr=0.2;         %su->r,E(|hjr|^2)=Sigma_jr
Theta_jb=1;
Theta_je=0.7;

%次要中继参数
R=8;                   %次要中继数（大于或等于2）
P=10^0.2;              %次要中继的分配功率P
Sigma_hbd=1;           %sr->d,E(|hbd|^2)=Sigma_bd
Sigma_hbr=0.2;         %sr->r,E(|hbr|^2)=Sigma_br
Theta_bd=1;
Theta_be=0.7;

%其他参数
L=10000;               %实验次数
N0=0.1;                %加性高斯白噪声功率
k=1;

for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio             
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir*Theta_ie/(Theta_ir*lambda);
    Sigma_hje=Sigma_hjb*Theta_je/(Theta_jb*lambda);
    Sigma_hbe=Sigma_hbd*Theta_be/(Theta_bd*lambda);

    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次要接收端的信道
    hib=sqrt(Sigma_hib/2)*randn(M,L)+sqrt(-Sigma_hib/2)*randn(M,L);        %主用户到次要中继的信道
      
    %次用户信道衰落系数%%%%
    hjb=sqrt(Sigma_hjb/2)*randn(N,L)+sqrt(-Sigma_hjb/2)*randn(N,L);        %次用户到最佳次要中继的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道
  
    %次要中继信道衰落系数%%%%
    hbd=sqrt(Sigma_hbd/2)*randn(R,L)+sqrt(-Sigma_hbd/2)*randn(R,L);        %次要中继到次要接收端的信道
    hbe=sqrt(Sigma_hbe/2)*randn(R,L)+sqrt(-Sigma_hbe/2)*randn(R,L);        %次要中继到窃听者的信道
    hbr=sqrt(Sigma_hbr/2)*randn(R,L)+sqrt(-Sigma_hbr/2)*randn(R,L);        %次要中继到主接收端的信道
    
    %%%%第二个时隙时主用户及次要中继的保密容量%%%%
    Temp2=zeros(M*R,L);
    Temp3=zeros(M*R,L);
    Temp4=zeros(R*M,L);
    Temp5=zeros(R*M,L);
    Temp_pu_sc_2=zeros(M*R,L);
    Temp_sr_sc=zeros(M*R,L);
    M_pu_R_sr_joint_sc=zeros(M*R,L);
    Pb=zeros(R,L);
   
    %%%%次要中继发射功率的限制%%%%
    for r=1:R
        for m=1:L 
            Pb(r,m)=min(P,I/(abs(hbr(r,m)^2)));                            %主要受P限制
        end
    end
        
    %%%%进行M*R次搜索对最佳次要中继进行选取%%%%
    a=1;
    for i=1:1:M
        for r=1:1:R
           Temp2(r,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Pb(r,:).*(abs(hbr(r,:)).^2)+N0));                                                   %主用户至主接收端的主信道容量
           Temp3(r,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Pb(r,:).*(abs(hbe(r,:)).^2)+N0));                                                   %主用户的窃听信道容量
           Temp4(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbd(r,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));                                               %次要中继至次接收端的主信道容量
           Temp5(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbe(r,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));                                               %次要中继的窃听信道容量
           Temp_pu_sc_2(a,:)=(Temp2(r,:)-Temp3(r,:)<0).*0+(Temp2(r,:)-Temp3(r,:)>=0).*(Temp2(r,:)-Temp3(r,:));                          %主用户的保密容量：若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_sr_sc(a,:)=(Temp4(r,:)-Temp5(r,:)<0).*0+(Temp4(r,:)-Temp5(r,:)>=0).*(Temp4(r,:)-Temp5(r,:));                            %次要中继的保密容量
           M_pu_R_sr_joint_sc(a,:)=Temp_pu_sc_2(a,:)+Temp_sr_sc(a,:);                                                                 %第i个主用户与第r个次要中继的保密容量之和
           a=a+1; 
        end     
    end
        
    %%%%求最大值及对应行数%%%%
    [Temp_2 Location_pu_sr]=max(M_pu_R_sr_joint_sc);                       %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp,Location_pu_sr记录最大值Temp所在行数
   
    %%%%确定max对应的最佳主用户%%%%
    Location_pu=zeros(1,L);
    for m=1:L
        if(mod(Location_pu_sr(1,m),R)~=0)
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R)+1;                 %max对应行数除以主用户数，若被整除，则第fix(Location_pu_sr(1,m)/R)+1个主用户为最佳主用户
        else
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R);                   %max对应行数除以主用户数，若不被整除，则第fix(Location_pu_sr(1,m)/R)个主用户为最佳主用户
        end
    end
            
    %%%%确定max对应的最佳次要中继%%%%
    Location_sr=rem(Location_pu_sr,R);                                     %max对应行数除以主用户数，若余数为零，则第R个次要中继为对应最佳次要中继
    Location_sr(Location_sr==0)=R;                                         %若余数不为零，则第rem(Location_pu_sr,M)
    
    %%%%分别计算max对应的最佳主用户及最佳次要中继的安全容量%%%%
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_2=zeros(1,L);
    Temp_sr_sc=zeros(1,L);
    for m=1:L
        Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbr(Location_sr(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2)+N0));         %Location_pu(1,m)表示第m次实验时第Location_pu(1,m)个主用户为最佳主用户 
        Temp_sr_sc(1,a)=0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbd(Location_sr(1,m),m))^2))/(Pp*abs(hid(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));
        Temp_pu_sc_2(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a);
        Temp_sr_sc(1,a)=(Temp_sr_sc(1,a)<0)*0+(Temp_sr_sc(1,a)>=0)*Temp_sr_sc(1,a);
        a=a+1;
    end
    
    Cs_M_pu_sc_2(k)=sum(Temp_pu_sc_2)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_R_sr_sc(k)=sum(Temp_sr_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    
    
    %%%%第一个时隙时最佳主用户及次要用户的保密容量%%%%
    %%%%注:主用户及次要中继已选定%%%%
    Temp6=zeros(1*N,L);
    Temp7=zeros(1*N,L);
    Temp8=zeros(N*1,L);
    Temp9=zeros(N*1,L);
    Temp_pu_sc_1=zeros(1*N,L);
    Temp_su_sc=zeros(1*N,L);
    M_pu_N_su_joint_sc=zeros(1*N,L);
    Ps=zeros(N,L);
    
    %%%%次要用户发射功率的限制%%%%
    for j=1:N
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m)^2)));
        end
    end
    
    %%%%根据确定好的最佳次要中继,进行N次搜索对次用户进行选取%%%%
    %%%%注:第m次实验的最佳主用户：Location_pu(1,m)%%%%
    %%%%注:第m次实验的最佳次要中继：Location_sr(1,m)%%%%
    
    a=1;
    for j=1:1:N
        for m=1:1:L
           Temp6(j,m)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hjr(j,m))^2)+N0));                            %主用户至主接收端的主信道容量
           Temp7(j,m)=log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hje(j,m))^2)+N0));                            %主用户的窃听信道容量
           Temp8(j,m)=0.5*log(1+(Ps(j,m)*(abs(hjb(j,m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0));                        %次用户至最佳次要中继的主信道容量
           Temp9(j,m)=0.5*log(1+(Ps(j,m)*(abs(hje(j,m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                        %次用户的窃听信道容量
   
           Temp_pu_sc_1(j,m)=(Temp6(j,m)-Temp7(j,m)<0).*0+(Temp6(j,m)-Temp7(j,m)>=0).*(Temp6(j,m)-Temp7(j,m));   %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(j,m)=(Temp8(j,m)-Temp9(j,m)<0).*0+(Temp8(j,m)-Temp9(j,m)>=0).*(Temp8(j,m)-Temp9(j,m));
           M_pu_N_su_joint_sc(j,m)= Temp_pu_sc_1(j,m)+ Temp_su_sc(j,m);      
        end
    end
         
    %%%%求最大值及对应行数%%%% 
    [Temp_1 Location_su]=max(M_pu_N_su_joint_sc);                          %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp_1,Location_su记录最大值Temp所在行数
    Cs_M_pu_N_su_joint_sc(k)=sum(Temp_1)/L;                                %MER为某一值时，主次用户保密容量之和的均值
    
    %%%%确定max对应的最佳次要用户%%%%      
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_1=zeros(1,L);
    Temp_su_sc=zeros(1,L);
    
    for m=1:L
           Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hjr(Location_su(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2)+N0));                            
           Temp_su_sc(1,a)=0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hjb(Location_su(1,m),m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                      
           Temp_pu_sc_1(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a); 
           Temp_su_sc(1,a)=(Temp_su_sc(1,a)<0)*0+(Temp_su_sc(1,a)>=0)*Temp_su_sc(1,a);
         a=a+1;
    end
        
    Cs_M_pu_sc_1(k)=sum(Temp_pu_sc_1)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_N_su_sc(k)=sum(Temp_su_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    %%%%最佳主用户安全容量与最佳次用户（或最佳次要中继）安全容量之和%%%%
    Cs_primary_sc(k)=min(Cs_M_pu_sc_1(k),Cs_M_pu_sc_2(k));                 %主用户安全容量取相邻两个时隙内的较小值
    Cs_secondary_sc(k)=min(Cs_R_sr_sc(k),Cs_N_su_sc(k));                   %次要网络安全容量取相邻两个时隙内次要用户与次要中继安全容量的较小值
    Cs_pri_second_joint_sc(k)=Cs_primary_sc(k)+Cs_secondary_sc(k);         %主网络与次网络安全容量之和
    
    clear hir hie hid hjd hje hjr
    k=k+1;
end




%%%%基于最佳中继的多用户调度%%%%
%%%%十个主用户,八个次用户,八个次要中继%%%%

%%%%各类参数（pu为主用户,su为次用户,sr为次要中继,r为主接收端,d为次接收端）%%%%
%%%%I=10dBm%%%%
%主用户参数
M=10;                  %主用户数（大于或等于2）
I=0.01;                 %主用户的干扰可容值  
Pp=10;                 %主用户功率  
Sigma_hir=0.6;         %pu->r,E(|hir|^2)=Sigma_ir
Sigma_hid=0.2;         %pu->d,E(|hid|^2)=Sigma_id
Sigma_hib=0.2;         %pu->sr,E(|hib|^2)=Sigma_ib
Theta_ir=1;
Theta_ie=0.7;

%次用户参数
N=8;                   %次用户数（大于或等于2）
P=10^0.2;              %次用户的分配功率P
Sigma_hjb=1;           %su->sr,E(|hjb|^2)=Sigma_jb
Sigma_hjr=0.2;         %su->r,E(|hjr|^2)=Sigma_jr
Theta_jb=1;
Theta_je=0.7;

%次要中继参数
R=8;                   %次要中继数（大于或等于2）
P=10^0.2;              %次要中继的分配功率P
Sigma_hbd=1;           %sr->d,E(|hbd|^2)=Sigma_bd
Sigma_hbr=0.2;         %sr->r,E(|hbr|^2)=Sigma_br
Theta_bd=1;
Theta_be=0.7;

%其他参数
L=10000;               %实验次数
N0=0.1;                %加性高斯白噪声功率
k=1;

for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio             
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir*Theta_ie/(Theta_ir*lambda);
    Sigma_hje=Sigma_hjb*Theta_je/(Theta_jb*lambda);
    Sigma_hbe=Sigma_hbd*Theta_be/(Theta_bd*lambda);

    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次要接收端的信道
    hib=sqrt(Sigma_hib/2)*randn(M,L)+sqrt(-Sigma_hib/2)*randn(M,L);        %主用户到次要中继的信道
      
    %次用户信道衰落系数%%%%
    hjb=sqrt(Sigma_hjb/2)*randn(N,L)+sqrt(-Sigma_hjb/2)*randn(N,L);        %次用户到最佳次要中继的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道
  
    %次要中继信道衰落系数%%%%
    hbd=sqrt(Sigma_hbd/2)*randn(R,L)+sqrt(-Sigma_hbd/2)*randn(R,L);        %次要中继到次要接收端的信道
    hbe=sqrt(Sigma_hbe/2)*randn(R,L)+sqrt(-Sigma_hbe/2)*randn(R,L);        %次要中继到窃听者的信道
    hbr=sqrt(Sigma_hbr/2)*randn(R,L)+sqrt(-Sigma_hbr/2)*randn(R,L);        %次要中继到主接收端的信道
    
    %%%%第二个时隙时主用户及次要中继的保密容量%%%%
    Temp2=zeros(M*R,L);
    Temp3=zeros(M*R,L);
    Temp4=zeros(R*M,L);
    Temp5=zeros(R*M,L);
    Temp_pu_sc_2=zeros(M*R,L);
    Temp_sr_sc=zeros(M*R,L);
    M_pu_R_sr_joint_sc=zeros(M*R,L);
    Pb=zeros(R,L);
   
    %%%%次要中继发射功率的限制%%%%
    for r=1:R
        for m=1:L 
            Pb(r,m)=min(P,I/(abs(hbr(r,m)^2)));                            %主要受P限制
        end
    end
        
    %%%%进行M*R次搜索对最佳次要中继进行选取%%%%
    a=1;
    for i=1:1:M
        for r=1:1:R
           Temp2(r,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Pb(r,:).*(abs(hbr(r,:)).^2)+N0));                                                   %主用户至主接收端的主信道容量
           Temp3(r,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Pb(r,:).*(abs(hbe(r,:)).^2)+N0));                                                   %主用户的窃听信道容量
           Temp4(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbd(r,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));                                               %次要中继至次接收端的主信道容量
           Temp5(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbe(r,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));                                               %次要中继的窃听信道容量
           Temp_pu_sc_2(a,:)=(Temp2(r,:)-Temp3(r,:)<0).*0+(Temp2(r,:)-Temp3(r,:)>=0).*(Temp2(r,:)-Temp3(r,:));                          %主用户的保密容量：若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_sr_sc(a,:)=(Temp4(r,:)-Temp5(r,:)<0).*0+(Temp4(r,:)-Temp5(r,:)>=0).*(Temp4(r,:)-Temp5(r,:));                            %次要中继的保密容量
           M_pu_R_sr_joint_sc(a,:)=Temp_pu_sc_2(a,:)+Temp_sr_sc(a,:);                                                                 %第i个主用户与第r个次要中继的保密容量之和
           a=a+1; 
        end     
    end
        
    %%%%求最大值及对应行数%%%%
    [Temp_2 Location_pu_sr]=max(M_pu_R_sr_joint_sc);                       %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp,Location_pu_sr记录最大值Temp所在行数
   
    %%%%确定max对应的最佳主用户%%%%
    Location_pu=zeros(1,L);
    for m=1:L
        if(mod(Location_pu_sr(1,m),R)~=0)
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R)+1;                 %max对应行数除以主用户数，若被整除，则第fix(Location_pu_sr(1,m)/R)+1个主用户为最佳主用户
        else
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R);                   %max对应行数除以主用户数，若不被整除，则第fix(Location_pu_sr(1,m)/R)个主用户为最佳主用户
        end
    end
            
    %%%%确定max对应的最佳次要中继%%%%
    Location_sr=rem(Location_pu_sr,R);                                     %max对应行数除以主用户数，若余数为零，则第R个次要中继为对应最佳次要中继
    Location_sr(Location_sr==0)=R;                                         %若余数不为零，则第rem(Location_pu_sr,M)
    
    %%%%分别计算max对应的最佳主用户及最佳次要中继的安全容量%%%%
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_2=zeros(1,L);
    Temp_sr_sc=zeros(1,L);
    for m=1:L
        Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbr(Location_sr(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2)+N0));         %Location_pu(1,m)表示第m次实验时第Location_pu(1,m)个主用户为最佳主用户 
        Temp_sr_sc(1,a)=0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbd(Location_sr(1,m),m))^2))/(Pp*abs(hid(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));
        Temp_pu_sc_2(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a);
        Temp_sr_sc(1,a)=(Temp_sr_sc(1,a)<0)*0+(Temp_sr_sc(1,a)>=0)*Temp_sr_sc(1,a);
        a=a+1;
    end
    
    Cs_M_pu_sc_2(k)=sum(Temp_pu_sc_2)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_R_sr_sc(k)=sum(Temp_sr_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    
    
    %%%%第一个时隙时最佳主用户及次要用户的保密容量%%%%
    %%%%注:主用户及次要中继已选定%%%%
    Temp6=zeros(1*N,L);
    Temp7=zeros(1*N,L);
    Temp8=zeros(N*1,L);
    Temp9=zeros(N*1,L);
    Temp_pu_sc_1=zeros(1*N,L);
    Temp_su_sc=zeros(1*N,L);
    M_pu_N_su_joint_sc=zeros(1*N,L);
    Ps=zeros(N,L);
    
    %%%%次要用户发射功率的限制%%%%
    for j=1:N
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m)^2)));
        end
    end
    
    %%%%根据确定好的最佳次要中继,进行N次搜索对次用户进行选取%%%%
    %%%%注:第m次实验的最佳主用户：Location_pu(1,m)%%%%
    %%%%注:第m次实验的最佳次要中继：Location_sr(1,m)%%%%
    
    a=1;
    for j=1:1:N
        for m=1:1:L
           Temp6(j,m)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hjr(j,m))^2)+N0));                            %主用户至主接收端的主信道容量
           Temp7(j,m)=log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hje(j,m))^2)+N0));                            %主用户的窃听信道容量
           Temp8(j,m)=0.5*log(1+(Ps(j,m)*(abs(hjb(j,m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0));                        %次用户至最佳次要中继的主信道容量
           Temp9(j,m)=0.5*log(1+(Ps(j,m)*(abs(hje(j,m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                        %次用户的窃听信道容量
   
           Temp_pu_sc_1(j,m)=(Temp6(j,m)-Temp7(j,m)<0).*0+(Temp6(j,m)-Temp7(j,m)>=0).*(Temp6(j,m)-Temp7(j,m));   %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(j,m)=(Temp8(j,m)-Temp9(j,m)<0).*0+(Temp8(j,m)-Temp9(j,m)>=0).*(Temp8(j,m)-Temp9(j,m));
           M_pu_N_su_joint_sc(j,m)= Temp_pu_sc_1(j,m)+ Temp_su_sc(j,m);      
        end
    end
         
    %%%%求最大值及对应行数%%%% 
    [Temp_1 Location_su]=max(M_pu_N_su_joint_sc);                          %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp_1,Location_su记录最大值Temp所在行数
    Cs_M_pu_N_su_joint_sc(k)=sum(Temp_1)/L;                                %MER为某一值时，主次用户保密容量之和的均值
    
    %%%%确定max对应的最佳次要用户%%%%      
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_1=zeros(1,L);
    Temp_su_sc=zeros(1,L);
    
    for m=1:L
           Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hjr(Location_su(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2)+N0));                            
           Temp_su_sc(1,a)=0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hjb(Location_su(1,m),m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                      
           Temp_pu_sc_1(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a); 
           Temp_su_sc(1,a)=(Temp_su_sc(1,a)<0)*0+(Temp_su_sc(1,a)>=0)*Temp_su_sc(1,a);
         a=a+1;
    end
        
    Cs_M_pu_sc_1(k)=sum(Temp_pu_sc_1)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_N_su_sc(k)=sum(Temp_su_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    %%%%最佳主用户安全容量与最佳次用户（或最佳次要中继）安全容量之和%%%%
    Cs_primary_sc_2(k)=min(Cs_M_pu_sc_1(k),Cs_M_pu_sc_2(k));                 %主用户安全容量取相邻两个时隙内的较小值
    Cs_secondary_sc_2(k)=min(Cs_R_sr_sc(k),Cs_N_su_sc(k));                   %次要网络安全容量取相邻两个时隙内次要用户与次要中继安全容量的较小值
    Cs_pri_second_joint_sc_2(k)=Cs_primary_sc(k)+Cs_secondary_sc(k);         %主网络与次网络安全容量之和
    
    clear hir hie hid hjd hje hjr
    k=k+1;
end

%%%%基于最佳中继的多用户调度%%%%
%%%%十个主用户,八个次用户,八个次要中继%%%%

%%%%各类参数（pu为主用户,su为次用户,sr为次要中继,r为主接收端,d为次接收端）%%%%
%%%%I=30dBm%%%%
%主用户参数
M=10;                  %主用户数（大于或等于2）
I=1;                   %主用户的干扰可容值  
Pp=10;                 %主用户功率  
Sigma_hir=0.6;         %pu->r,E(|hir|^2)=Sigma_ir
Sigma_hid=0.2;         %pu->d,E(|hid|^2)=Sigma_id
Sigma_hib=0.2;         %pu->sr,E(|hib|^2)=Sigma_ib
Theta_ir=1;
Theta_ie=0.7;

%次用户参数
N=8;                   %次用户数（大于或等于2）
P=10^0.2;              %次用户的分配功率P
Sigma_hjb=1;           %su->sr,E(|hjb|^2)=Sigma_jb
Sigma_hjr=0.2;         %su->r,E(|hjr|^2)=Sigma_jr
Theta_jb=1;
Theta_je=0.7;

%次要中继参数
R=8;                   %次要中继数（大于或等于2）
P=10^0.2;              %次要中继的分配功率P
Sigma_hbd=1;           %sr->d,E(|hbd|^2)=Sigma_bd
Sigma_hbr=0.2;         %sr->r,E(|hbr|^2)=Sigma_br
Theta_bd=1;
Theta_be=0.7;

%其他参数
L=10000;               %实验次数
N0=0.1;                %加性高斯白噪声功率
k=1;

for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio             
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir*Theta_ie/(Theta_ir*lambda);
    Sigma_hje=Sigma_hjb*Theta_je/(Theta_jb*lambda);
    Sigma_hbe=Sigma_hbd*Theta_be/(Theta_bd*lambda);

    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次要接收端的信道
    hib=sqrt(Sigma_hib/2)*randn(M,L)+sqrt(-Sigma_hib/2)*randn(M,L);        %主用户到次要中继的信道
      
    %次用户信道衰落系数%%%%
    hjb=sqrt(Sigma_hjb/2)*randn(N,L)+sqrt(-Sigma_hjb/2)*randn(N,L);        %次用户到最佳次要中继的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道
  
    %次要中继信道衰落系数%%%%
    hbd=sqrt(Sigma_hbd/2)*randn(R,L)+sqrt(-Sigma_hbd/2)*randn(R,L);        %次要中继到次要接收端的信道
    hbe=sqrt(Sigma_hbe/2)*randn(R,L)+sqrt(-Sigma_hbe/2)*randn(R,L);        %次要中继到窃听者的信道
    hbr=sqrt(Sigma_hbr/2)*randn(R,L)+sqrt(-Sigma_hbr/2)*randn(R,L);        %次要中继到主接收端的信道
    
    %%%%第二个时隙时主用户及次要中继的保密容量%%%%
    Temp2=zeros(M*R,L);
    Temp3=zeros(M*R,L);
    Temp4=zeros(R*M,L);
    Temp5=zeros(R*M,L);
    Temp_pu_sc_2=zeros(M*R,L);
    Temp_sr_sc=zeros(M*R,L);
    M_pu_R_sr_joint_sc=zeros(M*R,L);
    Pb=zeros(R,L);
   
    %%%%次要中继发射功率的限制%%%%
    for r=1:R
        for m=1:L 
            Pb(r,m)=min(P,I/(abs(hbr(r,m)^2)));                            %主要受P限制
        end
    end
        
    %%%%进行M*R次搜索对最佳次要中继进行选取%%%%
    a=1;
    for i=1:1:M
        for r=1:1:R
           Temp2(r,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Pb(r,:).*(abs(hbr(r,:)).^2)+N0));                                                   %主用户至主接收端的主信道容量
           Temp3(r,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Pb(r,:).*(abs(hbe(r,:)).^2)+N0));                                                   %主用户的窃听信道容量
           Temp4(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbd(r,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));                                               %次要中继至次接收端的主信道容量
           Temp5(r,:)=0.5*log(1+(Pb(r,:).*(abs(hbe(r,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));                                               %次要中继的窃听信道容量
           Temp_pu_sc_2(a,:)=(Temp2(r,:)-Temp3(r,:)<0).*0+(Temp2(r,:)-Temp3(r,:)>=0).*(Temp2(r,:)-Temp3(r,:));                          %主用户的保密容量：若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_sr_sc(a,:)=(Temp4(r,:)-Temp5(r,:)<0).*0+(Temp4(r,:)-Temp5(r,:)>=0).*(Temp4(r,:)-Temp5(r,:));                            %次要中继的保密容量
           M_pu_R_sr_joint_sc(a,:)=Temp_pu_sc_2(a,:)+Temp_sr_sc(a,:);                                                                 %第i个主用户与第r个次要中继的保密容量之和
           a=a+1; 
        end     
    end
        
    %%%%求最大值及对应行数%%%%
    [Temp_2 Location_pu_sr]=max(M_pu_R_sr_joint_sc);                       %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp,Location_pu_sr记录最大值Temp所在行数
   
    %%%%确定max对应的最佳主用户%%%%
    Location_pu=zeros(1,L);
    for m=1:L
        if(mod(Location_pu_sr(1,m),R)~=0)
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R)+1;                 %max对应行数除以主用户数，若被整除，则第fix(Location_pu_sr(1,m)/R)+1个主用户为最佳主用户
        else
            Location_pu(1,m)=fix(Location_pu_sr(1,m)/R);                   %max对应行数除以主用户数，若不被整除，则第fix(Location_pu_sr(1,m)/R)个主用户为最佳主用户
        end
    end
            
    %%%%确定max对应的最佳次要中继%%%%
    Location_sr=rem(Location_pu_sr,R);                                     %max对应行数除以主用户数，若余数为零，则第R个次要中继为对应最佳次要中继
    Location_sr(Location_sr==0)=R;                                         %若余数不为零，则第rem(Location_pu_sr,M)
    
    %%%%分别计算max对应的最佳主用户及最佳次要中继的安全容量%%%%
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_2=zeros(1,L);
    Temp_sr_sc=zeros(1,L);
    for m=1:L
        Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbr(Location_sr(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2)+N0));         %Location_pu(1,m)表示第m次实验时第Location_pu(1,m)个主用户为最佳主用户 
        Temp_sr_sc(1,a)=0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbd(Location_sr(1,m),m))^2))/(Pp*abs(hid(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Pb(Location_sr(1,m),m)*(abs(hbe(Location_sr(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));
        Temp_pu_sc_2(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a);
        Temp_sr_sc(1,a)=(Temp_sr_sc(1,a)<0)*0+(Temp_sr_sc(1,a)>=0)*Temp_sr_sc(1,a);
        a=a+1;
    end
    
    Cs_M_pu_sc_2(k)=sum(Temp_pu_sc_2)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_R_sr_sc(k)=sum(Temp_sr_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    
    
    %%%%第一个时隙时最佳主用户及次要用户的保密容量%%%%
    %%%%注:主用户及次要中继已选定%%%%
    Temp6=zeros(1*N,L);
    Temp7=zeros(1*N,L);
    Temp8=zeros(N*1,L);
    Temp9=zeros(N*1,L);
    Temp_pu_sc_1=zeros(1*N,L);
    Temp_su_sc=zeros(1*N,L);
    M_pu_N_su_joint_sc=zeros(1*N,L);
    Ps=zeros(N,L);
    
    %%%%次要用户发射功率的限制%%%%
    for j=1:N
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m)^2)));
        end
    end
    
    %%%%根据确定好的最佳次要中继,进行N次搜索对次用户进行选取%%%%
    %%%%注:第m次实验的最佳主用户：Location_pu(1,m)%%%%
    %%%%注:第m次实验的最佳次要中继：Location_sr(1,m)%%%%
    
    a=1;
    for j=1:1:N
        for m=1:1:L
           Temp6(j,m)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hjr(j,m))^2)+N0));                            %主用户至主接收端的主信道容量
           Temp7(j,m)=log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(j,m)*(abs(hje(j,m))^2)+N0));                            %主用户的窃听信道容量
           Temp8(j,m)=0.5*log(1+(Ps(j,m)*(abs(hjb(j,m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0));                        %次用户至最佳次要中继的主信道容量
           Temp9(j,m)=0.5*log(1+(Ps(j,m)*(abs(hje(j,m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                        %次用户的窃听信道容量
   
           Temp_pu_sc_1(j,m)=(Temp6(j,m)-Temp7(j,m)<0).*0+(Temp6(j,m)-Temp7(j,m)>=0).*(Temp6(j,m)-Temp7(j,m));   %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(j,m)=(Temp8(j,m)-Temp9(j,m)<0).*0+(Temp8(j,m)-Temp9(j,m)>=0).*(Temp8(j,m)-Temp9(j,m));
           M_pu_N_su_joint_sc(j,m)= Temp_pu_sc_1(j,m)+ Temp_su_sc(j,m);      
        end
    end
         
    %%%%求最大值及对应行数%%%% 
    [Temp_1 Location_su]=max(M_pu_N_su_joint_sc);                          %M*R次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值Temp_1,Location_su记录最大值Temp所在行数
    Cs_M_pu_N_su_joint_sc(k)=sum(Temp_1)/L;                                %MER为某一值时，主次用户保密容量之和的均值
    
    %%%%确定max对应的最佳次要用户%%%%      
    a=1;
    Temp_pu_sc=zeros(1,L);
    Temp_pu_sc_1=zeros(1,L);
    Temp_su_sc=zeros(1,L);
    
    for m=1:L
           Temp_pu_sc(1,a)=log(1+(Pp*abs(hir(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hjr(Location_su(1,m),m))^2)+N0))-log(1+(Pp*abs(hie(Location_pu(1,m),m))^2)/(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2)+N0));                            
           Temp_su_sc(1,a)=0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hjb(Location_su(1,m),m))^2))/(Pp*abs(hib(Location_pu(1,m),m))^2+N0))-0.5*log(1+(Ps(Location_su(1,m),m)*(abs(hje(Location_su(1,m),m))^2))/(Pp*abs(hie(Location_pu(1,m),m))^2+N0));                      
           Temp_pu_sc_1(1,a)=(Temp_pu_sc(1,a)<0)*0+(Temp_pu_sc(1,a)>=0)*Temp_pu_sc(1,a); 
           Temp_su_sc(1,a)=(Temp_su_sc(1,a)<0)*0+(Temp_su_sc(1,a)>=0)*Temp_su_sc(1,a);
         a=a+1;
    end
        
    Cs_M_pu_sc_1(k)=sum(Temp_pu_sc_1)/L;                                   %MER为某一值时，最佳主用户的安全容量均值
    Cs_N_su_sc(k)=sum(Temp_su_sc)/L;                                       %MER为某一值时，最佳次要中继的安全容量均值
 
    %%%%最佳主用户安全容量与最佳次用户（或最佳次要中继）安全容量之和%%%%
    Cs_primary_sc_3(k)=min(Cs_M_pu_sc_1(k),Cs_M_pu_sc_2(k));                 %主用户安全容量取相邻两个时隙内的较小值
    Cs_secondary_sc_3(k)=min(Cs_R_sr_sc(k),Cs_N_su_sc(k));                   %次要网络安全容量取相邻两个时隙内次要用户与次要中继安全容量的较小值
    Cs_pri_second_joint_sc_3(k)=Cs_primary_sc(k)+Cs_secondary_sc(k);         %主网络与次网络安全容量之和
    
    clear hir hie hid hjd hje hjr
    k=k+1;
end

MER=-10:2:30;
plot(MER,Cs_secondary_sc_3,'k-v');
hold on
plot(MER,Cs_secondary_sc,'k-*');
hold on
plot(MER,Cs_secondary_sc_2,'k-o');
hold on
set(gcf,'color','white');
xlabel('MER(dB)')
ylabel('次用户安全容量(bit/s/Hz)')
legend('location','NorthWest','I=30dBm (M=10 N=8 R=8)','I=20dBm (M=10 N=8 R=8)','I=10dBm (M=10 N=8 R=8)')