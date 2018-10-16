clear
clc

%%%%各类参数（p为主用户,s为次用户,r为主接收端,d为次接收端）%%%%
%%%%四个主用户与三个次用户%%%%
%主用户参数
Sigma_hir=0.6;         %p->r,E(|hir|^2)=Sigma_ir
Sigma_hid=0.2;         %p->d,E(|hid|^2)=Sigma_id
Theta_ir=1;
Theta_ie=0.7;
M=4;                   %主用户数（大于或等于2）
I=0.1;                 %主用户的干扰可容值  
Pp=10;                 %主用户功率P 

%次用户参数
Sigma_hjd=0.6;         %s->d,E(|hjd|^2)=Sigma_jd
Sigma_hjr=0.2;         %s->r,E(|hjr|^2)=Sigma_jr
Theta_jd=1;
Theta_je=0.7;
N=3;                   %次用户数（大于或等于2）
P=10^0.2;              %次用户的分配功率P

%其他参数
L=10000;               %实验次数
N0=0.1;                %加性高斯白噪声功率
k=1;

for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir*Theta_ie/(Theta_ir*lambda);
    Sigma_hje=Sigma_hjd*Theta_je/(Theta_jd*lambda);
    
    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次接收端的信道
    
    %%%%次用户信道衰落系数%%%%
    hjd=sqrt(Sigma_hjd/2)*randn(N,L)+sqrt(-Sigma_hjd/2)*randn(N,L);        %次用户到次接收端的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道

   
    %%%%M个主用户和对应的N个次用户的保密容量之和%%%%
    Temp2=zeros(M*N,L);
    Temp3=zeros(M*N,L);
    Temp4=zeros(N*M,L);
    Temp5=zeros(N*M,L);
    Temp_pu_sc=zeros(M*N,L);
    Temp_su_sc=zeros(M*N,L);
    M_pu_N_su_joint_sc=zeros(M*N,L);
    Ps=zeros(N,L);
    
    %%%%次用户的发射功率限制%%%%
    for j=1:N 
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m))^2));
        end
    end
    
    a=1;
    %%%%进行M*N次搜索%%%%
    for i=1:1:M
         for j=1:1:N
           Temp2(j,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Ps(j,:).*(abs(hjr(j,:)).^2)+N0));     %主用户的主信道容量
           Temp3(j,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Ps(j,:).*(abs(hje(j,:)).^2)+N0));     %主用户的窃听信道容量
           Temp4(j,:)=log(1+(Ps(j,:).*(abs(hjd(j,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));     %次用户的主信道容量
           Temp5(j,:)=log(1+(Ps(j,:).*(abs(hje(j,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));     %次用户的窃听信道容量
   
           Temp_pu_sc(a,:)=(Temp2(j,:)-Temp3(j,:)<0).*0+(Temp2(j,:)-Temp3(j,:)>=0).*(Temp2(j,:)-Temp3(j,:));                  %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(a,:)=(Temp4(j,:)-Temp5(j,:)<0).*0+(Temp4(j,:)-Temp5(j,:)>=0).*(Temp4(j,:)-Temp5(j,:));
           M_pu_N_su_joint_sc(a,:)= Temp_pu_sc(a,:)+ Temp_su_sc(a,:); 
           a=a+1;
         end
        
    end
        
    Temp=max(M_pu_N_su_joint_sc);                                          %M*N次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值
    Cs_M_pu_N_su(k)=sum(Temp)/L;                                           %MER为某一值时，主次用户保密容量之和的均值
    
    clear hir hie hid hjd hje hjr
    
    %%%%单个主用户与单个次用户的保密容量之和%%%% 
    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(1,L)+sqrt(-Sigma_hir/2)*randn(1,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(1,L)+sqrt(-Sigma_hie/2)*randn(1,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(1,L)+sqrt(-Sigma_hid/2)*randn(1,L);        %主用户到次接收端的信道
    
    %%%%次用户信道衰落系数%%%%
    hjd=sqrt(Sigma_hjd/2)*randn(1,L)+sqrt(-Sigma_hjd/2)*randn(1,L);        %次用户到次接收端的信道
    hje=sqrt(Sigma_hje/2)*randn(1,L)+sqrt(-Sigma_hje/2)*randn(1,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(1,L)+sqrt(-Sigma_hjr/2)*randn(1,L);        %次用户到主接收端的信道

      Cpm=zeros(1,L);
      Cpe=zeros(1,L);
      Temp2=0;
      Temp3=0;
      Ps=zeros(1,L);
      for j=1:L
           Ps(1,j)=min(P,I/(abs(hjr(1,j))^2));
      end
      for i=1:1:L         
          Cpm_one=log2(1+(Pp*abs(hir(1,i))^2)/(Ps(1,i).*abs(hjd(1,i))^2+N0));
          Cpe_one=log2(1+(Pp*abs(hie(1,i))^2)/(Ps(1,i).*abs(hje(1,i))^2+N0));
          Temp2=Temp2+([Cpm_one-Cpe_one<0].*0+[Cpm_one-Cpe_one].*[Cpm_one-Cpe_one>=0]);
          Csm_one=log2(1+(Ps(1,i).*abs(hjd(1,i))^2)/(Pp*abs(hid(1,i))^2+N0));
          Cse_one=log2(1+(Ps(1,i).*abs(hje(1,i))^2)/(Pp*abs(hie(1,i))^2+N0));
          Temp3=Temp3+([Csm_one-Cse_one<0].*0+[Csm_one-Cse_one].*[Csm_one-Cse_one>=0]);
      end
      Cs_one_pu_one_su(k)=(Temp2+Temp3)/L;  
      k=k+1;
end   


%%%%六个主用户与四个次用户%%%%
%主用户参数
M=6;                   %主用户数（大于或等于2）

%次用户参数
N=4;                   %次用户数（大于或等于2）

k=1;
for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir/lambda;
    Sigma_hje=Sigma_hjd/lambda;
    
    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次接收端的信道
    
    %%%%次用户信道衰落系数%%%%
    hjd=sqrt(Sigma_hjd/2)*randn(N,L)+sqrt(-Sigma_hjd/2)*randn(N,L);        %次用户到次接收端的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道

   
    %%%%M个主用户和对应的N个次用户的保密容量之和%%%%
    Temp2=zeros(M*N,L);
    Temp3=zeros(M*N,L);
    Temp4=zeros(N*M,L);
    Temp5=zeros(N*M,L);
    Temp_pu_sc=zeros(M*N,L);
    Temp_su_sc=zeros(M*N,L);
    M_pu_N_su_joint_sc=zeros(M*N,L);
    Ps=zeros(N,L);
    
    %%%%次用户的发射功率限制%%%%
    for j=1:N 
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m))^2));
        end
    end
    
    a=1;
    %%%%进行M*N次搜索%%%%
    for i=1:1:M
         for j=1:1:N
           Temp2(j,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Ps(j,:).*(abs(hjr(j,:)).^2)+N0));     %主用户的主信道容量
           Temp3(j,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Ps(j,:).*(abs(hje(j,:)).^2)+N0));     %主用户的窃听信道容量
           Temp4(j,:)=log(1+(Ps(j,:).*(abs(hjd(j,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));     %次用户的主信道容量
           Temp5(j,:)=log(1+(Ps(j,:).*(abs(hje(j,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));     %次用户的窃听信道容量
   
           Temp_pu_sc(a,:)=(Temp2(j,:)-Temp3(j,:)<0).*0+(Temp2(j,:)-Temp3(j,:)>=0).*(Temp2(j,:)-Temp3(j,:));                  %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(a,:)=(Temp4(j,:)-Temp5(j,:)<0).*0+(Temp4(j,:)-Temp5(j,:)>=0).*(Temp4(j,:)-Temp5(j,:));
           M_pu_N_su_joint_sc(a,:)= Temp_pu_sc(a,:)+ Temp_su_sc(a,:); 
           a=a+1;
         end
        
    end
        
    Temp=max(M_pu_N_su_joint_sc);                                          %M*N次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值
    Cs_M_pu_N_su_2(k)=sum(Temp)/L;                                         %MER为某一值时，主次用户保密容量之和的均值
    k=k+1;
end   

%%%%十个主用户与八个次用户%%%%
%主用户参数
M=10;                   %主用户数（大于或等于2）

%次用户参数
N=8;                   %次用户数（大于或等于2）

k=1;
for MER=-10:2:30                                                           %MER：main-to-eavesdropper ratio
    lambda=10^(MER/10);                                                    %dB化为十进制数
    Sigma_hie=Sigma_hir/lambda;
    Sigma_hje=Sigma_hjd/lambda;
    
    %%%%主用户信道衰落系数%%%%
    hir=sqrt(Sigma_hir/2)*randn(M,L)+sqrt(-Sigma_hir/2)*randn(M,L);        %主用户到主接收端的信道
    hie=sqrt(Sigma_hie/2)*randn(M,L)+sqrt(-Sigma_hie/2)*randn(M,L);        %主用户到窃听者的信道
    hid=sqrt(Sigma_hid/2)*randn(M,L)+sqrt(-Sigma_hid/2)*randn(M,L);        %主用户到次接收端的信道
    
    %%%%次用户信道衰落系数%%%%
    hjd=sqrt(Sigma_hjd/2)*randn(N,L)+sqrt(-Sigma_hjd/2)*randn(N,L);        %次用户到次接收端的信道
    hje=sqrt(Sigma_hje/2)*randn(N,L)+sqrt(-Sigma_hje/2)*randn(N,L);        %次用户到窃听者的信道
    hjr=sqrt(Sigma_hjr/2)*randn(N,L)+sqrt(-Sigma_hjr/2)*randn(N,L);        %次用户到主接收端的信道

   
    %%%%M个主用户和对应的N个次用户的保密容量之和%%%%
    Temp2=zeros(M*N,L);
    Temp3=zeros(M*N,L);
    Temp4=zeros(N*M,L);
    Temp5=zeros(N*M,L);
    Temp_pu_sc=zeros(M*N,L);
    Temp_su_sc=zeros(M*N,L);
    M_pu_N_su_joint_sc=zeros(M*N,L);
    Ps=zeros(N,L);
    
    %%%%次用户的发射功率限制%%%%
    for j=1:N 
        for m=1:L 
            Ps(j,m)=min(P,I/(abs(hjr(j,m))^2));
        end
    end
    
    a=1;
    %%%%进行M*N次搜索%%%%
    for i=1:1:M
         for j=1:1:N
           Temp2(j,:)=log(1+(Pp*abs(hir(i,:)).^2)./(Ps(j,:).*(abs(hjr(j,:)).^2)+N0));     %主用户的主信道容量
           Temp3(j,:)=log(1+(Pp*abs(hie(i,:)).^2)./(Ps(j,:).*(abs(hje(j,:)).^2)+N0));     %主用户的窃听信道容量
           Temp4(j,:)=log(1+(Ps(j,:).*(abs(hjd(j,:)).^2))./(Pp*abs(hid(i,:)).^2+N0));     %次用户的主信道容量
           Temp5(j,:)=log(1+(Ps(j,:).*(abs(hje(j,:)).^2))./(Pp*abs(hie(i,:)).^2+N0));     %次用户的窃听信道容量
   
           Temp_pu_sc(a,:)=(Temp2(j,:)-Temp3(j,:)<0).*0+(Temp2(j,:)-Temp3(j,:)>=0).*(Temp2(j,:)-Temp3(j,:));                  %若主信道与窃听信道的信道容量差值小于零则结果作为零处理
           Temp_su_sc(a,:)=(Temp4(j,:)-Temp5(j,:)<0).*0+(Temp4(j,:)-Temp5(j,:)>=0).*(Temp4(j,:)-Temp5(j,:));
           M_pu_N_su_joint_sc(a,:)= Temp_pu_sc(a,:)+ Temp_su_sc(a,:); 
           a=a+1;
         end
        
    end
        
    Temp=max(M_pu_N_su_joint_sc);                                          %M*N次搜索出的最大值,在M_pu_N_su_joint_sc的每列找出最大值
    Cs_M_pu_N_su_3(k)=sum(Temp)/L;                                         %MER为某一值时，主次用户保密容量之和的均值
    k=k+1;
end  

MER=-10:2:30;
plot(MER, Cs_M_pu_N_su_3,'k-d');
hold on
plot(MER, Cs_M_pu_N_su_2,'k-o');
hold on
plot(MER, Cs_M_pu_N_su,'k-v');
hold on
plot(MER, Cs_one_pu_one_su,'k-*');
set(gcf,'color','white');
xlabel('MER(dB)')
ylabel('主次用户安全容量之和(bit/s/Hz)')
legend('location','southeast','M=10 N=8','M=6 N=4','M=4 N=3','M=1 N=1')