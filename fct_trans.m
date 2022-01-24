clear all; close all; clc;

%------------------------------Cindy Iskandar------------------------------
%----------------------------------161862----------------------------------
%---------------------------Responsable: R. Ghosn--------------------------
%--------------------------------TC MAS 2019-------------------------------

%% Parametres de la MAS 

Ns=3000;                                                %en tr/min
Nn=2954;                                                %en tr/min
Pn= 90e3;                                               %puissance mecanique sur l'arbre en W
J=0.51;                                                 %en Kg.m2
Wn=2*pi*Nn/60;                                          %en rad/s
Cn=Pn/Wn;                                               %en N.m
P=1;                                                    %Nombre de paires de Poles
fn=50;                                                  %en Hz
Un=400;                                                 %en V
    
%facteur de puissance
cos_phi_n_1=0.91;
cos_phi_n_2=0.89;
cos_phi_n_3=0.85;
    
%rendement de la MAS
R_n_1=0.949;
R_n_2=0.953;
R_n_3=0.952;

In=(Pn/(R_n_1*sqrt(3)*Un*cos_phi_n_1));                 %en A
Id=7.3*In;                                              %en A
Cd=2.4*Cn;                                              %en N.m
Cm=2.3*Cn;                                              %en N.m

alpha=0.015;
f=(alpha/100)*Pn/(Wn)^2;                                %coef de frottements en N.m.s/rad

wn=2*pi*fn;                                             %pulsation en rad/s
Vn=Un/sqrt(3);                                          %en V
Wb=wn/P;                                                %Vitesse de base en rad/s
Zn=Vn/In;                                               %en Ohm

gn=(Ns-Nn)/Ns;                                          %adimensionne

phi_n=acos(cos_phi_n_1);
sin_phi_n_1=sqrt(1-(cos_phi_n_1)^2);

kc=0;                                                   %On est a vide
%% Calcul des elements du schema monophase etoile eq

%Calcul de N2w
syms N2w;
N2w=solve((3*(Vn)^2)/(wn*2*N2w)==Cm,N2w);
N2w=double(N2w);
nw=(N2w/Zn);

%Calcul de R2
a=Cn*wn;
b=-(Un^2);
c=N2w^2;
p=[a b c];
r=roots(p);

if ((r(1)*gn/Zn >0.01) && (r(1)*gn/Zn < 0.05))
    R2=gn*r(1);
end 

if ((r(2)*gn/Zn >0.01) && (r(2)*gn/Zn < 0.05))
    R2=gn*r(2);
end

r2=(R2/Zn);

%Calcul de L1w et de Rpfe
Yn=(1/Zn)*exp(-1j*phi_n);
comp=Yn-(R2/gn-1j*N2w)/((R2/gn)^2+(N2w)^2);

Rpfe=double(1/real(comp));
rpfe=(Rpfe/Zn);
L1w=double(-1/imag(comp));
lw=(L1w/Zn);

%Calcul de Ls, Lr, Msr et sigma
Ls=(N2w/2+L1w)/wn;
Lr=(N2w/2+L1w)/wn;
Msr=L1w/wn;
sigma=1-(Msr^2)/(Ls*Lr);

%Calcul de Rs et Rr
Rr=double(R2);
Rs=0.9*Rr;
%% Boucle de flux 

Ts=Ls/Rs;
Tr=Lr/Rr;
p=tf('p');

H_phi_BO=(1/Rs)/(1+(Ts+Tr)*p+(sigma*Ts*Tr)*p^2);

%Un correcteur dimensionnee prend la forme de :
k_phi=1/(0.3*Tr);
poles=pole(H_phi_BO);
[zphi,pphi,kphi]=zpkdata(H_phi_BO);
for i=1:length(pole(H_phi_BO))
    if (abs(poles(i))<abs(poles(1)))
        C_phi=k_phi*(p-poles(i))/(kphi*p);
    else
        C_phi=k_phi*(p+poles(1))/(kphi*p);
    end
end
zero_phi=zero(C_phi);

%la fonction de transfert approximee s'ecrit de la forme:
H_phi_app=kphi/(p-zero_phi);
%temps de montee de l'ordre de Tr;

%% Boucle de couple

%phi_s_n_abc=Vn/wn;                                 %dans abc et on veut passer en dq
%phi_s_n=sqrt(3/2)*phi_s_n_abc*sqrt(2);             %on multiplie par sqrt(2) pour obtenir la valeur efficace
%donc on aura: phi_s_n=Un/wn donc toute machine ayant meme wn et meme Un on
%meme phi_s_n dans le repere tournant dq.

phi_s_n=Un/wn;
phi_r_n=phi_s_n;                                    %on neglige les pertes
phi_rd_n=phi_r_n;                                   %car on a aligne l'axe d avec celui du flux par commande vectorielle
Imr=phi_rd_n/Msr;
T_phi_BF=0.3*Tr;
H_couple_BO=P*Msr^2*Imr/(Lr*(Rs+sigma*Ls*p)*(1+T_phi_BF*p));

%On trouve H_couple_BO = 1311.2/(p+146.8)*(p+33.3)
%on elimine le pole le plus lent (le plus loin de 0):
%H_phi_BO_3=1311.2/(p+33.3);

k_couple=1/(0.3*sigma*Ts);                          %temps de montee de l'ordre de sigma*Ts
poles=pole(H_couple_BO);
[zcouple,pcouple,kcouple]=zpkdata(H_couple_BO);
for i=1:length(pole(H_couple_BO))
    if (abs(poles(i))<abs(poles(1)))
        C_couple=k_couple*(p-poles(i))/(kcouple*p);
    else
        C_couple=k_couple*(p+poles(1))/(kcouple*p);
    end
end
zero_couple=zero(C_couple);

%la fonction de transfert approximee s'ecrit de la forme:
H_couple_app=kcouple/(p-zero_couple);

%le temps de montee doit etre de l'ordre de sigma*Ts
%la presence d'un depassement est normal

%% Boucle de Vitesse

T_C_BF=0.0445;                                      %temps de montee de la courbe du couple
H_omega_BO=1/((J*p+f)*(1+T_C_BF*p));

k_omega=1/(0.3*0.5);
poles=pole(H_omega_BO);
[zomega,pomega,komega]=zpkdata(H_omega_BO);
for i=1:length(pole(H_omega_BO))
    if (abs(poles(i))<abs(poles(1)))
        C_omega=k_omega*(p-poles(i))/(komega*p);
    else
        C_omega=k_omega*(p+poles(1))/(komega*p);
    end
end
zero_omega=zero(C_omega);

%la fonction de transfert approximee s'ecrit de la forme:
H_omega_app=komega/(p-zero_omega);
%temps de montee egal a 0.45 s il doit d'etre de l'ordre de 0.5 s

T_s=10;
%% Figures

sim('Fct_app.slx');
figure()
subplot(3,1,1)
plot(temps,flux,'linewidth',2),grid;
hold on
plot(temps,flux_ref,'linewidth',2);
axis([0.9 4 0 1.25]);
xlabel('Temps (s)');
ylabel('Imr (A)');
title('Variation du courant Imr en fonction du temps');
legend('Imr','Imr^*');
subplot(3,1,2)
plot(temps,couple,'linewidth',2),grid;
hold on
plot(temps,couple_ref,'linewidth',2);
axis([0.9 4 0 1.25]);
xlabel('Temps (s)');
ylabel('Couple (N.m)');
title('Variation du couple en fonction du temps');
legend('Couple','Couple Ref');
subplot(3,1,3)
plot(temps,omega,'linewidth',2),grid;
hold on
plot(temps,omega_ref,'linewidth',2);
axis([0.9 4 0 1.25]);
xlabel('Temps (s)');
ylabel('Omega (rad/s)');
title('Variation de la vitesse mecanique en fonction du temps');
legend('Omega','Omega Ref');
