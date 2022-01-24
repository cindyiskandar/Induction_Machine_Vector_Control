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

%% Modelisation de la MAS

A=sqrt(2/3)*[1 -0.5 -0.5; 0 sqrt(3)/2 -sqrt(3)/2; sqrt(0.5) sqrt(0.5) sqrt(0.5)];

Ai=inv(A);

B=[Ls 0 Msr 0; 0 Ls 0 Msr; Msr 0 Lr 0; 0 Msr 0 Lr];

Binv=inv(B);

Bi = [1/(sigma*Ls) 0 -Msr/(sigma*Ls*Lr) 0; 0 1/(sigma*Ls) 0 -Msr/(sigma*Ls*Lr); -Msr/(sigma*Ls*Lr) 0 1/(sigma*Ls) 0; 0 -Msr/(sigma*Ls*Lr) 0 1/(sigma*Ls) ];

%% Intro 

disp('Nous allons effectuer la commande vectorielle sur la machine asynchrone FLSES 280 M');
disp('Calculons les parametres du schema monophase etoile equivalent relatif a la machine.');
disp('       ');
disp('_____________________________________________________________________');
disp('       ');

%% Affichage des parametres de la MAS

disp('     ');
disp('Les Valeurs reelles des elements du Schema monophase etoile equivalent sont les suivantes:  ');
fprintf('Rpfe = %0.5f Ohms\n', Rpfe);
fprintf('L1w = %0.5f Ohms\n', L1w);
fprintf('N2w = %0.5f Ohms\n', N2w);
fprintf('R2 = %0.5f Ohms\n', R2);
fprintf('Rs = %0.5f Ohms\n', Rs);
disp('     ');
disp('Les Valeurs reduites des elements du schema monophase etoile equivalent sont les suivantes:  ');
fprintf('rpfe = %0.5f pu\n', rpfe);
fprintf('lw = %0.5f pu\n', lw);
fprintf('nw = %0.5f pu\n', nw);
fprintf('r2 = %0.5f pu\n',r2);
fprintf('rs = %0.5f pu\n', Rs/Zn);
disp('     ');
disp('Les valeurs des inductances cycliques referees au stator');
fprintf('Ls = %0.5f H\n', Ls);
fprintf('Lr = %0.5f H\n', Lr);
fprintf('Msr = %0.5f H\n', Msr);
disp('     ');
disp('La valeur du coefficient de dispersion');
fprintf('sigma = %0.5f \n', sigma);
disp ('      ');
disp('_____________________________________________________________________');
disp('       ');
disp('Choisissez le temps de simulation que vous voulez.');
disp('       ');
disp('Il est preferable que le temps de simulation soit plus grand que 7 secs afin de pouvoir observer le regime transitoire et permanent de la machine.');
T_s=input('Ts = ');
while(isempty(T_s) || T_s<0 || T_s>60)
    T_s=input('Ts = ');
end 
disp('     ');
disp('_____________________________________________________________________');
disp('     ');

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

%% Resistance Rotorique Variable

alfa=3.98e-3;                               %coefficient de temperature du cuivre (en K^-1)
delta_T=80;                                 %en K
to=1;                                       %Pour atteindre le regime permanent en 5s
R0=Rr/(1+alfa*delta_T);                     %Resistance a froid

%% Onduleur modelise par un gain

epsn=5;
ktachy=90/157;
kshunt=1/240;
kom=epsn/(Wn*ktachy);
kcurr=epsn/(In*kshunt);
ktension=epsn/Vn;

%% Onduleur 3ph

gam=[2 -1 -1; -1 2 -1; -1 -1 2];


%% Couple Resistant

cte=1;
disp('Voulez vous : ');
disp('        ');
disp('1. Simuler la machine en Boucle Ouverte ? ');
disp('        ');
disp('2. Simuler la machine en BF et avec commande vectorielle ?');
disp('        ');
disp('3. Simuler la machine en BF et avec onduleur ?');
disp('        ');
mot1=input('Reponse : ');
while(isempty(mot1)|| (mot1~=1 && mot1~=2 && mot1~=3))
    disp('        ');
    mot1=input('Reponse : ');
end 
disp('        ');
disp('_____________________________________________________________________');
while(cte==1)
        disp('     ');
        disp('Choisissez la forme de votre couple : ');
        disp('     ');
        disp('Voulez vous : ');
        disp('     ');
        disp('1. Un couple resistant proportionnel au carre de la vitesse ? ');
        disp('     ');
        disp('2. Un couple resistant proportionnel a la vitesse ? ');
        disp('     ');
        disp('3. Un couple resistant constant ? ');
        disp('     ');
        disp('4. Un couple resistant inversement proportionnel a la vitesse ? ');
        disp('     ');
        disp('5. Un couple resistant quelconque (prop. au carre de la vitesse / prop. a la vitesse / cte / inv. prop. a la vitesse)  ?');
        disp('     ');
        q=input('Reponse : ');
        while(isempty(q)|| (q~=1 && q~=2 && q~=3 && q~=4 && q~=5))
            disp('     ');
            q=input('Reponse : '); 
        end 
        a1=0;b1=0;c1=0;d1=0;
        if (q==5)
            q=randi(4);
            disp('      ');
        end 
        if (q==1)
            aa1=(Cn/(Wn)^2);
            disp('Sachant que Cr = k * (omega)^2 : ');
            disp('        ');
            fprintf('Inserer k sachant qu''il n''excede pas : %0.5f SI ',aa1);
            disp('       ');
            a1=input('Reponse : ');
            while (a1>aa1 || isempty(a1) || a1<0)
               disp('       ');
               a1=input('Reponse : ');
            end
        end
        if q==2
            aa2=(Cn/(Wn)-f);
            disp('Sachant que Cr = k * (omega) : ');
            disp('        ');
            fprintf('Inserer k sachant qu''il n''excede pas : %0.5f N.m.s/rad ',aa2);
            disp('       ');
            b1=input('Reponse : ');
            while (b1>aa2 || isempty(b1) || b1<0)
               disp('       ');
               b1=input('Reponse : ');
            end
        end
        if q==3
           aa3=Cn;
           disp('Sachant que Cr = k : ');
           disp('      ');
           fprintf('Inserer k sachant qu''il n''excede pas : %0.5f N.m ',aa3);
           disp('       ');
           c1=input('Reponse : ');
           while (c1>aa3 || isempty(c1) || c1<0)
                disp('       ');
                c1=input('Reponse : ');
           end
        end
        if q==4
           aa4=0.01*Cn*Wn;
           disp('Sachant que Cr = k / (omega) : ');
           disp('       ');
           fprintf('Inserer k sachant qu''il n''excede pas : %0.5f N.m.rad/s ',aa4);
           disp('       ');
           d1=input('Reponse : ');
            while (d1>aa4 || isempty(d1) || d1<0)
                disp('       ');
                d1=input('Reponse : ');
            end
        end
        %% Resistance Rotorique Variable
        if (mot1==2)
            disp('_____________________________________________________________________');
            disp('     ');
            disp('Voulez vous : ');
            disp('     ');
            disp('1. Une resistance rotorique variable avec la temperature ?');
            disp('     ');
            disp('2. Une resistance rotorique constante ?');
            disp('     ');
            x=input('Reponse : ');
            while(isempty(x) || (x~=1 && x~=2))
                disp('     ');
                x=input('Reponse : ');
            end 
            if (x==1)
                variable = 0;
            end
            if (x==2)
                variable = 1;
            end
        else
            variable = 1;
        end 
        
    %% Resultats de la simulation  
    if (mot1==2 || mot1==3)
            if (mot1==2)
                Tp=[0];
                Omg=[0];
                consi=0;
                disp('     ');
                disp('_____________________________________________________________________');
                disp('     ');
                disp('Voulez vous :');
                disp('     ');
                disp('1. Une consigne de vitesse sous forme d''echelon ?');
                disp('     ');
                disp('2. Une consigne de vitesse en escaliers ?');
                disp('     ');
                cons=input('Reponse : ');
                while(isempty(cons) || (cons~=1 && cons~=2))
                    disp('     ');
                    cons=input('Reponse : ');
                end
            end
            
            if (mot1==3)
                disp('     ');
                disp('_____________________________________________________________________');
                disp('     ');
                disp('Choisissez la consigne que vous desirez: ');
                disp('     ');
                disp('Par exemple: consigne = Wn/10. ');
                disp('     ');
                cons=input('Consigne = ');
                siz=size(cons);
                while(isempty(cons) || abs(cons) > Wb || (siz(1)~=1 && siz(2)~=1))
                    disp('     ');
                    disp('La consigne doit etre un seul entier relatif plus petit que (+/-) la vitesse de base !');
                    cons=input('Consigne = ');
                    siz=size(cons);
                end
            end 

            if (cons==1)
                consigne=1;
                disp('     ');
                disp('_____________________________________________________________________');
                disp('     ');
                disp('Choisissez la consigne que vous desirez: ');
                disp('     ');
                disp('Par exemple: consigne = Wn/10. Avec consigne <= Wb ! ');
                disp('     ');
                consi=input('Consigne = ');
                siz=size(consi);
                while(isempty(consi) || abs(consi) > Wb || (siz(1)~=1 && siz(2)~=1))
                    disp('     ');
                    disp('La consigne doit etre un seul entier relatif plus petit que (+/-) la vitesse de base !');
                    consi=input('Consigne = ');
                    siz=size(consi);
                end 
            end
            if (cons==2)
                consigne=0;
                disp('     ');
                disp('_____________________________________________________________________');
                disp('     ');
                disp('Choisissez les vecteur de vitesse et de temps que vous desirez: ');
                disp('     ');
                disp('Par exemple: Tps=[1,3,5,10] ; Omg=[Wn/3, Wn/2, Wn/4, Wn]. ');
                disp('     ');
                disp('Attention les vecteurs Tps et Omg doivent etre de meme longueur. ');
                disp('     ');
                disp('Attention le premier element de Tp doit toujours etre plus grand ou egal a 1. ');
                disp('     ');
                Tp=input('Tps = ');
                s=size(Tp);
                while(isempty(Tp) || s(1)~=1)
                    disp('     ');
                    Tp=input('Tps = ');
                    s=size(Tp);
                end 
                disp('     ');
                Omg=input('Omg = ');
                so=size(Omg);
                while(isempty(Omg) || so(1)~=1)
                    disp('     ');
                    Omg=input('Omg = ');
                    so=size(Omg);
                end 
                while(length(Tp)~=length(Omg))
                    disp('     ');
                    disp('Les longueurs des vecteurs ne sont pas egales !');
                    disp('     ');
                    Tp=input('Tps = ');
                    disp('     ');
                    Omg=input('Omg = ');
                end
            end
            
        if (variable==0 && mot1==2)
            sim('Rr_Variable.slx');
        end
        if (variable ==1 && mot1==2)
            sim('MAS_BF.slx');
        end 
        
        if mot1==3
            fs=5001;
            sim('PWM_inverter.slx');
            
            figure()
            plot(t,Va_MLI,'linewidth',2);
            axis([3.2 3.225 -500 500]),grid;
            xlabel('Temps (s)');
            ylabel('Tension (V)');
            hold on
            plot(t,Va_Ref,'linewidth',2);
            legend('Va MLI','Va reference');
            
        end 
        
        repo=0;
        while(repo~=9)
            disp('         ');
            disp('__________________________________________________________________');
            disp('         ');
            disp('Voulez-vous : ');
            disp('         ');
            disp('1. Visualiser la vitesse mecanique de la MAS ?');
            disp('         ');
            disp('2. Visualiser les courants statoriques Isa, Isb et Isc ?');
            disp('         ');
            disp('3. Visualiser les grandeurs dq (Tensions statoriques, Courants statoriques, flux rotorique) ?');
            disp('         ');
            disp('4. Visualiser la variation de la resistance rotorique en fonction de la Temperature ?');
            disp('         ');
            disp('5. Visualiser la variation du flux en fonction du temps ?');
            disp('         ');
            disp('6. Visualiser la vitesse, le flux et le courant (statorique) sur la meme figure ?');
            disp('         ');
            disp('7. Visualiser les flux rotoriques dans le repere dq ?');
            disp('         ');
            disp('8. Visualier l''angle Ro ?');
            disp('         ');
            disp('9. Aucune courbe');
            disp('         ');
            repo=input('Reponse : ');
            while(isempty(repo) || (repo~=1 && repo~=2 && repo~=3 && repo~=4 && repo~=5 && repo~=6 && repo~=7 && repo~=8 && repo~=9))
                disp('         ');
                repo=input('Reponse : ');
            end
            if (repo==1)
                figure()
                plot(t,omega,'linewidth',2),grid;
                hold on
                plot(t,omega_ref,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Omega (rad/s)');
                title('Variation de la vitesse mecanique de la MAS en fonction du temps');
            end 

            if (repo==2)
                figure()
                subplot(3,1,1)
                plot(t,I_sa,'linewidth',2),grid;
                axis([5 6 -300 300]);
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                title('Variation du courant statorique I_a en fonction du temps');
                subplot(3,1,2)
                plot(t,I_sb,'linewidth',2),grid;
                axis([5 6 -300 300]);
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                title('Variation du courant statorique I_b en fonction du temps');
                subplot(3,1,3)
                plot(t,I_sc,'linewidth',2),grid;
                axis([5 6 -300 300]);
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                title('Variation du courant statorique I_c en fonction du temps');
            end 

            if (repo==3)
                figure()
                subplot(3,1,1)
                plot(t,I_sd,'linewidth',2),grid;
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                hold on
                plot(t,I_sq,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                title('Variation des courants statoriques dq en fonction du temps');
                legend('Isd','Isq');
                subplot(3,1,2)
                plot(t,Vd,'linewidth',2),grid;
                xlabel('Temps (s)');
                ylabel('Tension (V)');
                hold on
                plot(t,Vq,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Tension (V)');
                title('Variation des tensions statoriques dq en fonction du temps');
                legend('Vd','Vq');
                subplot(3,1,3)
                plot(t,phird,'linewidth',2),grid;
                xlabel('Temps (s)');
                ylabel('Flux rotorique suivant l''axe d (Wb)');
                hold on
                plot(t,phirq,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Flux rotorique suivant l''axe q (Wb)');
                title('Variation des flux rotoriques dq en fonction du temps');
                legend('PhiRd','PhiRq');
            end

            if (repo==4)
                figure()
                plot(t,R_variable,'linewidth',2),grid;
                hold on
                plot(t,R_ref,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Resistance (Ohms)');
                title('Variation de la resistance rotorique en fontion de la temperature et du temps');
            end

            if (repo==5)
                figure()
                plot(t,flux,'linewidth',2),grid;
                hold on
                plot(t,flux_ref,'linewidth',2);
                xlabel('Temps (s)');
                ylabel('Flux');
                title('Variation du flux en fonction du temps');
            end
            
            if (repo==6)
                figure()
                subplot(3,1,1)
                plot(t,omega,'linewidth',2),grid;
                hold on
                plot(t,omega_ref,'linewidth',2);
                axis([0.9 7 0 3100]);
                xlabel('Temps (s)');
                ylabel('Omega (tr/min)');
                title('Variation de la vitesse mecanique de la MAS en fonction du temps');
                legend('Omega','Omega Ref');
                subplot(3,1,2);
                plot(t,flux,'linewidth',2),grid;
                hold on
                plot(t,flux_ref,'linewidth',2);
                axis([0 5 0 62]);
                xlabel('Temps (s)');
                ylabel('Imr (A)');
                title('Variation du courant Imr en fonction du temps');
                legend('Imr','Imr^*');
                subplot(3,1,3)
                plot(t,Couple_ref,'linewidth',2),grid;
                hold on
                plot(t,Couple,'linewidth',2);
                axis([0 7 0 350]);
                xlabel('Temps (s)');
                ylabel('Couple (N.m)');
                title('Variation du couple en fonction du temps');
                legend('Couple','Couple Ref');
            end
            
            if(repo==7)
                figure()
                subplot(2,1,1)
                plot(t,phird,'linewidth',2),grid;
                xlabel('Temps (s)');
                ylabel('Flux rotorique suivant l''axe d (Wb)');
                title('Variation du flux rototrique suivant l''axe d en fonction du temps');
                subplot(2,1,2)
                plot(t,phirq,'linewidth',2),grid;
                axis([0 T_s -1 1]);
                xlabel('Temps (s)');
                ylabel('Flux rotorique suivant l''axe q (Wb)');
                title('Variation du flux rotorique suivant l''axe q en fonction du temps');
            end
            
            if (repo==8)
                figure()
                plot(t,Rho,'linewidth',2),grid;
                axis([0.9 1.5 0 7]);
                xlabel('Temps (s)');
                ylabel('Angle (rad)');
                title('Variation de l''angle Ro en fonction du temps.');
            end 

        end
    end 
    %% MAS en BO
    
    if (mot1==1)
        sim('MAS_180_BO.slx');
        figure()
        subplot(2,1,1)
        plot(tps,Omega,'linewidth',2),grid;
        xlabel('Temps (s)');
        ylabel('Vitesse (tr/min)');
        title('Variation de la vitesse en fonction du temps');
        subplot(2,1,2)
        plot(tps,cem,'linewidth',2),grid;
        xlabel('Temps (s)');
        ylabel('Couple (N.m)');
        title('Variation du couple electromagnetique en fonction du temps.');
    end 
    
    %% Boucle While 
    
    disp('       ');
    disp('_____________________________________________________________________');
    disp('       ');
    disp('Voulez vous : ');
    disp('        ');
    disp('1. Simuler la machine en BO ? ');
    disp('        ');
    disp('2. Simuler la machine en BF avec commande vectorielle ? ');
    disp('        ');
    disp('3. Simuler la machine en BF avec onduleur ? ');
    disp('        ');
    disp('4. Quitter le programme ? ');
    disp('        ');
    mot1=input('Reponse : ');
    while(isempty(mot1)|| (mot1~=1 && mot1~=2 && mot1~=3 && mot1~=4))
        disp('        ');
        mot1=input('Reponse : ');
    end 
    disp('        ');
    disp('_____________________________________________________________________');
    if (mot1==4)
        cte=2;
    end 
end 