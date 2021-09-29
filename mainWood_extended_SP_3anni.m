% Main Wood - Stima parametri per le funzioni di mortalita'
clear all
close all
clc

%% Carico i dati: temperature e rilievi 
% ammettiamo ora che i rilevi avvengano in giorni differenti
% uova, larve, pupe sono rilevati contemporaneamente
loadDati2008
loadDati2009
loadDati2011

loadDati2010
loadDati2012

%% Parametri da stimare, base e Guess iniziale
% base: Bspline
% coefficienti: 7 x n_stadi

% N=7;
% Knots=[0 0 0 0 10 20 30 40 40 40 40];
N=5;
Knots=[0 0 0 0 20 40 40 40 40];
for i=1:N
 BS(i)= bspline(Knots(i:i+4));
end



p=[ 0.5240    0.5240    0.5240    0.5240;
%    0.2268    0.2268    0.2268    0.2268;
   -0.0105   -0.0105   -0.0105   -0.0105;
    0.1    0.1    0.1    0.1;
   -0.0119   -0.0119   -0.0119   -0.0119;
%    0.2183    0.2183    0.2183    0.2183;
    1.2696    1.2696    1.2696    1.2696];

% p=[ 2.8    2.8   2.8    2.8;
%     0.2268    0.2268    0.2268    0.2268;
%    -0.0105   -0.0105   -0.0105   -0.0105;
%     0.01    0.01    0.01    0.01;
%    -0.0119   -0.0119   -0.0119   -0.0119;
%     0.2183    0.2183    0.2183    0.2183;
%     2.8    2.8    2.8    2.8];

% p=[ 1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
%    0.500000000000000   0.500000000000000   0.500000000000000   0.500000000000000
%    0.020000000000000   0.020000000000000   0.020000000000000   0.020000000000000
%    0.010000000000000   0.010000000000000   0.010000000000000   0.010000000000000
%    0.080000000000000   0.080000000000000   0.080000000000000   0.080000000000000
%    0.218300000000000   0.218300000000000   0.218300000000000   0.218300000000000
%    0.400000000000000   0.400000000000000   0.400000000000000   0.400000000000000];

% p=[ 0.1000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
%    0.1500000000000000   0.500000000000000   0.500000000000000   0.500000000000000
%    0.020000000000000   0.020000000000000   0.020000000000000   0.020000000000000
%    0.010000000000000   0.010000000000000   0.010000000000000   0.010000000000000
%    0.080000000000000   0.080000000000000   0.080000000000000   0.080000000000000
%    0.218300000000000   0.218300000000000   0.218300000000000   0.218300000000000
%    0.400000000000000   0.400000000000000   0.400000000000000   0.400000000000000];

% p=[49.2225112259184 5.74040736009983 31.6905168978969 3.27202419302683;
%     13.5465806810932 0.654439036732514 -0.651033615329903 23.7449046166775;
%     2.53811183710093 3.26226067611894 6.88960926349849 -0.333981093097549;
%     -1.13321332215428 2.80493328453480 1.19784063173859 -2.66924949866476;
%     67.5578875461499 4.46751176346099 2.49214721251958 119.538334356289;
%     -2.95788230771607 2.74527094580265 -0.456430622137304 -12.1828659074163;
%     113.190654803651 7.36193139204698 10.0796419375194 91.0483503014060];

% p=p*0.5;
% p=zeros(N,n_stadi);
np_s=numel(p(:,1));%numero di parametri per stadio
np=np_s*ns;%numero totale parametri

% Guess iniziale
% p=0.1*ones(N,n_stadi);




% Y_2009=Y_2009+randn(numel(dreg_2009),4);
% Y_2009=Y_2009.*(Y_2009>0);

% figure()
%     for j=1:ns
%         subplot(2,2,j)
%         plot(0:0.1:40,fm_Wood(0:0.1:40,p(:,j),BS))
%         switch j
%         case 1
%             title('Eggs')
%         case 2
%             title('Larvae')
%         case 3
%             title('Pupae')
%          case 4
%             title('Adults')
%         end
%     end 

% fid = fopen('p.txt','w');
% fprintf(fid,'%12.8f %12.8f %12.8f %12.8f\n',p);
% fclose(fid);
%% Ciclo
nit=1;
% Funzionale senza pesi
%funzionale=@(pp,y,J) (y-J*pp(:))'*(y-J*pp(:));

% Funzionale con peso sulle larve
W_2008=0.1*diag([ones(1,d1_2008) 1*ones(1,d2_2008) ones(1,d3_2008) ones(1,d4_2008)]);
W_2009=0.1*diag([ones(1,d1_2009) 1*ones(1,d2_2009) ones(1,d3_2009) ones(1,d4_2009)]);
W_2011=0.1*diag([ones(1,d1_2011) 1*ones(1,d2_2011) ones(1,d3_2011) ones(1,d4_2011)]);
funzionale=@(pp,y,J,W) (y-J*pp(:))'*W*(y-J*pp(:));


%options=optimset('MaxFunEvals',10000,'MaxIter',10000);
delta=1e-6; % perturbazione

% Constraints
A=zeros(np,np);
A(1,1)=-1;
A(np_s+1,np_s+1)=-1;
A(2*np_s+1,2*np_s+1)=-1;
A(3*np_s+1,3*np_s+1)=-1;
b=zeros(np,1);

%% Integrazione
% calcolo la mortalita' per ogni stadio ad ogni temperatura per il p
% considerato
tic;
MU_2008=Giro2008(p,T_2008,BS,ns,d_imm_2008,d_ad_2008);
MU_2009=Giro2009(p,T_2009,BS,ns,d_imm_2009,d_ad_2009);
MU_2011=Giro2011(p,T_2011,BS,ns,d_imm_2011,d_ad_2011); 


%%
cond=ExitCondEval(Y_2008,MU_2008,d1_2008,d2_2008)...
    +ExitCondEval(Y_2009,MU_2009,d1_2009,d2_2009)...
    +ExitCondEval(Y_2011,MU_2011,d1_2011,d2_2011);
COND(nit)=cond;
P(:,:,nit)=p;
disp(cond)
%% 
opts = optimoptions('fmincon', 'MaxFunctionEvaluations', 2000000,'MaxIterations',5000);
CC=1;

while (cond>1 && CC~=0)
    
    % estimate of J
   
    for k=1:np
            disp(k)
            tic;
            pert=zeros(np_s,ns);
            pert(k)=delta;
            %MUp1=Kolmogorov_Wood(n_stadi,ngen,nt,y0_2009,T,Tmed,td,p+pert,BS,dreg_2009);
            ppp=p+pert;
            MUp1_2008=Giro2008(ppp,T_2008,BS,ns,d_imm_2008,d_ad_2008);
            MUp1_2009=Giro2009(ppp,T_2009,BS,ns,d_imm_2009,d_ad_2009);
            MUp1_2011=Giro2011(ppp,T_2011,BS,ns,d_imm_2011,d_ad_2011);
%             MUp1_2012=Giro2012(ppp,T_2012,BS,ns,d_imm_2012,d_ad_2012);
            
            ppp=p-pert;
            MUm1_2008=Giro2008(ppp,T_2008,BS,ns,d_imm_2008,d_ad_2008);
            MUm1_2009=Giro2009(ppp,T_2009,BS,ns,d_imm_2009,d_ad_2009);
            MUm1_2011=Giro2011(ppp,T_2011,BS,ns,d_imm_2011,d_ad_2011);
%             MUm1_2012=Giro2012(ppp,T_2012,BS,ns,d_imm_2012,d_ad_2012);
            
            Jk_2008(:,k)=(MUp1_2008(:)-MUm1_2008(:))/(2*delta);
            Jk_2009(:,k)=(MUp1_2009(:)-MUm1_2009(:))/(2*delta);
            Jk_2011(:,k)=(MUp1_2011(:)-MUm1_2011(:))/(2*delta);
%             Jk_2012(:,k)=(MUp1_2012(:)-MUm1_2012(:))/(2*delta);
            
            toc;
    end
    %YY=Y_2009';
    %YY=Y_2009;
    %[~,Y_2009]=MatriceRilievi(dreg_imm_2009,dreg_ad_2009,rilievi_imm,rilievi_adulti_2009,MU);
    %yk=Y_2009(:)-MU(:)+Jk*p(:);
    yk_2008=Y_2008-MU_2008+Jk_2008*p(:);
    yk_2009=Y_2009-MU_2009+Jk_2009*p(:);
    yk_2011=Y_2011-MU_2011+Jk_2011*p(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % minimize the objective
    %p=fminsearch(@(pp) funzionale(pp,yk,Jk),p);
    lb=zeros(N,4);
%      lb=-1*ones(7,4);
    ub=2.*ones(N,4);
    %p=fmincon(@(pp) funzionale(pp,yk,Jk),p,A,b);
    p=fmincon(@(pp) funzionaleCONincertezza(pp,yk_2008,incertezza_2008,Jk_2008,W_2008)...
        +funzionaleCONincertezza(pp,yk_2009,incertezza_2009,Jk_2009,W_2009)...        
        +funzionaleCONincertezza(pp,yk_2011,incertezza_2011,Jk_2011,W_2011),... %         +funzionale(pp,yk_2012,Jk_2012,W_2012),...
        p,[],[],[],[],lb,ub,@(pp) positivity(pp,BS),opts);     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     p=fmincon(@(pp) funzionale(pp,yk_2008,Jk_2008,W_2008)...
%         +funzionale(pp,yk_2009,Jk_2009,W_2009)...        
%         +funzionale(pp,yk_2011,Jk_2011,W_2011),... %         +funzionale(pp,yk_2012,Jk_2012,W_2012),...
%         p,[],[],[],[],lb,ub,@(pp) positivity(pp,BS),opts);     %%
    disp(p)
    % numerically solve the model equations
    MU_2008=Giro2008(p,T_2008,BS,ns,d_imm_2008,d_ad_2008);
    MU_2009=Giro2009(p,T_2009,BS,ns,d_imm_2009,d_ad_2009);
    MU_2011=Giro2011(p,T_2011,BS,ns,d_imm_2011,d_ad_2011);
    
    cond=ExitCondEval(Y_2008,MU_2008,d1_2008,d2_2008)...
        +ExitCondEval(Y_2009,MU_2009,d1_2009,d2_2009)...
         +ExitCondEval(Y_2011,MU_2011,d1_2011,d2_2011);
%         +ExitCondEval(Y_2012,MU_2012,d1_2012,d2_2012);
    disp(cond)
    COND(nit+1)=cond;
    CC=COND(nit+1)-COND(nit)
    P(:,:,nit+1)=p;
    nit=nit+1;
    save('Sim_2008_2009_2011_pesiuguali_N5_nuovofunzionale')
end
%%

[~,Dyn_2008]=Giro2008(p,T_2008,BS,ns,d_imm_2008,d_ad_2008);
time_2008=linspace(dreg_iniz_2008/24,dreg_fin_2008/24,length(Dyn_2008(:,1)));
[~,Dyn_2009]=Giro2009(p,T_2009,BS,ns,d_imm_2009,d_ad_2009);
time_2009=linspace(dreg_iniz_2009/24,dreg_fin_2009/24,length(Dyn_2009(:,1)));
[~,Dyn_2011]=Giro2011(p,T_2011,BS,ns,d_imm_2011,d_ad_2011);
time_2011=linspace(dreg_iniz_2011/24,dreg_fin_2011/24,length(Dyn_2011(:,1)));
%%
[~,Dyn_2010]=Giro2010(p,T_2010,BS,ns,d_imm_2010,d_ad_2011);
time_2010=linspace(dreg_iniz_2010/24,dreg_fin_2010/24,length(Dyn_2010(:,1)));
[~,Dyn_2012]=Giro2012(p,T_2012,BS,ns,d_imm_2012,d_ad_2012);
time_2012=linspace(dreg_iniz_2012/24,dreg_fin_2012/24,length(Dyn_2012(:,1)));


%% Plot dinamica
%     close all
         le2008=length(rilievi_uova_2008);
         ll2008=length(rilievi_larve_2008);
         lp2008=length(rilievi_pupe_2008);
         la2008=length(rilievi_adulti_2008);
         
         le2009=length(rilievi_uova_2009);
         ll2009=length(rilievi_larve_2009);
         lp2009=length(rilievi_pupe_2009);
         la2009=length(rilievi_adulti_2009);
         
         le2011=length(rilievi_uova_2011);
         ll2011=length(rilievi_larve_2011);
         lp2011=length(rilievi_pupe_2011);
         la2011=length(rilievi_adulti_2011);
         
         le2010=length(rilievi_uova_2010);
         ll2010=length(rilievi_larve_2010);
         lp2010=length(rilievi_pupe_2010);
         la2010=length(rilievi_adulti_2010);
         
         le2012=length(rilievi_uova_2012);
         ll2012=length(rilievi_larve_2012);
         lp2012=length(rilievi_pupe_2012);
         la2012=length(rilievi_adulti_2012);
         
        figure(2008)
    for j=1:ns
        subplot(2,2,j)
        hold on
        plot(time_2008,Dyn_2008(:,j),'-b')
        %plot(dreg_2008/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            errorbar(dreg_imm_2008/24,rilievi_uova_2008,0*incertezza_2008(1:le2008),incertezza_2008(1:le2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_uova_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            errorbar(dreg_imm_2008/24,rilievi_larve_2008,0*incertezza_2008(le2008+1:ll2008+le2008),incertezza_2008(le2008+1:ll2008+le2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_larve_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            errorbar(dreg_imm_2008/24,rilievi_pupe_2008,0*incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_pupe_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            errorbar(dreg_ad_2008/24,rilievi_adulti_2008,0*incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2008/24,rilievi_adulti_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    title('2008')
    
    figure(2009)
    for j=1:ns
        subplot(2,2,j)
        hold on
        plot(time_2009,Dyn_2009(:,j),'-b')
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            errorbar(dreg_imm_2009/24,rilievi_uova_2009,0*incertezza_2009(1:le2009),incertezza_2009(1:le2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_uova_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            errorbar(dreg_imm_2009/24,rilievi_larve_2009,0*incertezza_2009(le2009+1:ll2009+le2009),incertezza_2009(le2009+1:ll2009+le2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_larve_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            errorbar(dreg_imm_2009/24,rilievi_pupe_2009,0*incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_pupe_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            errorbar(dreg_ad_2009/24,rilievi_adulti_2009,0*incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2009/24,rilievi_adulti_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    title('2009')
    
    %2011
    figure(2011)
    for j=1:ns
        subplot(2,2,j)
        hold on
        plot(time_2011,Dyn_2011(:,j),'-b')
        %plot(dreg/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            errorbar(dreg_imm_2011/24,rilievi_uova_2011,0*incertezza_2011(1:le2011),incertezza_2011(1:le2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_uova_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            errorbar(dreg_imm_2011/24,rilievi_larve_2011,0*incertezza_2011(le2011+1:ll2011+le2011),incertezza_2011(le2011+1:ll2011+le2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_larve_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            errorbar(dreg_imm_2011/24,rilievi_pupe_2011,0*incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_pupe_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            errorbar(dreg_ad_2011/24,rilievi_adulti_2011,0*incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2011/24,rilievi_adulti_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    title('2011')
    
    
    %2010
    figure(2010)
    for j=1:ns
        subplot(2,2,j)
        hold on
        plot(time_2010,Dyn_2010(:,j),'-k')
        %plot(dreg/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            errorbar(dreg_imm_2010/24,rilievi_uova_2010,0*incertezza_2010(1:le2010),incertezza_2010(1:le2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_uova_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            errorbar(dreg_imm_2010/24,rilievi_larve_2010,0*incertezza_2010(le2010+1:ll2010+le2010),incertezza_2010(le2010+1:ll2010+le2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_larve_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            errorbar(dreg_imm_2010/24,rilievi_pupe_2010,0*incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_pupe_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            errorbar(dreg_ad_2010/24,rilievi_adulti_2010,0*incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2010/24,rilievi_adulti_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    title('2010')
    
    %2011
    figure(2012)
    for j=1:ns
        subplot(2,2,j)
        hold on
        plot(time_2012,Dyn_2012(:,j),'-k')
        %plot(dreg/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            errorbar(dreg_imm_2012/24,rilievi_uova_2012,0*incertezza_2012(1:le2012),incertezza_2012(1:le2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_uova_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            errorbar(dreg_imm_2012/24,rilievi_larve_2012,0*incertezza_2012(le2012+1:ll2012+le2012),incertezza_2012(le2012+1:ll2012+le2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_larve_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            errorbar(dreg_imm_2012/24,rilievi_pupe_2012,0*incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_pupe_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            errorbar(dreg_ad_2012/24,rilievi_adulti_2012,0*incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2012/24,rilievi_adulti_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    title('2012')
    
    
%%


% Plot dinamica
%     close all
    figure(2008)
    for j=1:ns-1
        subplot(3,3,j)
        hold on
        plot(time_2008,Dyn_2008(:,j),'-b','LineWidth',1)
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('time(days)')
        ax=gca;
        ax.FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2008/24,rilievi_uova_2008,0*incertezza_2008(1:le2008),incertezza_2008(1:le2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_uova_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2008/24,rilievi_larve_2008,0*incertezza_2008(le2008+1:ll2008+le2008),incertezza_2008(le2008+1:ll2008+le2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_larve_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2008')
        case 3
            errorbar(dreg_imm_2008/24,rilievi_pupe_2008,0*incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2008/24,rilievi_pupe_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2008/24,rilievi_adulti_2008,0*incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2008/24,rilievi_adulti_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end


% Plot dinamica
%     close all
%     figure(2009)
    for j=1:ns-1
        subplot(3,3,j+3)
        hold on
        plot(time_2009,Dyn_2009(:,j),'-b','LineWidth',1)
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('days')
        ax=gca;
        ax.FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2009/24,rilievi_uova_2009,0*incertezza_2009(1:le2009),incertezza_2009(1:le2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_uova_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2009/24,rilievi_larve_2009,0*incertezza_2009(le2009+1:ll2009+le2009),incertezza_2009(le2009+1:ll2009+le2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_larve_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2009')
        case 3
            errorbar(dreg_imm_2009/24,rilievi_pupe_2009,0*incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2009/24,rilievi_pupe_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2009/24,rilievi_adulti_2009,0*incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2009/24,rilievi_adulti_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end


% Plot dinamica
%     close all
%     figure(2011)
    for j=1:ns-1
        subplot(3,3,j+6)
        hold on
        plot(time_2011,Dyn_2011(:,j),'-b','LineWidth',1)

        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('time (days)')
        ax=gca;
        ax.FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2011/24,rilievi_uova_2011,0*incertezza_2011(1:le2011),incertezza_2011(1:le2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_uova_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2011/24,rilievi_larve_2011,0*incertezza_2011(le2011+1:ll2011+le2011),incertezza_2011(le2011+1:ll2011+le2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_larve_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2011')
        case 3
            errorbar(dreg_imm_2011/24,rilievi_pupe_2011,0*incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2011/24,rilievi_pupe_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2011/24,rilievi_adulti_2011,0*incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2011/24,rilievi_adulti_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
    
    
    
    % Plot dinamica
%     close all
    figure(2010)
    for j=1:ns-1
        subplot(3,3,j)
        hold on
        plot(time_2010,Dyn_2010(:,j),'-b','LineWidth',1)
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('time(days)')
        ax=gca;
        ax.FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2010/24,rilievi_uova_2010,0*incertezza_2010(1:le2010),incertezza_2010(1:le2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_uova_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
             ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2010/24,rilievi_larve_2010,0*incertezza_2010(le2010+1:ll2010+le2010),incertezza_2010(le2010+1:ll2010+le2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_larve_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2010')
        case 3
            errorbar(dreg_imm_2010/24,rilievi_pupe_2010,0*incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2010/24,rilievi_pupe_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2010/24,rilievi_adulti_2010,0*incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2010/24,rilievi_adulti_2010,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end


% Plot dinamica
%     close all
%     figure(2009)
    for j=1:ns-1
        subplot(3,3,j+3)
        hold on
        plot(time_2012,Dyn_2012(:,j),'-b','LineWidth',1)
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('days')
        ax=gca;
        ax.FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2012/24,rilievi_uova_2012,0*incertezza_2012(1:le2012),incertezza_2012(1:le2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_uova_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2012/24,rilievi_larve_2012,0*incertezza_2012(le2012+1:ll2012+le2012),incertezza_2012(le2012+1:ll2012+le2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_larve_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2012')
        case 3
            errorbar(dreg_imm_2012/24,rilievi_pupe_2012,0*incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_imm_2012/24,rilievi_pupe_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2012/24,rilievi_adulti_2012,0*incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),'o-r','MarkerSize',2,'MarkerFaceColor','g')
%             plot(dreg_ad_2012/24,rilievi_adulti_2012,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
    
 

  
    
%% Plot della mortalita stimata
    figure(1)
    for j=1:ns
        subplot(2,2,j)
        plot(0:0.1:40,fm_Wood(0:0.1:40,p(:,j),BS),'b','LineWidth',1)

%         hold on
%         plot(0:0.1:40,fm(j,0:0.1:40),'k') % per confronti
        switch j
        case 1
            title('Eggs')
        case 2
            title('Larvae')
        case 3
            title('Pupae')
         case 4
            title('Adults')
        end
    end   
    
    figure()
    plot(0:0.1:40,fm_Wood(0:0.1:40,p(:,1),BS),'b',0:0.1:40,fm_Wood(0:0.1:40,p(:,2),BS),':r',0:0.1:40,fm_Wood(0:0.1:40,p(:,3),BS),'--g',0:0.1:40,fm_Wood(0:0.1:40,p(:,4),BS),'-.c','LineWidth',1)
    xlabel('T (^oC)')
    ylabel('mortality rate (d^{-1})')
    legend('eggs','larvae','pupae','adults')
    axis([0 40 0 2])

    
 %% Plot delle mortalita stimate
    figure()
    for kk=1:numel(COND)
    for j=1:ns
        subplot(2,2,j)
        plot(0:0.1:40,fm_Wood(0:0.1:40,P(:,j,kk),BS))
         hold on
%         plot(0:0.1:40,fm(j,0:0.1:40),'k') % per confronti
        switch j
        case 1
            title('Eggs')
        case 2
            title('Larvae')
        case 3
            title('Pupae')
         case 4
            title('Adults')
        end
    end     
    end
 
    
%% Seleziono il miglior valore
[a,b]=min(COND);
b=100  %%risultati migliori con 6 (9 per Cinzia)
p=P(:,:,b);
% Ora cancello Jk per evitare problemi con le dimensioni
clear Jk

%% Jacobiano come Marsili-Libelli
for k=1:np
            disp(k)
            tic;
            pert=zeros(np_s,ns);
            pert(k)=delta;
            %MUp1=Kolmogorov_Wood(n_stadi,ngen,nt,y0_2009,T,Tmed,td,p+pert,BS,dreg_2009);
            ppp=p+pert;
            MUp1_2008=Giro2008(ppp,T_2008,BS,ns,d_imm_2008,d_ad_2008);
            MUp1_2009=Giro2009(ppp,T_2009,BS,ns,d_imm_2009,d_ad_2009);
            MUp1_2011=Giro2011(ppp,T_2011,BS,ns,d_imm_2011,d_ad_2011);
            
            %MUm1=Kolmogorov_Wood(n_stadi,ngen,nt,y0_2009,T,Tmed,td,p-pert,BS,dreg_2009);
            ppp=p-pert;
            MUm1_2008=Giro2008(ppp,T_2008,BS,ns,d_imm_2008,d_ad_2008);
            MUm1_2009=Giro2009(ppp,T_2009,BS,ns,d_imm_2009,d_ad_2009);
            MUm1_2011=Giro2011(ppp,T_2011,BS,ns,d_imm_2011,d_ad_2011);
            
%             MUp1-MUm1
            
            % Usando la matrice del Wood
            Jk_2008(:,k)=(MUp1_2008(:)-MUm1_2008(:))/(2*delta);
            Jk_2009(:,k)=(MUp1_2009(:)-MUm1_2009(:))/(2*delta);
            Jk_2011(:,k)=(MUp1_2011(:)-MUm1_2011(:))/(2*delta);
            % Usando quella di Marsili-Libelli
            %Jk(:,k)=(sum(MUp1-MUm1)/(2*delta));
            toc;
end
Jk=[Jk_2008;Jk_2009;Jk_2011];

%% Tolgo le colonne nulle
JJk=[];
Inulli=[];
Iok=[];
% Tolgo le colonne tutte nulle
for i=1:np
    if  ne(sum(ne(Jk(:,i), zeros(size(Jk(:,i))))),0)  
        JJk=[JJk Jk(:,i)];
        Iok=[Iok i];
    else
        Inulli=[Inulli i];
    end
end
        
nnp=np-numel(Inulli);

%% Matrice Covarianze
MU=[MU_2008;MU_2009;MU_2011];
W=blkdiag(W_2008,W_2009,W_2011);
YY=[Y_2008;Y_2009;Y_2011];
Ep=((YY(:)-MU(:))')*W*(YY(:)-MU(:));
nreg=numel(dreg_2008)+numel(dreg_2009)+numel(dreg_2011);
% FIM (Jk matrice 4x28)
%CJp=Ep/(nreg-nnp)*inv(((JJk')*((diag([0.1 1 0.1 0.1]))^2\eye(4))*JJk));
% Wood (Jk matrice (dregx4)x28)
CJp=(Ep/(nreg-nnp))*((((JJk')*(W\eye(nreg))*JJk))\eye(nnp));

% for i=1:length(CJp)
%     if abs(CJp(i,i))>10
%         CJp(i,i)=CJp(i,i)/100;
%     end
% end

alfa=0.01;
tstud=tinv(1-alfa/2,nreg-nnp);
% %confInt=tstud*sqrt(diag(CJp));
tstud=1;
VarP=tstud.*sqrt(diag(CJp));


% for kk=1:27
%     if VarP(kk)>0.5
%         VarP(kk)=0.5;
%     end
% end

% VarP(1)=1.5;
% VarP(7)=1.5;
% VarP(8)=1.5;
% VarP(27)=1.5;

%% Prove per Bande
% simulated dynamics within the confidence interval

nprove=500;
%p
%var

pmean=P(:,:,b)
pvar=zeros(1,np);
pvar(Iok)=VarP;

%Pint=pmean(:)+pvar/1.96.*randn(28,nprove);

pmean=reshape(pmean,4*N,1);
CCJpp=(CJp+CJp')/2; %trucchetto per evitare il problema del det basso

TT=0:0.01:40;

for i=1:nprove
    disp(i)
    %
    %%%% estrazione di ogni parametro da una normal
%     Pint(:,i)=normrnd(pmean(:),pvar');
    %%%%%
    
%     M(1,1)=-1;
%     lbb=zeros(length(Iok),1)-pmean(Iok);
%     ubb=zeros(length(Iok),1)+10*pmean(Iok);
%     lbb=zeros(length(Iok),1)-(ones(length(Iok),1)-mvncdf(-inv(CCJpp)*pmean(Iok)))*pmean(Iok);
%     ubb=ones(length(Iok),1)*Inf;
%     while min(min(M))<0
    
%    
%      CCJpp=CJp-diag(diag(CJp))+diag(ones(26,1))
%      Pint(Iok,i)=mvrandn(zeros(26,1)-pmean(Iok),ones(26,1)*Inf,CJp,1)+pmean(Iok);

%%% estrazione da una normale multivariata
%     CCJpp=nearestSPD(CCJpp); %prende la matrice semidef. pos. piu'
%     vicina a CCJpp
    ctest1=-2;
    ctest2=3;
    while ctest1<-0.06 || ctest2>2
       Pint(Iok,i)=mvnrnd(pmean(Iok),CCJpp);
       ctest1=min(Pint(Iok,i));
       ctest2=max(Pint(Iok,i));
    end
    Pint(Inulli,i)=pmean(Inulli);
%%%%%
    
    p=reshape(Pint(:,i),N,4)
    
%     for t=1:numel(TT)
%     for is=1:ns
%         M(t,is)=fm_Wood(TT(t),p(:,is),BS);
%     end
%     end
%     
%     end
%     
%     p=reshape(Pint(:,i),N,4);
    
    for t=1:numel(T_2008)
    for is=1:ns
        M(t,is)=fm_Wood(T_2008(t),p(:,is),BS);
        if M(t,is)<0
            M(t,is)=0;
        end
    end
    end
 
% stampo su file le mortalita
fid = fopen('Meggs_2008.txt','w');
fprintf(fid,'%12.8f\n',M(:,1));
fclose(fid);

fid = fopen('Mlarvae_2008.txt','w');
fprintf(fid,'%12.8f\n',M(:,2));
fclose(fid);

fid = fopen('Mpupae_2008.txt','w');
fprintf(fid,'%12.8f\n',M(:,3));
fclose(fid);

fid = fopen('Mad_2008.txt','w');
fprintf(fid,'%12.8f\n',M(:,4));
fclose(fid);
   
    for t=1:numel(T_2009)
    for is=1:ns
        M(t,is)=fm_Wood(T_2009(t),p(:,is),BS);
        if M(t,is)<0
            M(t,is)=0;
        end
    end
    end
% stampo su file le mortalita
fid = fopen('Meggs_2009.txt','w');
fprintf(fid,'%12.8f\n',M(:,1));
fclose(fid);

fid = fopen('Mlarvae_2009.txt','w');
fprintf(fid,'%12.8f\n',M(:,2));
fclose(fid);

fid = fopen('Mpupae_2009.txt','w');
fprintf(fid,'%12.8f\n',M(:,3));
fclose(fid);

fid = fopen('Mad_2009.txt','w');
fprintf(fid,'%12.8f\n',M(:,4));
fclose(fid);

    for t=1:numel(T_2011)
    for is=1:ns
        M(t,is)=fm_Wood(T_2011(t),p(:,is),BS);
        if M(t,is)<0
            M(t,is)=0;
        end
    end
    end
% stampo su file le mortalita
fid = fopen('Meggs_2011.txt','w');
fprintf(fid,'%12.8f\n',M(:,1));
fclose(fid);

fid = fopen('Mlarvae_2011.txt','w');
fprintf(fid,'%12.8f\n',M(:,2));
fclose(fid);

fid = fopen('Mpupae_2011.txt','w');
fprintf(fid,'%12.8f\n',M(:,3));
fclose(fid);

fid = fopen('Mad_2011.txt','w');
fprintf(fid,'%12.8f\n',M(:,4));
fclose(fid);


    for t=1:numel(T_2010)
    for is=1:ns
        M(t,is)=fm_Wood(T_2010(t),p(:,is),BS);
        if M(t,is)<0
            M(t,is)=0;
        end
    end
    end
% stampo su file le mortalita
fid = fopen('Meggs_2010.txt','w');
fprintf(fid,'%12.8f\n',M(:,1));
fclose(fid);

fid = fopen('Mlarvae_2010.txt','w');
fprintf(fid,'%12.8f\n',M(:,2));
fclose(fid);

fid = fopen('Mpupae_2010.txt','w');
fprintf(fid,'%12.8f\n',M(:,3));
fclose(fid);

fid = fopen('Mad_2010.txt','w');
fprintf(fid,'%12.8f\n',M(:,4));
fclose(fid);

    for t=1:numel(T_2012)
    for is=1:ns
        M(t,is)=fm_Wood(T_2012(t),p(:,is),BS);
        if M(t,is)<0
            M(t,is)=0;
        end
    end
    end
% stampo su file le mortalita
fid = fopen('Meggs_2012.txt','w');
fprintf(fid,'%12.8f\n',M(:,1));
fclose(fid);

fid = fopen('Mlarvae_2012.txt','w');
fprintf(fid,'%12.8f\n',M(:,2));
fclose(fid);

fid = fopen('Mpupae_2012.txt','w');
fprintf(fid,'%12.8f\n',M(:,3));
fclose(fid);

fid = fopen('Mad_2012.txt','w');
fprintf(fid,'%12.8f\n',M(:,4));
fclose(fid);


tic;  
%MU=Kolmogorov_Wood(n_stadi,ngen,nt,y0_2009,T,Tmed,td,p,BS,dreg_2009);
system('./DemograficoLobesia2008.out');
system('./DemograficoLobesia2009.out');
system('./DemograficoLobesia2011.out');
system('./DemograficoLobesia2010.out');
system('./DemograficoLobesia2012.out');
toc;
load res_fortran_2008.txt
load res_fortran_2009.txt
load res_fortran_2011.txt
load res_fortran_2010.txt
load res_fortran_2012.txt
Dyn08(:,:,i)=res_fortran_2008;
Dyn09(:,:,i)=res_fortran_2009;
Dyn11(:,:,i)=res_fortran_2011;
Dyn10(:,:,i)=res_fortran_2010;
Dyn12(:,:,i)=res_fortran_2012;
end
%%
%dreg_iniz_2009=dreg_2009(1);
% dreg_fin_2008=dreg_2008(end);
time08=linspace(dreg_iniz_2008/24,dreg_fin_2008/24,length(res_fortran_2008(:,1)));
% dreg_fin_2009=dreg_2009(end);
time09=linspace(dreg_iniz_2009/24,dreg_fin_2009/24,length(res_fortran_2009(:,1)));
% dreg_fin_2011=dreg_2011(end);
time11=linspace(dreg_iniz_2011/24,dreg_fin_2011/24,length(res_fortran_2011(:,1)));
% dreg_fin_2010=dreg_2010(end);
time10=linspace(dreg_iniz_2010/24,dreg_fin_2010/24,length(res_fortran_2010(:,1)));
% dreg_fin_2012=dreg_2012(end);
time12=linspace(dreg_iniz_2012/24,dreg_fin_2012/24,length(res_fortran_2012(:,1)));


%% Plot dinamica
    close all
    figure(2008)
    
        for j=1:ns
        subplot(2,2,j)
        hold on
        for i=1:nprove
        plot(time08,Dyn08(:,j,i),'-k')
        end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            plot(dreg_imm_2008/24,rilievi_uova_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            plot(dreg_imm_2008/24,rilievi_larve_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            plot(dreg_imm_2008/24,rilievi_pupe_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            plot(dreg_ad_2008/24,rilievi_adulti_2008,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
        end
    
        figure(2009)
    for j=1:ns
        subplot(2,2,j)
        hold on
        for i=1:nprove
        plot(time09,Dyn09(:,j,i),'-k')
        end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            plot(dreg_imm_2009/24,rilievi_uova_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            plot(dreg_imm_2009/24,rilievi_larve_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            plot(dreg_imm_2009/24,rilievi_pupe_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            plot(dreg_ad_2009/24,rilievi_adulti_2009,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
    end
    
    figure(2011)
        for j=1:ns
        subplot(2,2,j)
        hold on
        for i=1:nprove
        plot(time11,Dyn11(:,j,i),'-k')
        end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('days')
        switch j
        case 1
            plot(dreg_imm_2011/24,rilievi_uova_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Eggs')
        case 2
            plot(dreg_imm_2011/24,rilievi_larve_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Larvae')
        case 3
            plot(dreg_imm_2011/24,rilievi_pupe_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Pupae')
        case 4
            plot(dreg_ad_2011/24,rilievi_adulti_2011,'o-r','MarkerSize',2,'MarkerFaceColor','r')
            title('Adults')
        end
        end
    
                
        
    
    
%% Plot della mortalita stimata
    figure()
    for j=1:ns
        subplot(2,2,j)
        plot(0:0.1:40,fm_Wood(0:0.1:40,p(:,j),BS))
%         hold on
%         plot(0:0.1:40,fm(j,0:0.1:40),'k') % per confronti
        switch j
        case 1
            title('Eggs')
        case 2
            title('Larvae')
        case 3
            title('Pupae')
         case 4
            title('Adults')
        end
    end   
    figure()
    plot(0:0.1:40,fm_Wood(0:0.1:40,p(:,1),BS),0:0.1:40,fm_Wood(0:0.1:40,p(:,2),BS),0:0.1:40,fm_Wood(0:0.1:40,p(:,3),BS),0:0.1:40,fm_Wood(0:0.1:40,p(:,4),BS),'LineWidth',1)
    legend('eggs','larvae','pupae','adults')
    
 %% Plot delle mortalita 
    figure()
    TT=0:0.1:40;
    nTT=length(TT);
    for kk=1:nprove
    for j=1:ns
        subplot(2,2,j)
        PP=reshape(Pint(:,kk),N,4);
%         plot(0:0.1:40,fm_Wood(0:0.1:40,PP(:,j),BS))
         mp(kk,j,1:nTT)=fm_Wood(TT,PP(:,j),BS);
         for hk=1:nTT
             if mp(kk,j,hk)<0
                 mp(kk,j,hk)=0;
             end
         end
%          hold on
%         plot(0:0.1:40,fm(j,0:0.1:40),'k') % per confronti
%         switch j
%         case 1
%             title('Eggs')
%         case 2
%             title('Larvae')
%         case 3
%             title('Pupae')
%          case 4
%             title('Adults')
%         end
    end     
    end

    for kj=1:nTT
        for j=1:ns

            
            %         plot(0:0.1:40,fm_Wood(0:0.1:40,PP(:,j),BS))
%             mp(kk,1:nTT)=fm_Wood(TT,PP(:,j),BS);
            bandemort=prctile(mp(:,j,kj),[2.5,97.5,50]);
            bande_min(j,kj)=bandemort(1);
            bande_max(j,kj)=bandemort(2);
            bande_med(j,kj)=bandemort(3);
        end
    end
    
%             
         for j=1:ns   
          subplot(2,2,j)
          hold on
        AREA1=area(TT,bande_max(j,:),'LineStyle','none','BaseValue',0.0001);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(TT,bande_min(j,:),'LineStyle','none','BaseValue',0.0001);
        set(AREA6,'FaceColor','w');

        plot(TT,bande_min(j,:),'-k');
        plot(TT,bande_max(j,:),'-k')
        plot(TT,bande_med(j,:),'-k')
            
         
            %         plot(0:0.1:40,fm(j,0:0.1:40),'k') % per confronti
            switch j
                case 1
                    title('Eggs');axis([0 40 0 2]);
                case 2
                    title('Larvae');axis([0 40 0 2]);
                case 3
                    title('Pupae');axis([0 40 0 2]);
                case 4
                    title('Adults');axis([0 40 0 2]);
            end
        end
  
       
%%    
save('prova_pesiuguali100_2020_N5_500simulazioni_nuovofunzionale')    
%%
% bande confidence intervals
% ho fatto le varie prove con le perturbazioni di p, ora devo fare le bande
% al 99% ad agni tempo, per ogni stadio

for i=1:numel(time08)
    for j=1:ns
        CIbande=prctile(Dyn08(i,j,:),[2.5,97.5,50]);
        DynMin(i,j)=CIbande(1);
        DynMed(i,j)=CIbande(3);
        DynMax(i,j)=CIbande(2);
    end
end

% Plot dinamica
%     close all
    figure(2008)
    for j=1:ns-1
        subplot(3,3,j)
        hold on
        %for i=1:nprove

        AREA1=area(time08,DynMax(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(time08,DynMin(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA6,'FaceColor','w');

        plot(time08,DynMin(:,j),'-k')
        plot(time08,DynMax(:,j),'-k')
        plot(time08,DynMed(:,j),'-k')
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('time(days)')
%         ax=gca;
        FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2008/24,rilievi_uova_2008,0*incertezza_2008(1:le2008),incertezza_2008(1:le2008),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2008/24,rilievi_uova_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2008/24,rilievi_larve_2008,0*incertezza_2008(le2008+1:ll2008+le2008),incertezza_2008(le2008+1:ll2008+le2008),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2008/24,rilievi_larve_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2008')
        case 3
            errorbar(dreg_imm_2008/24,rilievi_pupe_2008,0*incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),incertezza_2008(le2008+ll2008+1:ll2008+le2008+lp2008),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2008/24,rilievi_pupe_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2008/24,rilievi_adulti_2008,0*incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),incertezza_2008(lp2008+le2008+ll2008+1:ll2008+le2008+lp2008+la2008),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_ad_2008/24,rilievi_adulti_2008,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
%%
clear DynMin
clear DynMax
clear DynMed
for i=1:numel(time09)
    for j=1:ns
        CIbande=prctile(Dyn09(i,j,:),[2.5,97.5,50]);
        DynMin(i,j)=CIbande(1);
        DynMed(i,j)=CIbande(3);
        DynMax(i,j)=CIbande(2);
    end
end


% Plot dinamica
%     close all
%     figure(2009)
    for j=1:ns-1
        subplot(3,3,j+3)
        hold on
        %for i=1:nprove
        
        AREA1=area(time09,DynMax(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(time09,DynMin(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA6,'FaceColor','w');
        plot(time09,DynMin(:,j),'-k')
        plot(time09,DynMax(:,j),'-k')
        plot(time09,DynMed(:,j),'-k')
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('days')
%         ax=gca;
        FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2009/24,rilievi_uova_2009,0*incertezza_2009(1:le2009),incertezza_2009(1:le2009),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2009/24,rilievi_uova_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2009/24,rilievi_larve_2009,0*incertezza_2009(le2009+1:ll2009+le2009),incertezza_2009(le2009+1:ll2009+le2009),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2009/24,rilievi_larve_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2009')
        case 3
            errorbar(dreg_imm_2009/24,rilievi_pupe_2009,0*incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),incertezza_2009(le2009+ll2009+1:ll2009+le2009+lp2009),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2009/24,rilievi_pupe_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2009/24,rilievi_adulti_2009,0*incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),incertezza_2009(lp2009+le2009+ll2009+1:ll2009+le2009+lp2009+la2009),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_ad_2009/24,rilievi_adulti_2009,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
    %%
 clear DynMin
clear DynMax
clear DynMed
 for i=1:numel(time11)
    for j=1:ns
        CIbande=prctile(Dyn11(i,j,:),[2.5,97.5,50]);
        DynMin(i,j)=CIbande(1);
        DynMed(i,j)=CIbande(3);
        DynMax(i,j)=CIbande(2);
    end
end

% Plot dinamica
%     close all
%     figure(2011)
    for j=1:ns-1
        subplot(3,3,j+6)
        hold on
        %for i=1:nprove
        
        AREA1=area(time11,DynMax(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(time11,DynMin(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA6,'FaceColor','w');
        plot(time11,DynMin(:,j),'-k')
        plot(time11,DynMax(:,j),'-k')
        plot(time11,DynMed(:,j),'-k')
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('time (days)')
%         ax=gca;
        FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2011/24,rilievi_uova_2011,0*incertezza_2011(1:le2011),incertezza_2011(1:le2011),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2011/24,rilievi_uova_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2011/24,rilievi_larve_2011,0*incertezza_2011(le2011+1:ll2011+le2011),incertezza_2011(le2011+1:ll2011+le2011),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2011/24,rilievi_larve_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2011')
        case 3
            errorbar(dreg_imm_2011/24,rilievi_pupe_2011,0*incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),incertezza_2011(le2011+ll2011+1:ll2011+le2011+lp2011),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2011/24,rilievi_pupe_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2011/24,rilievi_adulti_2011,0*incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),incertezza_2011(lp2011+le2011+ll2011+1:ll2011+le2011+lp2011+la2011),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_ad_2011/24,rilievi_adulti_2011,'o:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
    
    
 %%
 % validazione su due anni
 

% bande confidence intervals
% ho fatto le varie prove con le perturbazioni di p, ora devo fare le bande
% al 99% ad agni tempo, per ogni stadio

 clear DynMin
clear DynMax
clear DynMed

for i=1:numel(time10)
    for j=1:ns
        CIbande=prctile(Dyn10(i,j,:),[0.25,97.5,50]);
        DynMin(i,j)=CIbande(1);
        DynMed(i,j)=CIbande(3);
        DynMax(i,j)=CIbande(2);
    end
end

% Plot dinamica
%     close all
    figure(2010)
    for j=1:ns-1
        subplot(3,3,j)
        hold on
        %for i=1:nprove

        AREA1=area(time10,DynMax(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(time10,DynMin(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA6,'FaceColor','w');

        plot(time10,DynMin(:,j),'-k')
        plot(time10,DynMax(:,j),'-k')
        plot(time10,DynMed(:,j),'-k')
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
%         xlabel('time(days)')
%         ax=gca;
        FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2010/24,rilievi_uova_2010,0*incertezza_2010(1:le2010),incertezza_2010(1:le2010),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2010/24,rilievi_uova_2010,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2010/24,rilievi_larve_2010,0*incertezza_2010(le2010+1:ll2010+le2010),incertezza_2010(le2010+1:ll2010+le2010),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2010/24,rilievi_larve_2010,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2010')
        case 3
            errorbar(dreg_imm_2010/24,rilievi_pupe_2010,0*incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),incertezza_2010(le2010+ll2010+1:ll2010+le2010+lp2010),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2010/24,rilievi_pupe_2010,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2010/24,rilievi_adulti_2010,0*incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),incertezza_2010(lp2010+le2010+ll2010+1:ll2010+le2010+lp2010+la2010),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_ad_2010/24,rilievi_adulti_2010,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end
%%
clear DynMin
clear DynMax
clear DynMed
for i=1:numel(time12)
    for j=1:ns-1
        CIbande=prctile(Dyn12(i,j,:),[0.25,97.5,50]);
        DynMin(i,j)=CIbande(1);
        DynMed(i,j)=CIbande(3);
        DynMax(i,j)=CIbande(2);
    end
end


% Plot dinamica
%     close all
%     figure(2009)
    for j=1:ns-1
        subplot(3,3,j+3)
        hold on
        %for i=1:nprove
        
        AREA1=area(time12,DynMax(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA1,'FaceColor',[0.8,0.8,0.8]);
        AREA6=area(time12,DynMin(:,j),'LineStyle','none','BaseValue',0.03);
        set(AREA6,'FaceColor','w');
        plot(time12,DynMin(:,j),'-k')
        plot(time12,DynMax(:,j),'-k')
        plot(time12,DynMed(:,j),'-k')
        %end
        %plot(dreg_2009/24,MU(:,j),'o-k','MarkerSize',2,'MarkerFaceColor','k')
        xlabel('time (days)')
%         ax=gca;
        FontSize=12;
        switch j
        case 1
            errorbar(dreg_imm_2012/24,rilievi_uova_2012,0*incertezza_2012(1:le2012),incertezza_2012(1:le2012),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2012/24,rilievi_uova_2012,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of eggs')
        case 2
            errorbar(dreg_imm_2012/24,rilievi_larve_2012,0*incertezza_2012(le2012+1:ll2012+le2012),incertezza_2012(le2012+1:ll2012+le2012),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2012/24,rilievi_larve_2012,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of larvae')
            title('2012')
        case 3
            errorbar(dreg_imm_2012/24,rilievi_pupe_2012,0*incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),incertezza_2012(le2012+ll2012+1:ll2012+le2012+lp2012),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_imm_2012/24,rilievi_pupe_2012,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of pupae')
        case 4
            errorbar(dreg_ad_2012/24,rilievi_adulti_2012,0*incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),incertezza_2012(lp2012+le2012+ll2012+1:ll2012+le2012+lp2012+la2012),'o:r','MarkerSize',4,'MarkerFaceColor','r')
%             plot(dreg_ad_2012/24,rilievi_adulti_2012,'*:r','MarkerSize',4,'MarkerFaceColor','r')
            ylabel('n. of adults')
        end
    end