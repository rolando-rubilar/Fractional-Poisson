%% Distribucion Poisson Fraccionaria
clear; clc;
%rng default;
max_error= 10^(-2); %3
max_error2=10^(-1); %1
s=10000;
tiempo=10;
NPoint=300;
decaimiento=1;
alpha=0.7;
mu=5;

%% Definiendo grillas
Ngrid=       2000;
alphaMin=    0.1;
alphaMax=    0.975;
MuMin=       0.2;   
MuMax=       15;


Alphagrid=   alphaMin: (alphaMax-alphaMin)/(Ngrid-1):alphaMax;
Mugrid=      MuMin: (MuMax-MuMin)/(Ngrid-1):MuMax;

gibs=      zeros(2,1);


%% Simular 
prior=  ones(1,Ngrid)/Ngrid; % inicializando prior
priorM=  ones(1,Ngrid)/Ngrid; % inicializando prior

AA=SimFracVec(50,alpha,mu);
N=length(AA); % numero de muestras AA(0.5,10)

Alpha0=pi/(3*var(log(AA))  + pi^2/6)^0.5;%Alphagrid(915); %0.5
mu0=exp(-Alpha0*( mean(log(AA)) + 0.57721566490153286060));%Mugrid(1325); %10

PRIOR(1,:)=prior;
PRIORm(1,:)=priorM;

for t=1:tiempo 
tic
    AA=SimFracVec(200,alpha,mu);
    N=length(AA); % numero de muestras AA(0.5,10)

    x=AA;
    X(:,t)=x;
    accepted  = [];
    acceptedM = [];
for z=1:100  
    %computar alpha
    for i=1:s
        proposed_A(i)= randsample(Alphagrid,1,true,prior); % data simulada;
        MU(i)= mu0;
    end
    T=SimFracVec(N,proposed_A,MU); 
    d_A=  pi./(3*var(log(T)) + pi^2/6).^0.5; 
    Va=   pi/(3*var(log(x))  + pi^2/6)^0.5;
    Dif_A= abs(d_A-Va); 
    
    for i=1:s
        if Dif_A(i) < max_error*decaimiento^(t-1)
            accepted(end+1) = proposed_A(i);
        end
        if length(accepted)>NPoint-1
            break;end
    end
    if length(accepted)>NPoint-1
            break;end
end       
   
    %VarAlpha=var(accepted)*2;
    [f1,Agrid] = ksdensity(accepted, alphaMin: (alphaMax-alphaMin)/(Ngrid-1):alphaMax);%,'Bandwidth',VarAlpha); 
    priorA=  f1/sum(f1);
    PRIOR(t+1,:)=priorA; %% guardand las prior
    [maxValueA ind]= max(f1);
    prior=priorA;
    gibs(1,t+1)= Agrid(ind);
    Alpha0= Agrid(ind);
    
  
for z=1:100         % computar Mu
    for i=1:s      
        proposed_M(i)= randsample(Mugrid,1,true,priorM); % data simulada;
        ALPHA(i)=Alpha0;
    end   
    T=SimFracVec(N,ALPHA,proposed_M); 
    d_M= exp(-Alpha0*( mean(log(T)) + 0.57721566490153286060)); 
    Vm=     exp(-Alpha0*( mean(log(x)) + 0.57721566490153286060));
    Dif_M= abs(d_M-Vm);
    for i=1:s
        if Dif_M(i) < max_error2*decaimiento^(t-1)
            acceptedM(end+1) = proposed_M(i);
        end
        if length(acceptedM)>NPoint-1
            break;end
    end
    if length(acceptedM)>NPoint-1
            break;end
end
     
     %max_error2=max_error2/2;
   

    %VarMu=var(acceptedM)*2;
    [f,Mgrid,bw] = ksdensity(acceptedM, MuMin: (MuMax-MuMin)/(Ngrid-1):MuMax);%,'Bandwidth',VarMu); 
    priorMM=  f/sum(f);
    PRIORm(t+1,:)=priorMM; %% guardand las prior
    [maxValueM ind]= max(f);
    prioM=priorMM;
    gibs(2,t+1)= Mgrid(ind);
    mu0=  Mgrid(ind);
 toc        
end


% for i=1:length(PRIOR(:,1))
%     figure(1)
%     plot(Alphagrid,PRIOR(i,:))
%     title('\alpha')
%     
%     figure(2)
%     plot(Mugrid,PRIORm(i,:))
%     title('\mu')
%     pause
% end


cdfA(1)=0;
cdfM(1)=0;
for i=1:length(priorMM)
    cdfA(i+1)=cdfA(i)+priorA(i);
    cdfM(i+1)=cdfM(i)+priorMM(i); 
    if cdfA(i) > 0.0001
        InflimCDFa(i)=i;
    else
        InflimCDFa(i)=10000;
    end
    if cdfA(i) > 0.9999
        SuplimCDFa(i)=i;
    else
        SuplimCDFa(i)=10000;
    end
    if cdfM(i) > 0.0001
        InflimCDFm(i)=i;
    else
        InflimCDFm(i)=10000;
    end
    if cdfM(i) > 0.9999
        SuplimCDFm(i)=i;
    else
        SuplimCDFm(i)=10000;
    end
end
    
if min(SuplimCDFm)==10000
    SuplimCDFm(end)=2000;
end

figure(10)
subplot(2,1,1);
plot(Alphagrid,PRIOR(end,:),'k','linewidth',4)
hold on 
scatter(0.6745,0,150,'filled','d','r')

scatter(alpha,0,150,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
          
xlabel('\alpha')
ylabel('P(\alpha)')
legend('Posterior','Kullback-Leibler','Real Value')
ylim([0 0.03])
xlim([Alphagrid(min(InflimCDFa)) Alphagrid(min(SuplimCDFa))])

subplot(2,1,2);
plot(Mugrid,PRIORm(end,:),'k','linewidth',4)
hold on
scatter(4.53,0,150,'filled','d','r')
scatter(mu,0,150,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
xlabel('\mu')
ylabel('P(\mu)')
legend('Posterior','Kullback-Leibler','Real Value')
ylim([0 0.03])
xlim([Mugrid(min(InflimCDFm)) Mugrid(min(SuplimCDFm))])

%% Tests





%% plots alternativas
% 
% % 
% figure(10)
% subplot(2,1,1);
% plot(Alphagrid,ZZ1(1,:),'k','linewidth',2)
% hold on
% plot(Alphagrid,ZZ1(2,:),'--k','linewidth',2)
% scatter(alpha,0,150,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5)
% 
% xlabel('\alpha')
% ylabel('P(\alpha)')
% legend('\rho evolutive','\rho fix')
% 
% subplot(2,1,2);
% plot(Mugrid,ZZ2(1,:),'k','linewidth',2)
% hold on
% plot(Mugrid,ZZ2(2,:),'--k','linewidth',2)
% scatter(mu,0,150,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5)
% xlabel('\mu')
% ylabel('P(\mu)')
% legend('\rho evolutive','\rho fix')
% 
% % 
% % 
