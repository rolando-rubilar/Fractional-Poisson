%% Distribucion Poisson Fraccionaria
clear; clc;
rng default;
max_error= 10^(-2); %3
max_error2=10^(-1); %1
s=10000;
tiempo=10;
NPoint=300;
decaimiento=1;
alpha=0.7;
mu=5;
n_simulado = 200;

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

AA=SimFracVec(n_simulado,alpha,mu);
N=length(AA); % numero de muestras AA(0.5,10)

Alpha0=pi/(3*var(log(AA))  + pi^2/6)^0.5;%Alphagrid(915); %0.5
mu0=exp(-Alpha0*( mean(log(AA)) + 0.57721566490153286060));%Mugrid(1325); %10

PRIOR(1,:)=prior;
PRIORm(1,:)=priorM;

tic
for tt=1:100

for t=1:tiempo 

    AA=SimFracVec(n_simulado,alpha,mu);
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
    

    
    [f,Mgrid,bw] = ksdensity(acceptedM, MuMin: (MuMax-MuMin)/(Ngrid-1):MuMax);%,'Bandwidth',VarMu); 
    priorMM=  f/sum(f);
    PRIORm(t+1,:)=priorMM; %% guardand las prior
    [maxValueM ind]= max(f);
    prioM=priorMM;
    gibs(2,t+1)= Mgrid(ind);
    mu0=  Mgrid(ind);
         
end
disp(tt)
media_alpha(tt) = Alpha0;
media_mu(tt) = mu0;
end
toc
