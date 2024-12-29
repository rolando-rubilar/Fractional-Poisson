%% Experimiento KL
clear; clc; 

%% Definiendo grillas
Ngrid=       100;
alphaMin=    0.1;
alphaMax=    0.975;
MuMin=       0.2;   
MuMax=       15;
muestra=     50;

Alphagrid=   alphaMin: (alphaMax-alphaMin)/(Ngrid-1):alphaMax;
Mugrid=      MuMin: (MuMax-MuMin)/(Ngrid-1):MuMax;

alpha=0.7;
mu=5;
X=SimFracVec(muestra*1000,alpha,mu);
LogX=log(X);

Xmuestra=SimFracVec(muestra,alpha,mu);

[f_data,LogXgridMuestra] = ksdensity(log(Xmuestra), min(LogX): (max(LogX)-min(LogX))/(Ngrid-1):max(LogX));%,'Bandwidth',VarAlpha); 
DensidadX=  f_data/sum(f_data);
CDFx=0;
for z=1:length(DensidadX)
    CDFx(z+1)=CDFx(z)+DensidadX(z);
end
%plot(LogXgridMuestra,DensidadX)
%plot(LogXgridMuestra,CDFx(2:end))

for i=1:length(Mugrid)
    for s=1:length(Alphagrid)
        X_sim=SimFracVec(muestra,Alphagrid(s),Mugrid(i));
        
        [f_data_muestra,LogXgridMuestra] = ksdensity(log(X_sim), min(LogX):...
            (max(LogX)-min(LogX))/(Ngrid-1):max(LogX));
        DensidadXmuestra=  f_data_muestra/sum(f_data_muestra);
        %CDFxMuestra=0;
        %for z=1:length(DensidadX)
         %   CDFxMuestra(z+1)=CDFxMuestra(z)+DensidadXmuestra(z);
        %end
    DENSIDAD(s,:,i)= DensidadXmuestra;
    Distancia(s,i)=sum(DensidadX.*log(DensidadX./DensidadXmuestra));
    end
    
    disp(i)
end

[MinMu posMu]=min(Distancia)
[MinAlpha posAlpha]=min(MinMu);

Alpha_est=Alphagrid(posMu(posAlpha));
[MinMu posMuB]=min(Distancia(posMu(posAlpha),:))
Mu_est=Mugrid(posMuB);

figure(1)
plot(LogXgridMuestra,DENSIDAD(1,:,1))
hold on
plot(LogXgridMuestra,DENSIDAD(end,:,end))


minx = min(min(Distancia));
Distancia(Distancia==Inf) = NaN;
maxx = max(nanmax(Distancia));
niveles=100;
levels =  minx:(maxx-minx)/niveles:maxx;

figure(2)
contour(Mugrid,Alphagrid,Distancia,levels)
hold on 
scatter(mu, alpha)
xlabel('\alpha')
ylabel('\mu')

