%% FPP Simulacion
function SimFracVec= SimFracVec(N,v,mu)
Nrep=length(v);
V=repmat(v,N,1);
MU=repmat(mu,N,1);
U1=rand(N,Nrep);
U2=rand(N,Nrep);
U3=rand(N,Nrep);
pp= (abs(log(U1)).^(1./V))./(MU.^(1./V)); %check
S1=sin(V*pi.*U2);
S2=(sin((1-V)*pi.*U2)).^(1./V-1);
S3=(sin(pi*U2)).^(1./V);
S4=abs(log(U3)).^(1./V-1);

SimFracVec=pp.*S1.*S2./(S3.*S4);

