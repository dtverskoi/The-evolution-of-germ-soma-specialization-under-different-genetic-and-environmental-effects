%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is a Matlab implementation of the model for the case of
% several microenvironmentsl gradients.
% This function automatically simulates the model 20 runs and shows the
% results on the average equilibrium characteristics of colonies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q] = Main_gradients(lN,lG,lb,lk)
Runs=20;                               %number of runs
sch=0;

% Parameters
N=4*2^lN;                              %colony size
G=2*4^(lG-1);                          %number of genes
beta=1+2*(lb-1)+(lb-1)*(lb-2);         %the strength of competition among colonies
if lb==3
    beta=6;
end
k1=0.1*lk;                             %cost of fecundity
d=10;                                  %gradient constant
%rates for different gradients
r1=0.2;
r2=0.4;
r3=0.6;
r4=0.8;
k2=1-k1;                               %cost of activity
M0=1;                                  %initial number of groups
mu=0.001;                              %mutation rate
sigma=0.05;                            %mutation st. dv.
T=1000000;                             %number of time steps
skip=10000;                            %for sampling/saving data
tau=T/skip;
eksp=4;                                %parameter in sigmoid function
s=4;                                   %number of prototypes
C=1000;                                %the amount of the resource

%Curvature of the trade-off
if N==8
    gamma=2;
else
    if N==16
        gamma=3/2;
    else
        gamma=2/3;
    end
end

Q=NaN(Runs,2*s+8);                     %output
run=1;

while run<=Runs
    % initialize random number generator
    if 1
        rng('shuffle');
        seed=rng;
        save seed;
    else
        load seed;
        rng(seed)
    end
    
    %microenvironment
    genematrix=zeros(s,G);
    for si=1:s
        for sd=1:G/8
            genematrix(si,sd)=d*r1^(si-1);
        end
        for sd=G/8+1:2*G/8
            genematrix(si,sd)=-d*r1^(si-1);
        end
        for sd=2*G/8+1:3*G/8
            genematrix(si,sd)=d*r2^(si-1);
        end
        for sd=3*G/8+1:4*G/8
            genematrix(si,sd)=-d*r2^(si-1);
        end
        for sd=4*G/8+1:5*G/8
            genematrix(si,sd)=d*r3^(si-1);
        end
        for sd=5*G/8+1:6*G/8
            genematrix(si,sd)=-d*r3^(si-1);
        end
        for sd=6*G/8+1:7*G/8
            genematrix(si,sd)=d*r4^(si-1);
        end
        for sd=7*G/8+1:8*G/8
            genematrix(si,sd)=-d*r4^(si-1);
        end
    end
    
    % Program
    % Initial values
    M=M0;                                            % initial number of colonies
    genes=0.5*ones(G,M);                             % initial genes
    b=1./(1.+exp(-eksp*genematrix*genes/sqrt(G)));   % initial fertilities
    v=(1-b.^gamma).^(1/gamma);                       % initial viabilities
    btypes=zeros(s,tau);
    vtypes=zeros(s,tau);
    MMb=zeros(1,tau);
    MM=zeros(1,tau);
    BB=zeros(1,tau);
    VV=zeros(1,tau);
    p=zeros(1,tau);
    pg=zeros(1,tau);
    Suv=zeros(1,tau);
    
    for t=1:T 
        % Competition betwen colonies
        B=(N/s)*sum(b,1);
        V=(N/s)*sum(v,1);
        V1=V.^beta;
        c=C*V1/sum(V1);                              % resources per colony
        X=k1*B+k2*V;                                 % constraints per colony
        
        % Statistics
        if  rem(t,skip)==0
            t/skip;
            MMb(t/skip)=M;
            Suv(t/skip)=mean(min(1,c./X));
        end
        
        % viability selection of colonies
        R=rand(1,M);
        b(:,(R<(X-c)./X))=[];                        % killing colonies
        genes(:,(R<(X-c)./X))=[];                    % killing genes of killed colonies
        M=size(b,2);                                 % counting surviving colonies
        
        v=(1-b.^gamma).^(1/gamma);                   % new viabilities
        avB=(N/s)*mean(sum(b,1));
        avV=(N/s)*mean(sum(v,1));
        
        % Statistics
        if  rem(t,skip)==0
            btypes(:,t/skip)=mean(b,2);
            vtypes(:,t/skip)=mean(v,2);
            MM(t/skip)=M;
            BB(t/skip)=avB;
            VV(t/skip)=avV;
            h=zeros(s,M);
            h1=zeros(s,M);
            h(b<=0.1)=1;
            h(v<=0.1)=1;
            h1(v<=0.1)=1;
            hh=sum(h,1);
            hh1=sum(h1,1);
            p(t/skip)=mean(hh)/s;
            pg(t/skip)=mean(hh1)/s;
        end
        
        % fertility selection of cells
        genesnew=[];
        Mnew=0;
        for kk=1:N/s
            R=rand(s,M);
            f=find(R<b);                          % fertile cells
            
            % fertilities of cells in new colonies
            [~,J]=ind2sub([s,M],f);               % colonies to which fertile cells belong
            Mnew=Mnew+length(f);                  % number of fertile cells (= number of new colonies)
            genesnew=[genesnew genes(:,J)];       % new genes matrix before mutations
        end
        M=Mnew;
        genes=genesnew;
        
        % mutations
        R=rand(G,M);
        mutants=find(R<mu);                                             % mutant genes
        genes(mutants)=genes(mutants)+sigma*randn(length(mutants),1);   % new genes
        genes(genes<0)=0;
        genes(genes>1)=1;
        
        b=1./(1.+exp(-eksp*genematrix*genes/sqrt(G)));                  % new fertilities
        v=(1-b.^gamma).^(1/gamma);                                      % new viabilities
        
    end
    
    if M>0
        % types clasterisation
        Linemin=(min((btypes(:,[3*tau/4:tau]))'))';
        Linemax=(max((btypes(:,[3*tau/4:tau]))'))';
        P=zeros(s,s);
        numbertypes=s;
        for i=1:s
            for j=1:s
                if (~((Linemin(i)>Linemax(j)+0.05)||(Linemax(i)<Linemin(j)-0.05)))&&(i~=j)
                    P(i,j)=1;
                end
            end
        end
        PP=P+P^2+P^3;
        for i=1:s
            PP(i,i)=0;
        end
        i=1;
        while sum(sum(PP))>0
            for j=1:s
                if PP(i,j)>0
                    numbertypes=numbertypes-1;
                    PP(j,:)=zeros(1,s);
                end
            end
            PP(i,:)=zeros(1,s);
            i=i+1;
        end
        
        Q(run,:)=[mean(btypes(:,3*tau/4:tau),2)' mean(vtypes(:,3*tau/4:tau),2)' mean(MM(3*tau/4:tau)) mean(MMb(3*tau/4:tau)) mean(BB(3*tau/4:tau)) mean(VV(3*tau/4:tau)) mean(Suv(3*tau/4:tau)) mean(p(3*tau/4:tau)) mean(pg(3*tau/4:tau)) numbertypes];
        run=run+1;
        sch=sch+1;
    end
    
end

dlmwrite(['G4N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) '.txt'],Q)
end