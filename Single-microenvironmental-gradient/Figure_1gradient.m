%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code visualize the results of our numerical simulations for the case
% of a single microenvironmental gradient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% phenotypes b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    if lN==1
        ga=2;
    else
        ga=2/3;
    end
    
    for lG=1:2
        G=2*16^(lG-1);
        
        figure
        set(gcf, 'Position',  [200, 0, 1200, 1000])
        
        for lb=1:3
            beta=1+2*(lb-1)+(lb-1)*(lb-2);
            if lb==3
                beta=6;
            end
            
            for lk=1:3
                k1=0.2+0.3*(lk-1);
                
                subplot(3,3,3*(lk-1)+lb)
                avp=zeros(2,4,4);
                stdp=zeros(2,4,4);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);
                    p=Q(:,1:4);
                    avp(lG,lR,:)=mean(p,1);
                    stdp(lG,lR,:)=std(p);
                end
                bl=0.1*ones(1,9);
                bh=(1-0.1^ga)^(1/ga)*ones(1,9);
                errorbar(linspace(0.2,0.8,4),avp(lG,:,1),stdp(lG,:,1),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
                hold on
                errorbar(linspace(0.2,0.8,4),avp(lG,:,2),stdp(lG,:,2),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
                errorbar(linspace(0.2,0.8,4),avp(lG,:,3),stdp(lG,:,3),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
                errorbar(linspace(0.2,0.8,4),avp(lG,:,4),stdp(lG,:,4),'-s','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',2,'Color',[0.5 0.5 0.5])
                plot(linspace(0.1,0.9,9),bl,'--','color',[0.5 0.5 0.5],'LineWidth',2)
                plot(linspace(0.1,0.9,9),bh,'--','color',[0.5 0.5 0.5],'LineWidth',2)
                hold off

                set(gca,'FontSize',20)
                xlim([0 1])
                ylim([0 1])
                
                if lk==3
                    xlabel('r','Fontsize',30)
                end
                if lb==1
                    ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}b'})
                end
                if lk==1
                    title(['\beta=', num2str(beta)],'Fontsize',30)
                end
            end
        end
        
        print(['bS=',num2str(N), 'G=', num2str(G)],'-dpng')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% total fecundity B and total activity A %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        
        for lk=1:3
            k1=0.2+0.3*(lk-1);
            
            subplot(3,3,3*(lk-1)+lb)
            avp=zeros(2,4);
            stdp=zeros(2,4);
            for lG=1:2
                G=2*16^(lG-1);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);   
                    p=Q(:,11);
                    pn=Q(:,12);
                    avp(lG,lR)=mean(p);
                    stdp(lG,lR)=std(p);
                    avpn(lG,lR)=mean(pn);
                    stdpn(lG,lR)=std(pn);
                end
            end
            errorbar(linspace(0.2,0.8,4),avp(1,:),stdp(1,:),'-s','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            hold on
            errorbar(linspace(0.2,0.8,4),avp(2,:),stdp(2,:),'-s','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            errorbar(linspace(0.2,0.8,4),avpn(1,:),stdpn(1,:),'-o','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            errorbar(linspace(0.2,0.8,4),avpn(2,:),stdpn(2,:),'-o','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            hold off
            
            set(gca,'FontSize',20)
            xlim([0 1])
            
            if lk==3
                xlabel('r','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}A,B'})
            end
            if lk==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['ABS=',num2str(N)],'-dpng')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% share of specialized cells p %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        
        for lk=1:3
            k1=0.2+0.3*(lk-1);
            
            subplot(3,3,3*(lk-1)+lb)
            avp=zeros(2,4);
            stdp=zeros(2,4);
            for lG=1:2
                G=2*16^(lG-1);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);   
                    p=Q(:,14);
                    avp(lG,lR)=mean(p);
                    stdp(lG,lR)=std(p);
                end
            end
            errorbar(linspace(0.2,0.8,4),avp(1,:),stdp(1,:),'-s','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            hold on
            errorbar(linspace(0.2,0.8,4),avp(2,:),stdp(2,:),'-s','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            hold off
            
            if (lb==2)&&(lk==1)
                legend('G=2','G=32','FontSize',25,'Location','northeast')
            end
            
            set(gca,'FontSize',20)
            xlim([0 1])
            ylim([0 1])
            
            if lk==3
                xlabel('r','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}p'})
            end
            if lk==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['pS=',num2str(N)],'-dpng')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% share of reproductive cells pg %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        
        for lk=1:3
            k1=0.2+0.3*(lk-1);
            
            subplot(3,3,3*(lk-1)+lb)
            avp=zeros(2,4);
            stdp=zeros(2,4);
            for lG=1:2
                G=2*16^(lG-1);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);
                    p=Q(:,15);
                    avp(lG,lR)=mean(p);
                    stdp(lG,lR)=std(p);
                end
            end
            errorbar(linspace(0.2,0.8,4),avp(1,:),stdp(1,:),'-s','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            hold on
            errorbar(linspace(0.2,0.8,4),avp(2,:),stdp(2,:),'-s','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            hold off
            
            if (lb==2)&&(lk==1)
                legend('G=2','G=32','FontSize',25,'Location','northeast')
            end
            
            set(gca,'FontSize',20)
            xlim([0 1])
            ylim([0 1])
            
            if lk==3
                xlabel('r','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}p_g'})
            end
            if lk==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['pgS=',num2str(N)],'-dpng')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% equilibrium number of colonies %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        
        for lk=1:3
            k1=0.2+0.3*(lk-1);
            
            subplot(3,3,3*(lk-1)+lb)
            avp=zeros(2,4);
            stdp=zeros(2,4);
            for lG=1:2
                G=2*16^(lG-1);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);

                    p=Q(:,9);
                    avp(lG,lR)=mean(p);
                    stdp(lG,lR)=std(p);
                end
            end
            errorbar(linspace(0.2,0.8,4),avp(1,:),stdp(1,:),'-s','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            hold on
            errorbar(linspace(0.2,0.8,4),avp(2,:),stdp(2,:),'-s','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            hold off
            
            if (lb==2)&&(lk==1)
                legend('G=2','G=32','FontSize',25,'Location','northeast')
            end
            
            set(gca,'FontSize',20)
            xlim([0 1])
            if lN==1
                ylim([150 310])
            else
                ylim([45 160])
            end
            
            if lk==3
                xlabel('r','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}N'})
            end
            if lk==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['MS=',num2str(N)],'-dpng')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% number of different cell types %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lN=1:2
    N=2*4^lN;
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        
        for lk=1:3
            k1=0.2+0.3*(lk-1);
            
            subplot(3,3,3*(lk-1)+lb)
            avp=zeros(2,4);
            stdp=zeros(2,4);
            for lG=1:2
                G=2*16^(lG-1);
                
                for lR=1:4
                    r=0.2*lR;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1) 'r' num2str(10*r)]);
                    p=Q(:,16);
                    avp(lG,lR)=mean(p);
                    stdp(lG,lR)=std(p);
                end
            end
            errorbar(linspace(0.2,0.8,4),avp(1,:),stdp(1,:),'-s','MarkerSize',15,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3,'Color','blue')
            hold on
            errorbar(linspace(0.2,0.8,4),avp(2,:),stdp(2,:),'-s','MarkerSize',15,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',3,'Color',[0 230/255 115/255])
            hold off
            
            if (lb==2)&&(lk==1)
                legend('G=2','G=32','FontSize',25,'Location','northeast')
            end
            
            set(gca,'FontSize',20)
            xlim([0 1])
            ylim([1 4])
            
            if lk==3
                xlabel('r','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}k=',num2str(k1)],'\fontsize{30}M'})
            end
            if lk==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['NS=',num2str(N)],'-dpng')
end