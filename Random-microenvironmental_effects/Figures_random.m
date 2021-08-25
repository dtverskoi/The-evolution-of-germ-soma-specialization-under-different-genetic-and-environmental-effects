%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code visualize the results of our numerical simulations for the case
% of random microenvironmental effects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% phenotypes b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gg=1:3
    
    figure
    set(gcf, 'Position',  [200, 0, 1200, 1000])
    
    for lN=1:3
        N=4*2^lN;
        if lN==1
            ga=2;
        else
            if lN==2
                ga=3/2;
            else
                ga=2/3;
            end
        end
        for lb=1:3
            beta=1+2*(lb-1)+(lb-1)*(lb-2);
            if lb==3
                beta=6;
            end
            subplot(3,3,3*(lN-1)+lb)
            avp=zeros(3,9,4);
            stdp=zeros(3,9,4);
            for lG=1:gg
                G=2*4^(lG-1);
                
                for lk=1:9
                    k1=0.1*lk;
                    Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                    p=(sort(Q(:,1:4)'))';
                    avp(lG,lk,:)=mean(p,1);
                    stdp(lG,lk,:)=std(p);
                end
            end
            bl=0.1*ones(1,9);
            bh=(1-0.1^ga)^(1/ga)*ones(1,9);
            errorbar(linspace(0.1,0.9,9),avp(gg,:,1),stdp(gg,:,1),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
            hold on
            errorbar(linspace(0.1,0.9,9),avp(gg,:,2),stdp(gg,:,2),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
            errorbar(linspace(0.1,0.9,9),avp(gg,:,3),stdp(gg,:,3),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
            errorbar(linspace(0.1,0.9,9),avp(gg,:,4),stdp(gg,:,4),'-s','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',2,'Color',[0.5 0.5 0.5])
            plot(linspace(0.1,0.9,9),bl,'--','color',[0.5 0.5 0.5],'LineWidth',2)
            plot(linspace(0.1,0.9,9),bh,'--','color',[0.5 0.5 0.5],'LineWidth',2)
            hold off
            
            set(gca,'FontSize',20)
            xlim([0 1])
            ylim([0 1])
            
            if lN==3
                xlabel('k','Fontsize',30)
            end
            if lb==1
                ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}b'})
            end
            if lN==1
                title(['\beta=', num2str(beta)],'Fontsize',30)
            end
        end
    end
    
    print(['bG=',num2str(G)],'-dpng')
end
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% total fecundity B and total activity A %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

figure
set(gcf, 'Position',  [200, 0, 1200, 1000])

for lN=1:3
    N=4*2^lN;
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        subplot(3,3,3*(lN-1)+lb)
        avp=zeros(3,9);
        stdp=zeros(3,9);
        avpn=zeros(3,9);
        stdpn=zeros(3,9);
        for lG=1:3
            G=2*4^(lG-1);
            
            for lk=1:9
                k1=0.1*lk;
                Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                p=Q(:,11);
                pn=Q(:,12);
                avp(lG,lk)=mean(p);
                stdp(lG,lk)=std(p);
                avpn(lG,lk)=mean(pn);
                stdpn(lG,lk)=std(pn);
            end
        end
        errorbar(linspace(0.1,0.9,9),avp(1,:),stdp(1,:),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        hold on
        errorbar(linspace(0.1,0.9,9),avp(2,:),stdp(2,:),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avp(3,:),stdp(3,:),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        errorbar(linspace(0.1,0.9,9),avpn(1,:),stdpn(1,:),'-o','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        errorbar(linspace(0.1,0.9,9),avpn(2,:),stdpn(2,:),'-o','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avpn(3,:),stdpn(3,:),'-o','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        hold off

        set(gca,'FontSize',20)
        xlim([0 1])
        
        if lN==3
            xlabel('k','Fontsize',30)
        end
        if lb==1
            ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}A, B'})
        end
        if lN==1
            title(['\beta=', num2str(beta)],'Fontsize',30)
        end
    end
end

print('AB','-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% share of specialized cells p %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Position',  [200, 0, 1200, 1000])

for lN=1:3
    N=4*2^lN;
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        subplot(3,3,3*(lN-1)+lb)
        avp=zeros(3,9);
        stdp=zeros(3,9);
        for lG=1:3
            G=2*4^(lG-1);
            
            for lk=1:9
                k1=0.1*lk;
                Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                p=Q(:,14);
                avp(lG,lk)=mean(p);
                stdp(lG,lk)=std(p);
            end
        end
        errorbar(linspace(0.1,0.9,9),avp(1,:),stdp(1,:),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        hold on
        errorbar(linspace(0.1,0.9,9),avp(2,:),stdp(2,:),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avp(3,:),stdp(3,:),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        hold off
        
        if (lb==2)&&(lN==1)
            legend('G=2','G=8','G=32','FontSize',25,'Location','northeast')
        end
        
        set(gca,'FontSize',20)
        xlim([0 1])
        ylim([0 1])
        
        if lN==3
            xlabel('k','Fontsize',30)
        end
        if lb==1
            ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}p'})
        end
        if lN==1
            title(['\beta=', num2str(beta)],'Fontsize',30)
        end
    end
end

print('p','-dpng') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% share of reproductive cells pg %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Position',  [200, 0, 1200, 1000])

for lN=1:3
    N=4*2^lN;
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        subplot(3,3,3*(lN-1)+lb)
        avp=zeros(3,9);
        stdp=zeros(3,9);
        for lG=1:3
            G=2*4^(lG-1);
            
            for lk=1:9
                k1=0.1*lk;
                Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                p=Q(:,15);
                avp(lG,lk)=mean(p);
                stdp(lG,lk)=std(p);
            end
        end
        errorbar(linspace(0.1,0.9,9),avp(1,:),stdp(1,:),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        hold on
        errorbar(linspace(0.1,0.9,9),avp(2,:),stdp(2,:),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avp(3,:),stdp(3,:),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        hold off

        if (lb==2)&&(lN==1)
            legend('G=2','G=8','G=32','FontSize',25,'Location','northeast')
        end
        
        set(gca,'FontSize',20)
        xlim([0 1])
        ylim([0 1])
        
        if lN==3
            xlabel('k','Fontsize',30)
        end
        if lb==1
            ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}p_g'})
        end
        if lN==1
            title(['\beta=', num2str(beta)],'Fontsize',30)
        end
    end
end

print('pg','-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% equilibrium number of colonies %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Position',  [200, 0, 1200, 1000])

for lN=1:3
    N=4*2^lN;
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        subplot(3,3,3*(lN-1)+lb)
        avp=zeros(3,9);
        stdp=zeros(3,9);
        for lG=1:3
            G=2*4^(lG-1);
            
            for lk=1:9
                k1=0.1*lk;
                Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                p=Q(:,9);
                avp(lG,lk)=mean(p);
                stdp(lG,lk)=std(p);
            end
        end
        errorbar(linspace(0.1,0.9,9),avp(1,:),stdp(1,:),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        hold on
        errorbar(linspace(0.1,0.9,9),avp(2,:),stdp(2,:),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avp(3,:),stdp(3,:),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        hold off
        
        if (lb==2)&&(lN==1)
            legend('G=2','G=8','G=32','FontSize',25,'Location','northwest')
        end
        
        set(gca,'FontSize',20)
        xlim([0 1])
        %ylim([1 4])
        
        if lN==3
            xlabel('k','Fontsize',30)
        end
        if lb==1
            ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}N'})
        end
        if lN==1
            title(['\beta=', num2str(beta)],'Fontsize',30)
        end
    end
end

print('M','-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% number of different cell types %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Position',  [200, 0, 1200, 1000])

for lN=1:3
    N=4*2^lN;
    for lb=1:3
        beta=1+2*(lb-1)+(lb-1)*(lb-2);
        if lb==3
            beta=6;
        end
        subplot(3,3,3*(lN-1)+lb)
        avp=zeros(3,9);
        stdp=zeros(3,9);
        for lG=1:3
            G=2*4^(lG-1);
            
            for lk=1:9
                k1=0.1*lk;
                Q=readmatrix(['N' num2str(N) 'G' num2str(G) 'be' num2str(beta) 'k' num2str(10*k1)]);
                p=Q(:,16);
                avp(lG,lk)=mean(p);
                stdp(lG,lk)=std(p);
            end
        end
        errorbar(linspace(0.1,0.9,9),avp(1,:),stdp(1,:),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2,'Color','blue')
        hold on
        errorbar(linspace(0.1,0.9,9),avp(2,:),stdp(2,:),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2,'Color','red')
        errorbar(linspace(0.1,0.9,9),avp(3,:),stdp(3,:),'-s','MarkerSize',10,'MarkerEdgeColor',[0 230/255 115/255],'MarkerFaceColor',[0 230/255 115/255],'LineWidth',2,'Color',[0 230/255 115/255])
        hold off
        
        if (lb==2)&&(lN==1)
            legend('G=2','G=8','G=32','FontSize',25,'Location','northeast')
        end
        
        set(gca,'FontSize',20)
        xlim([0 1])
        ylim([1 4])
        
        if lN==3
            xlabel('k','Fontsize',30)
        end
        if lb==1
            ylabel({['\fontsize{30}S=',num2str(N)],'\fontsize{30}M'})
        end
        if lN==1
            title(['\beta=', num2str(beta)],'Fontsize',30)
        end
    end
end

print('N','-dpng')
