function swarm2pic(chromosome,varargin)
    global numOfObj dynamic itrCounter colors window;
    if itrCounter == 1 && dynamic  == 1
        set(0,'units','centimeters');
        position    =[0 0 51 19.5];
        h           =figure;
        set(h,'PaperType','A4'); 
        set(h,'PaperUnits','centimeters'); 
        set(h,'paperpositionmode','auto');
        set(h,'PaperPosition',position);
        set(h,'units','centimeters');
        set(h,'position',position);
        hold off;
    end 
    if numel(varargin) == 0;
        if dynamic == 0
            subfunction(chromosome,numOfObj,'rx',6);
        else
            subplot(1,2,1);
            subfunction(chromosome,numOfObj,'rx',6);
            hold off;
            if mod(itrCounter,window) == 0
                subplot(1,2,2);
                subfunction(chromosome,numOfObj,colors{mod(floor(itrCounter/window),18)+1},4);
                hold on;
            end
        end
        title(['NSGA-II Demo',num2str(itrCounter),'Iteration']);
        drawnow;
    else
        hold off;
        figure;
        subfunction(chromosome,numOfObj,colors{mod(floor(itrCounter/window),18)+1},4);
    end
end

function subfunction(chromosome,numOfObj,colorStr,markersize) 
   if numOfObj == 2
        rep_costs=GetCosts(chromosome);
        plot(rep_costs(1,:),rep_costs(2,:),colorStr,'markersize',markersize);
        hold off;
    end
    if numOfObj == 3
        rep_costs=GetCosts(chromosome);
        plot3(rep_costs(1,:),rep_costs(2,:),...
                            rep_costs(3,:),...
                            colorStr,'markersize',markersize);
        view(-243,29);
        axis square;
        grid on;
        hold off;
    end 
end

function Cost = GetCosts(chromosome)
    global nVar numOfObj;
    Cost = chromosome(:,nVar+1:nVar+numOfObj)';
end