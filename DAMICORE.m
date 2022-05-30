classdef DAMICORE < FastNewman
    properties
        Clusters
        ColorList 
    end
    %%
    methods
        function obj = DAMICORE(Dataset,K)
            obj@FastNewman(Dataset,K)
            obj = obj.colorgen();
            obj = obj.setclusters();
            obj = obj.setnewickcolor();
            obj.printnewick();
        end
        function Clusters = getclusterlist(obj,Klist)
            ClusteList = zeros(size(obj.LeafList,1),numel(Klist));
            for k=1:numel(Klist)
                obj = obj.getCommunities(Klist(k));
                ClusteList(:,k) = obj.getclusters();
            end
            Clusters = ClusteList;
        end
        function getmaxmodularity(obj)
        end
        function plotclusters(obj)
            G = graph(obj.AdjacencyMatrix);
            obj.plotgraph(G,obj.NodeCommunities);
        end
        function plotnewick(obj)
        %% Call Newick    
        !"C:\FigTree v1.4.3\figtree.jar" -graphic PNG PhyTree.nexus PhyTree.png
        [y,z]=imread('PhyTree','png','BackgroundColor', [1,1,1]);
        %% Plot Tree PNG
        imshow(y);hold on
        Cmap = linspecer(max(obj.NodeCommunities));
        maxC = size(Cmap,1);
        for k=1:max(obj.NodeCommunities)
            getColor = Cmap(rem(k-1,maxC)+1,:);
            Nodes = (obj.NodeCommunities==k);
            h(k) = plot([nan,nan],[nan,nan],'o','MarkerEdgeColor',getColor,'MarkerFaceColor',getColor);
            leg(k) = string(['C_{',int2str(k),'}']);
        end
        %% Preparing Layout.
        T = title('Graph Communities','Units','normalized');
        hleg = legend(h,leg,'Location','northeastoutside','Orientation','vertical');
        htitle = get(hleg,'Title');
        set(htitle,'String','Clusters');
        end
        function plotmodularity(obj)
            Rep  = sortrows(obj.ReportOptimization,'N');
            N     = Rep.N;
            Q     = Rep.Modularity;
            Color = linspecer(5);
            h1=plot(N,Q,'-o','linewidth',2,'Color',Color(1,:),...
                'MarkerFaceColor',Color(1,:),'MarkerSize',4);grid on;hold on;
            h2=xline(N(Q==max(Q)),'-.',"N = "+int2str(N(Q==max(Q))),'Color',Color(3,:),'LineWidth',1,'LabelOrientation','horizontal');
            h3=plot(N(Q==max(Q)),max(Q),'s','Color',Color(2,:),'MarkerFaceColor',Color(2,:),'LineWidth',2,'MarkerSize',6);
            hleg = legend([h1,h3,h2],'Q','Q_{Max}','N_{opt}','location','northeastoutside');
            title(hleg,'arg max Q_{N}');
            ylabel('Modularity (Q)');
            xlim([0,2*N(Q==max(Q))]);
            ylim([0,1]);grid on
        end
        function Clusters = getclusters(obj)
            Clusters(obj.LeafList,1) = obj.NodeCommunities(obj.isLeaf);
        end
    end
    %% Auxiliar methods
    methods
        function obj= setclusters(obj)
            obj.Clusters = getclusters(obj);
        end
    end
    %% View Clusters Newick Methods
    methods
        function obj = setnewickcolor(obj)
            for k=1:size(obj.NodeCommunities,1)
                Cluster  = obj.NodeCommunities(k);
                Node     = string(obj.NodeId(k));
                Color    = obj.ColorList(Cluster);
                ColorNewick = obj.newickcolorformat(Color);
                NewString = "'"+Node+"'" + ColorNewick;
                Rep = regexp(obj.NewickFile,char(Node+'(\D)'));
                
                s1 = obj.NewickFile(1:Rep-1);
                s2 = obj.NewickFile(Rep+numel(char(Node)):end);
                obj.NewickFile = [s1,char(NewString),s2];
            end
        end
        function printnewick(obj)
            %%
            f= fopen('PhyTree.nexus','w+');
            fprintf(f,'#NEXUS\r\n');
            fprintf(f,'begin taxa;\r\n');
            fprintf(f,['\tdimensions ntax=' num2str(size(obj.NodeId,1)) ';\r\n']);
            fprintf(f,'\ttaxlabels\r\n');
            for k = 1: size(obj.NodeId,1)
                fprintf(f,['\t','''',obj.NodeId{k},'''','\r\n']);
            end
            fprintf(f,';\nend\n\r\nbegin trees;\r\n');
            fprintf(f,['\ttree tree_1=[&R]' char(obj.NewickFile) '\r\n']);
            fprintf(f,'end;\r\n\r\n');
            fprintf(f,'begin figtree;\r\n');
            fprintf(f,'\tset appearance.backgroundColorAttribute="Default";\r\n');
            fprintf(f,'\tset appearance.backgroundColour=#ffffff;\r\n');
            fprintf(f,'\tset appearance.branchColorAttribute="User selection";\r\n');
            fprintf(f,'\tset appearance.branchColorGradient=false;\r\n');
            fprintf(f,'\tset appearance.branchLineWidth=2.0;\r\n');
            fprintf(f,'\tset appearance.branchMinLineWidth=0.0;\r\n');
            fprintf(f,'\tset appearance.branchWidthAttribute="Fixed";\r\n');
            fprintf(f,'\tset appearance.foregroundColour=#000000;\r\n');
            fprintf(f,'\tset appearance.hilightingGradient=false;\r\n');
            fprintf(f,'\tset appearance.selectionColour=#2d3680;\r\n');
            fprintf(f,'\tset branchLabels.colorAttribute="User selection";\r\n');
            fprintf(f,'\tset branchLabels.displayAttribute="Branch times";\r\n');
            fprintf(f,'\tset branchLabels.fontName="Adobe Caslon Pro";\r\n');
            fprintf(f,'\tset branchLabels.fontSize=8;\r\n');
            fprintf(f,'\tset branchLabels.fontStyle=0;\r\n');
            fprintf(f,'\tset branchLabels.isShown=false;\r\n');
            fprintf(f,'\tset branchLabels.significantDigits=4;\r\n');
            fprintf(f,'\tset layout.expansion=0;\r\n');
            fprintf(f,'\tset layout.layoutType="POLAR";\r\n');
            fprintf(f,'\tset layout.zoom=0;\r\n');
            fprintf(f,'\tset legend.attribute=null;\r\n');
            fprintf(f,'\tset legend.fontSize=10.0;\r\n');
            fprintf(f,'\tset legend.isShown=false;\r\n');
            fprintf(f,'\tset legend.significantDigits=4;\r\n');
            fprintf(f,'\tset nodeBars.barWidth=4.0;\r\n');
            fprintf(f,'\tset nodeBars.displayAttribute=null;\r\n');
            fprintf(f,'\tset nodeBars.isShown=false;\r\n');
            
            fprintf(f,'\tset nodeShapeExternal.isShown=true;\r\n');
            fprintf(f,'\tset nodeShapeExternal.shapeType=Diamond;\r\n');
            fprintf(f,'\tset nodeShapeExternal.size=3.0;\r\n');
            
            
            fprintf(f,'\tset nodeLabels.colorAttribute="User selection";\r\n');
            fprintf(f,'\tset nodeLabels.displayAttribute="Node ages";\r\n');
            fprintf(f,'\tset nodeLabels.fontName="Adobe Caslon Pro";\r\n');
            fprintf(f,'\tset nodeLabels.fontSize=8;\r\n');
            fprintf(f,'\tset nodeLabels.fontStyle=0;\r\n');
            fprintf(f,'\tset nodeLabels.isShown=false;\r\n');
            fprintf(f,'\tset nodeLabels.significantDigits=4;\r\n');
            fprintf(f,'\tset nodeShape.colourAttribute="User selection";\r\n');
            fprintf(f,'\tset nodeShape.isShown=false;\r\n');
            fprintf(f,'\tset nodeShape.minSize=10.0;\r\n');
            fprintf(f,'\tset nodeShape.scaleType=Width;\r\n');
            fprintf(f,'\tset nodeShape.shapeType=Circle;\r\n');
            fprintf(f,'\tset nodeShape.size=4.0;\r\n');
            fprintf(f,'\tset nodeShape.sizeAttribute="Fixed";\r\n');
            fprintf(f,'\tset polarLayout.alignTipLabels=false;\r\n');
            fprintf(f,'\tset polarLayout.angularRange=0;\r\n');
            fprintf(f,'\tset polarLayout.rootAngle=0;\r\n');
            fprintf(f,'\tset polarLayout.rootLength=100;\r\n');
            fprintf(f,'\tset polarLayout.showRoot=false;\r\n');
            fprintf(f,'\tset radialLayout.spread=0.0;\r\n');
            fprintf(f,'\tset rectilinearLayout.alignTipLabels=false;\r\n');
            fprintf(f,'\tset rectilinearLayout.curvature=0;\r\n');
            fprintf(f,'\tset rectilinearLayout.rootLength=100;\r\n');
            fprintf(f,'\tset scale.offsetAge=0.0;\r\n');
            fprintf(f,'\tset scale.rootAge=1.0;\r\n');
            fprintf(f,'\tset scale.scaleFactor=1.0;\r\n');
            fprintf(f,'\tset scale.scaleRoot=false;\r\n');
            fprintf(f,'\tset scaleAxis.automaticScale=true;\r\n');
            fprintf(f,'\tset scaleAxis.fontSize=8.0;\r\n');
            fprintf(f,'\tset scaleAxis.isShown=false;\r\n');
            fprintf(f,'\tset scaleAxis.lineWidth=1.0;\r\n');
            fprintf(f,'\tset scaleAxis.majorTicks=0.1;\r\n');
            fprintf(f,'\tset scaleAxis.origin=0.0;\r\n');
            fprintf(f,'\tset scaleAxis.reverseAxis=false;\r\n');
            fprintf(f,'\tset scaleAxis.showGrid=true;\r\n');
            fprintf(f,'\tset scaleBar.automaticScale=true;\r\n');
            fprintf(f,'\tset scaleBar.fontSize=10.0;\r\n');
            fprintf(f,'\tset scaleBar.isShown=false;\r\n');
            fprintf(f,'\tset scaleBar.lineWidth=1.0;\r\n');
            fprintf(f,'\tset scaleBar.scaleRange=0.06;\r\n');
            fprintf(f,'\tset tipLabels.colorAttribute="User selection";\r\n');
            fprintf(f,'\tset tipLabels.displayAttribute="Names";\r\n');
            fprintf(f,'\tset tipLabels.fontName="Times New Roman";\r\n');
            fprintf(f,'\tset tipLabels.fontSize=12;\r\n');
            fprintf(f,'\tset tipLabels.fontStyle=1;\r\n');
            fprintf(f,'\tset tipLabels.isShown=false;\r\n');
            fprintf(f,'\tset tipLabels.significantDigits=4;\r\n');
            fprintf(f,'\tset trees.order=true;\r\n');
            fprintf(f,'\tset trees.orderType="increasing";\r\n');
            fprintf(f,'\tset trees.rooting=false;\r\n');
            fprintf(f,'\tset trees.rootingType="Midpoint";\r\n');
            fprintf(f,'\tset trees.transform=true;\r\n');
            fprintf(f,'\tset trees.transformType="cladogram";\r\n');
            fprintf(f,'end;\n\r\n');
            
            fclose(f);
        end
        function obj = colorgen(obj)
            obj.ColorList =  string(rgb2hex(linspecer(max(obj.NodeCommunities))));
        end
    end
    methods(Static)
        function color = newickcolorformat(hexColor)
            color = "[&!color=" +hexColor+ "]";
        end
    end
    %% View Cluster Matlab Built-in Methods
    methods(Static)
        function plotgraph(G,Com)
            %% Input Data
            Branch    = G.Edges;
            MaxWeight = max(G.Edges.Weight);
            Weights   = G.Edges.Weight;
            %% Plot the adjacency graph.
            %{
            %Layout options
            % 'Layout','force','WeightEffect','direct','Iterations',500,...
            % 'Layout','force3','WeightEffect','direct','Iterations',500,'UseGravity',true,...
            % 'Layout','layered'
            %}
            % [40,42,54]./255
            % [248 248 242]./255
            h1 = plot(G,...
                'Layout','layered',...
                'EdgeAlpha',1,...
                'EdgeColor',[40,42,54]./255,"LineWidth",2,...
                "MarkerSize",3,"NodeColor","k");
            %% Highlight Node Communities
            Cmap = get(0, 'DefaultAxesColorOrder');
            Cmap = linspecer(max(Com));
            maxC = size(Cmap,1);
            for k=1:max(Com)
                getColor = Cmap(rem(k-1,maxC)+1,:);
                Nodes = (Com==k);
                highlight(h1,Nodes,'NodeColor',getColor);hold on;
                h(k) = plot([nan,nan],[nan,nan],'o','MarkerEdgeColor',getColor,'MarkerFaceColor',getColor);
                leg(k) = string(['C_{',int2str(k),'}']);
            end
            %% Preparing Layout.
            T = title('Graph Communities','Units','normalized');
            hleg = legend(h,leg,'Location','northeastoutside','Orientation','vertical');
            htitle = get(hleg,'Title');
            set(htitle,'String','Clusters')
            axis on;
            set(hleg,'Color',[248 248 242]./255);
            set(gcf,'color','w');
            set(gca,'color',[248 248 242]./255);
            axis tight
        end
    end
end