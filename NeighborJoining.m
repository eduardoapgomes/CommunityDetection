classdef NeighborJoining < DistanceMatrix
    properties
        AdjacencyMatrix
        NewickFile
        LeafList
        NodeId        
    end
    properties%(Hidden,Access= 'protected')
        Tree
        isLeaf
    end
    methods
        function obj = NeighborJoining(Dataset)
            %% 1- Distance Matrix:
            obj@DistanceMatrix(Dataset);
            %% 2- SeqNeighJoin
            obj = obj.neighborjoining();
            %% 3- Extract Parameters
            obj = obj.getadjacencymatrix();
            obj = obj.getnodesidentification();
            obj = obj.getleaflist();
            obj = obj.getnewickfile();
        end
        function obj = neighborjoining(obj)
            obj.Tree = seqneighjoin(obj.Distances,'equivar');
        end
        function obj = getadjacencymatrix(obj)
            obj.AdjacencyMatrix = obj.symmetrize(getmatrix(obj.Tree));
        end
        function obj = getnodesidentification(obj)
            [~,obj.NodeId] = getmatrix(obj.Tree);
        end
        function obj = getleaflist(obj)
            obj  = obj.findleafs();
            obj  = obj.getleafs();
            obj  = obj.leaflist2array();
        end
        function obj = getnewickfile(obj)
            obj.NewickFile =  getnewickstr(obj.Tree,'BranchNames',true,'Distances',false);
        end
    end
    methods (Hidden,Access= 'protected')
        %% getleaflist Leaf List Methods
        function obj = leaflist2array(obj)
            obj.LeafList = strrep(obj.LeafList,'Leaf','');
            obj.LeafList = cellfun(@str2num,obj.LeafList);
        end
        function obj = getleafs(obj)
            obj.LeafList = obj.NodeId(obj.isLeaf);  
        end
        function obj = findleafs(obj)
            obj.isLeaf = contains(obj.NodeId,'Leaf');
        end
    end          
    methods (Static)
        function M = symmetrize(M)
            M = M+M';
        end
    end
end