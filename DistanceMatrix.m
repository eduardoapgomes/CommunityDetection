classdef DistanceMatrix
    properties
        Distances
    end
    methods
        function obj = DistanceMatrix(Dataset)
            %% 1- Distance Matrix
            % disp('distance_matrix...')
            [~,N]         = size(Dataset);
            MaxDistance   = pdist([zeros(1,N);ones(1,N)]);
            obj.Distances = squareform(pdist(Dataset)./MaxDistance);
        end
        
    end
end