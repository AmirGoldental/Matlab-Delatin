classdef EdgesClass < handle
    properties (Access = private)
        EdgesT1
        EdgesT2
    end
    methods
        % contstuctor of the queue
        function obj = Edges()
            obj.EdgesT1 = sparse([]);
            obj.EdgesT2 = sparse([]);
        end
        function ClearEdge(obj, Point1, Point2)
            Points = sort([Point1, Point2]);
            obj.EdgesT1(Points(1), Points(2)) = 0;
            obj.EdgesT2(Points(1), Points(2)) = 0;
        end
        function SetEdge(obj, Point1, Point2, TriangleIdxs)
            TriangleIdxs(isnan(TriangleIdxs)) = [];
            TriangleIdxs = sort(TriangleIdxs);
            if isempty(TriangleIdxs)
                obj.ClearEdge(Point1, Point2)
                return
            end
            Points = sort([Point1, Point2]);
            obj.EdgesT1(Points(1), Points(2)) = TriangleIdxs(1);
            if numel(TriangleIdxs) == 1
                obj.EdgesT2(Points(1), Points(2)) = nan;
            elseif numel(TriangleIdxs) == 2
                obj.EdgesT2(Points(1), Points(2)) = TriangleIdxs(2);
            else
                error('wtf?')
            end
        end
        function TriangleIdxs = GetEdge(obj, Point1, Point2)
            Points = sort([Point1, Point2]);
            TriangleIdxs = [obj.EdgesT1(Points(1), Points(2)), obj.EdgesT2(Points(1), Points(2))];
            TriangleIdxs(TriangleIdxs==0) = nan;
        end
        function N = numel(obj)
            N = nnz(obj.EdgesT1);
        end
        function EdgesTable = ToTable(obj)
            [P1, P2, T1] = find(obj.EdgesT1);
            [~, ~, T2] = find(obj.EdgesT2);
            EdgesTable = table(P1, P2, T1, T2);
        end
    end
end