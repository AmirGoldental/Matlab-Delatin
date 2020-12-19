function TriMesh = Raster2Mesh(Raster, XLim, YLim, MaxError)
%TriMesh = triangulation(TRI, X, Y) creates a 2D triangulation from the
%   Raster
ToPlot = 1;
if ~exist('Raster', 'var')
    Raster = peaks(1000);
    Raster = Raster(1:2:end,:);
end
if ~exist('XLim', 'var')
    XLim = [1, size(Raster,2)];
end
if ~exist('YLim', 'var')
    YLim = [1, size(Raster,1)];
end

if ~exist('MaxError', 'var')
    MaxError = 0.1;
end

this = struct();
TriMesh.ConnectivityList = []; % PointIdxs = TriMesh.ConnectivityList(TriangleIdx, 1:3)
TriMesh.Points = []; % PointCoordinates = TriMesh.Points(PointIdx, 1:2)
Meta.Candidates = []; % PointCoordinates = Meta.Candidates(TriangleIdx, 1:2)
Edges = EdgesClass();

this.coords = []; % vertex coordinates (x, y)
this.triangles = []; % mesh triangle indices

% additional triangle data
this.halfedges = [];
this.candidates = [];
this.queue = PriorityQueue();

this.rmsSum = 0;
TriMesh.Points = ...
    [1              ,1             ; ...
     size(Raster,2) ,1             ; ...
     1              ,size(Raster,1); ...
     size(Raster,2) ,size(Raster,1)];
TriMesh.ConnectivityList = [1,2,4; 1,3,4];
Edges.SetEdge(1, 4, [1,2]);
Edges.SetEdge(1, 2, [1,nan]);
Edges.SetEdge(4, 2, [1,nan]);
Edges.SetEdge(1, 3, [2,nan]);
Edges.SetEdge(3, 4, [2,nan]);

AddTriangleToQ(1);
AddTriangleToQ(2);

% refine the mesh until its maximum error gets below the given one
StepIdx = 0;
if ToPlot
    Hndl = figure;
    Hndl.Position = [560.6000 145 804.8000 439.2000];
    subplot(1,2,1)
    surf(Raster(1:10:end,1:10:end))
    subplot(1,2,2)    
end

% Test
close all
figure; DegbugPlot()
SplitEdgeToTwo(1,4,[500,250])
DegbugPlot()
SplitTriangleToThree(1,[600,100])
DegbugPlot()
SplitTriangleToThree(3,[900,200])
DegbugPlot()
SplitTriangleToThree(2,[100,200])
DegbugPlot()
SplitEdgeToTwo(1,5,[200,100])
DegbugPlot()


while getMaxError() > MaxError
    StepIdx = StepIdx + 1;
    disp(StepIdx)
    Step();
    if ToPlot
        TriMesh3D = create3DMesh();
        trisurf(TriMesh3D);
        drawnow;
    end
end
TriMesh3D = create3DMesh();

   

    function TriMesh3D = create3DMesh()
        Points = TriMesh.Points;
        ConMat = 1+reshape(this.triangles,3,[])';
        Points(:,end+1) = Raster(sub2ind(size(Raster), TriMesh.Points(:,2), TriMesh.Points(:,1)));
        TriMesh3D = triangulation(TriMesh.ConnectivityList, Points);
    end

    function TriMesh2D = create2DMesh()
        TriMesh2D = triangulation(TriMesh.ConnectivityList,TriMesh.Points);
    end
    function DegbugPlot()        
        triplot(triangulation(TriMesh.ConnectivityList,TriMesh.Points))
        for TriIdx = 1:size(TriMesh.ConnectivityList,1)
            Points = TriMesh.ConnectivityList(TriIdx,:);
            PointsCoor = TriMesh.Points(Points,:);
            text(mean(PointsCoor(:,1)), mean(PointsCoor(:,2)), ...
                 num2str(TriIdx), 'HorizontalAlignment', 'center',...
                'BackgroundColor', 'w')
        end
        for PointIdx = 1:size(TriMesh.Points,1)
            text(TriMesh.Points(PointIdx,1), TriMesh.Points(PointIdx,2),...
                ['(' num2str(PointIdx) ')'], 'HorizontalAlignment', 'center',...
                'BackgroundColor', 'w', 'FontWeight', 'bold')
        end
        EdgesTable = ToTable(Edges);
        for EdgeIdx = 1:size(EdgesTable,1)
            PointsIdxs = [EdgesTable.P1(EdgeIdx), EdgesTable.P2(EdgeIdx)];
            PointsCoor = TriMesh.Points(PointsIdxs,:);
            if ~isnan(EdgesTable.T2(EdgeIdx))
                text(mean(PointsCoor(:,1)), mean(PointsCoor(:,2)),...
                    [num2str(EdgesTable.T1(EdgeIdx)) '-' num2str(EdgesTable.T2(EdgeIdx))],...
                    'HorizontalAlignment', 'center', 'BackgroundColor', 'w')
            else
                text(mean(PointsCoor(:,1)), mean(PointsCoor(:,2)),...
                    num2str(EdgesTable.T1(EdgeIdx)),...
                    'HorizontalAlignment', 'center', 'BackgroundColor', 'w')
            end
        end
        text
    end

    function res = getMaxError()
        res = -this.queue.peek();
        res = res(1);
    end

% rasterize a triangle, find its max error, and queue it for processing
    function AddTriangleToQ(TriangleIdx)
        PointIdxs = TriMesh.ConnectivityList(TriangleIdx, :);
        PointCoordinates = TriMesh.Points(PointIdxs, :);
        [Xs, Ys] = meshgrid(...
            floor(min(PointCoordinates(:,1))):ceil(max(PointCoordinates(:,1))),...
            floor(min(PointCoordinates(:,2))):ceil(max(PointCoordinates(:,2))));
        B = cartesianToBarycentric(create2DMesh(), repmat(TriangleIdx, numel(Xs), 1), [Xs(:), Ys(:)]);
        Outside = any(B<0,2);
        Xs(Outside) = [];
        Ys(Outside) = [];
        B(Outside,:) = [];
        C = barycentricToCartesian(create3DMesh(), repmat(TriangleIdx, numel(Xs), 1), B);
        H = Raster(sub2ind(size(Raster), Ys(:), Xs(:)));
        [maxError, MaxIdx] = max(abs(H-C(:,3)));
        Y = Ys(MaxIdx);
        X = Xs(MaxIdx);
        rms = sum((H-C(:,3)).^2);
        
        % update triangle metadata
        Meta.Candidates(TriangleIdx,:) = [X,Y];
        this.rms(TriangleIdx) = rms;
        
        % add triangle to priority queue
        this.queue.insert([-maxError, TriangleIdx, rms]);
        %disp(['X: ' num2str(X) ', Y: ' num2str(Y) ' T: ' num2str(t) ', maxError: ' num2str(maxError)])
    end


    function Step()
        % pop triangle with highest error from priority queue and split it
        T = this.queue.remove();
        TriangleIdx = T(2);
        
        P1 = TriMesh.ConnectivityList(TriangleIdx, 1);
        P2 = TriMesh.ConnectivityList(TriangleIdx, 2);
        P3 = TriMesh.ConnectivityList(TriangleIdx, 3);
        
        P1_Coor = TriMesh.Points(P1, :);
        P2_Coor = TriMesh.Points(P2, :);
        P3_Coor = TriMesh.Points(P3, :);
        P_Coor = Meta.Candidates(TriangleIdx, 1:2);
        
        if IsColinear(P1_Coor, P2_Coor, P_Coor)
            EdgeIdx = GetEdgeIdx(TriangleIdx, P1_Coor, P2_Coor);
            SecondTriangleIdx = setdiff(Edges.Triangles(EdgeIdx,:), TriangleIdx);
            this.queue.remove(SecondTriangleIdx);
            NewTriangles = SplitEdgeToTwo(EdgeIdx, P_Coor);
        elseif IsColinear(P2_Coor, P3_Coor, P_Coor)
            EdgeIdx = GetEdgeIdx(TriangleIdx, P3_Coor, P2_Coor);
            SecondTriangleIdx = setdiff(Edges.Triangles(EdgeIdx,:), TriangleIdx);
            this.queue.remove(SecondTriangleIdx);
            NewTriangles = SplitEdgeToTwo(EdgeIdx, P_Coor);            
        elseif IsColinear(P3_Coor, P1_Coor, P_Coor)
            EdgeIdx = GetEdgeIdx(TriangleIdx, P3_Coor, P1_Coor);
            SecondTriangleIdx = setdiff(Edges.Triangles(EdgeIdx,:), TriangleIdx);
            this.queue.remove(SecondTriangleIdx);
            NewTriangles = SplitEdgeToTwo(EdgeIdx, P_Coor);
        else
            NewTriangles = SplitTriangleToThree(TriangleIdx, P_Coor);
        end
        for TriangleIdx = NewTriangles(:)'
            AddTriangleToQ(TriangleIdx);
        end
    end
    
    function EdgeIdx = GetEdgeIdx(TriangleIdx, PointA, PointB)
        TriangleEdges = TriMesh.TriangleEgdes(TriangleIdx, :);
        for EdgeIdxInTriangle = 1:3
            TempEdgeIdx = TriangleEdges(EdgeIdxInTriangle);
            if all(sort(Edges.PointIdxs(TempEdgeIdx,:)) == sort([PointA, PointB]))
                EdgeIdx = TempEdgeIdx;
                return
            end
        end
        error('Unexpected error')
    end

    function NewTriangles = SplitTriangleToThree(ABC, P_Coor)
        P = size(TriMesh.Points, 1) + 1;
        TriMesh.Points(P,:) = P_Coor; 
        
        ABP = ABC;
        BCP = size(TriMesh.ConnectivityList, 1) + 1;
        CAP = size(TriMesh.ConnectivityList, 1) + 2;
        A = TriMesh.ConnectivityList(ABC, 1);
        B = TriMesh.ConnectivityList(ABC, 2);
        C = TriMesh.ConnectivityList(ABC, 3);
       
        TriMesh.ConnectivityList(ABP, :) = [A, B, P];
        TriMesh.ConnectivityList(BCP, :) = [B, C, P];
        TriMesh.ConnectivityList(CAP, :) = [C, A, P];
        
        [T_Idxs] = Edges.GetEdge(A,B);
        T_Idxs(T_Idxs==ABC) = ABP;
        Edges.SetEdge(A,B,T_Idxs);
        
        [T_Idxs] = Edges.GetEdge(A,C);
        T_Idxs(T_Idxs==ABC) = CAP;
        Edges.SetEdge(A,C,T_Idxs);
        
        [T_Idxs] = Edges.GetEdge(B,C);
        T_Idxs(T_Idxs==ABC) = BCP;
        Edges.SetEdge(B,C,T_Idxs);
        
        Edges.SetEdge(A,P,[ABP, CAP]);
        Edges.SetEdge(B,P,[ABP, BCP]);
        Edges.SetEdge(C,P,[BCP, CAP]);
        
        LegalizeEdge(A,B)
        LegalizeEdge(B,C)
        LegalizeEdge(C,A)
        NewTriangles = [ABP, BCP, CAP];
    end

    function NewTriangles = SplitEdgeToTwo(A, B, NewPointCoor)
        % NewPoint = [x,y]
        %           A                     A
        %          /||\                  /||\
        %         / || \                / || \
        %        /  ||  \              /  ||  \
        %       / T1||T2 \    Split   / T1||T2 \
        %     P1\   ||   /P2   =>   P1\---C----/P2
        %        \  ||  /              \T3||T4/
        %     T1N \ || / T2N        T1N \ || / T2N
        %          \||/                  \||/
        %           B                     B
        T_Idxs = Edges.GetEdge(A, B);
        T1 = T_Idxs(1);
        T2 = T_Idxs(2);
        if ~IsColinear(TriMesh.Points(A,:),TriMesh.Points(B,:),NewPointCoor)
            error('WHAT')
        end
        % Create a new point C
        TriMesh.Points(end+1,:) = NewPointCoor;
        C = size(TriMesh.Points,1);
        
        Edges.ClearEdge(A, B);
        
        % Create the new edge AC
        Edges.SetEdge(A, C, [T1, T2]);
        
        if ~isnan(T1)
            % Find P1, T1N
            P1 = setdiff(TriMesh.ConnectivityList(T1,:), [A, B]);
            T1N = setdiff(Edges.GetEdge(B, P1), T1);
            
            % Update T1
            TriMesh.ConnectivityList(T1,:) = [P1, C, A];
            
            % Create new triangle
            T3 = size(TriMesh.ConnectivityList,1)+1;
            TriMesh.ConnectivityList(T3,:) = [P1, B, C];
            
            % Create C-P1 Edge
            Edges.SetEdge(C, P1, [T1, T3]);
            
            % Update B-P1 Edge
            Edges.SetEdge(B, P1, [T1N T3]);
        else
            T3 = nan;
        end
        if ~isnan(T2)
            % Find P2, T2N
            P2 = setdiff(TriMesh.ConnectivityList(T2,:), [A, B]);
            T2N = setdiff(Edges.GetEdge(B, P2), T2);
            
            % Update T1
            TriMesh.ConnectivityList(T2,:) = [P2, C, A];
            
            % Create new triangle
            T4 = size(TriMesh.ConnectivityList,1)+1;
            TriMesh.ConnectivityList(T4,:) = [P2, B, C];
            
            % Create C-P1 Edge
            Edges.SetEdge(C, P2, [T2, T4]);
            
            % Update B-P2 Edge
            Edges.SetEdge(B, P2, [T2N T4]);
        else
            T4 = nan;
        end
        
        % Update B-C Edge
        Edges.SetEdge(B, C, [T3, T4]);
        
        NewTriangles = [T1, T2, T3, T4];
        NewTriangles(isnan(NewTriangles)) = [];
        % Todo: legalize all edges
        if ~isnan(T1)
            LegalizeEdge(A,P1)
            LegalizeEdge(C,P1)
        end
        if ~isnan(T2)
            LegalizeEdge(A,P2)
            LegalizeEdge(C,P2)
        end
    end

    function LegalizeEdge(P1, P2)
        % if the pair of triangles doesn't satisfy the Delaunay condition
        % (p1 is inside the circumcircle of (p0, pl, pr)), flip them,
        % then do the same check/flip recursively for the new pair of triangles
        %
        %           U                     U
        %          /|\                  /  \
        %      TUL/ | \TUR          TUL/    \TUR
        %        /  |  \              /      \
        %       /   |   \    flip    /   T1   \
        %      L\ T1|T2 /R    =>    L\--------/R 
        %        \  |  /              \  T2  /
        %      TLD\ | /TRD          TLD\    /TRD
        %          \|/                  \  /
        %           D                     D 
        return
        b = this.halfedges(1 + EdgeIdx);
        
        if (b < 0)
            return;
        end
        
        a0 = EdgeIdx - mod(EdgeIdx, 3);
        b0 = b - mod(b, 3);
        al = a0 + mod(EdgeIdx + 1, 3);
        ar = a0 + mod(EdgeIdx + 2, 3);
        bl = b0 + mod(b + 2, 3);
        br = b0 + mod(b + 1, 3);
        p0 = this.triangles(1 + ar);
        pr = this.triangles(1 + EdgeIdx);
        pl = this.triangles(1 + al);
        p1 = this.triangles(1 + bl);
        coords = this.coords;
        
        if (~inCircle(...
                coords(1 + 2 * p0), coords(1 + 2 * p0 + 1),...
                coords(1 + 2 * pr), coords(1 + 2 * pr + 1),...
                coords(1 + 2 * pl), coords(1 + 2 * pl + 1),...
                coords(1 + 2 * p1), coords(1 + 2 * p1 + 1)))
            return;
        end
        
        hal = this.halfedges(1 + al);
        har = this.halfedges(1 + ar);
        hbl = this.halfedges(1 + bl);
        hbr = this.halfedges(1 + br);
        
        this.queue.remove(a0 / 3);
        this.queue.remove(b0 / 3);
        
        t0 = addTriangle(p0, p1, pl, -1, hbl, hal, a0);
        t1 = addTriangle(p1, p0, pr, t0, har, hbr, b0);
        
        LegalizeEdge(t0 + 1);
        LegalizeEdge(t1 + 2);
    end

% handle a case where new vertex is on the edge of a triangle
    function handleCollinear(pn, a)
        a0 = a - mod(a, 3);
        al = a0 + mod(a + 1, 3);
        ar = a0 + mod(a + 2, 3);
        p0 = this.triangles(1 + ar);
        pr = this.triangles(1 + a);
        pl = this.triangles(1 + al);
        hal = this.halfedges(1 + al);
        har = this.halfedges(1 + ar);
        
        b = this.halfedges(1 + a);
        
        if (b < 0)
            t0 = addTriangle(pn, p0, pr, -1, har, -1, a0);
            t1 = addTriangle(p0, pn, pl, t0, -1, hal);
            LegalizeEdge(t0 + 1);
            LegalizeEdge(t1 + 2);
            return;
        end
        
        b0 = b - mod(b, 3);
        bl = b0 + mod(b + 2, 3);
        br = b0 + mod(b + 1, 3);
        p1 = this.triangles(1 + bl);
        hbl = this.halfedges(1 + bl);
        hbr = this.halfedges(1 + br);
        
        this.queue.remove(b0 / 3);
        
        t0 = addTriangle(p0, pr, pn, har, -1, -1, a0);
        t1 = addTriangle(pr, p1, pn, hbr, -1, t0 + 1, b0);
        t2 = addTriangle(p1, pl, pn, hbl, -1, t1 + 1);
        t3 = addTriangle(pl, p0, pn, hal, t0 + 2, t2 + 1);
        
        LegalizeEdge(t0);
        LegalizeEdge(t1);
        LegalizeEdge(t2);
        LegalizeEdge(t3);
    end

    function res = inCircle(ax, ay, bx, by, cx, cy, px, py)
        dx = ax - px;
        dy = ay - py;
        ex = bx - px;
        ey = by - py;
        fx = cx - px;
        fy = cy - py;
        
        ap = dx * dx + dy * dy;
        bp = ex * ex + ey * ey;
        cp = fx * fx + fy * fy;
        
        res = dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0;
    end

    function res = orient(ax, ay, bx, by, cx, cy)
        res = (bx - cx) * (ay - cy) - (by - cy) * (ax - cx);
    end



end

function res = IsColinear(a, b, c)
ab = b(:) - a(:);
ab_hat = ab./norm(ab);
ac = c(:) - a(:);
h = norm((eye(2)-ab_hat*ab_hat')*ac);
res = h<1;
end

