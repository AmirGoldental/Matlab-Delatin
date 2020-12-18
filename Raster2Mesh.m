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
TriMesh.ConList = []; % PointIdxs = TriMesh.ConList(TriangleIdx, 1:3)
TriMesh.Points = []; % PointCoordinates = TriMesh.Points(PointIdx, 1:2)
TriMesh.TriangleEgdes = []; % EdgeIdxs = TriMesh.TriangleEgdes(TriangleIdx, 1:3) % 1.ab, 2.bc, 3.cd
Meta.Candidates = []; % PointCoordinates = Meta.Candidates(TriangleIdx, 1:2)
Edges.TriangleIdxs = []; % TriangleIdxs = Edges.TriangleIdxs(EdgeIdx, 1:2)
Edges.PointIdxs = []; % PointIdxs = Edges.PointIdxs(EdgeIdx, 1:2)

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
TriMesh.ConList = [1,2,4; 1,3,4];
TriMesh.TriangleEgdes = [1,2,3; 1,4,5];
Edges.TriangleIdxs = [1, 2;nan(4,2)];
Edges.PointIdxs = [...
    1, 4; ...
    1, 2; ...
    4, 2; ...
    1, 3; ...
    3, 4];

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
while getMaxError() > MaxError
    StepIdx = StepIdx + 1;
    disp(StepIdx)
    Step();
    if ToPlot
        TriMesh = create3DMesh();
        trisurf(TriMesh);
        drawnow;
    end
end
TriMesh = create3DMesh();
    function TriMesh3D = create3DMesh()
        Points = TriMesh.Points;
        ConMat = 1+reshape(this.triangles,3,[])';
        Points(:,end+1) = Raster(sub2ind(size(Raster), TriMesh.Points(:,2), TriMesh.Points(:,1)));
        TriMesh3D = triangulation(TriMesh.ConList, Points);
    end
    function TriMesh2D = create2DMesh()
        TriMesh2D = triangulation(TriMesh.ConList, TriMesh.Points);
    end
    function res = getMaxError()
        res = -this.queue.peek();
        res = res(1);
    end

% rasterize a triangle, find its max error, and queue it for processing
    function AddTriangleToQ(TriangleIdx)
        PointIdxs = TriMesh.ConList(TriangleIdx, :);
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
        
        P1_Idx = TriMesh.ConList(TriangleIdx, 1);
        P2_Idx = TriMesh.ConList(TriangleIdx, 2);
        P3_Idx = TriMesh.ConList(TriangleIdx, 3);
        
        P1_Coor = TriMesh.Points(P1_Idx, :);
        P2_Coor = TriMesh.Points(P2_Idx, :);
        P3_Coor = TriMesh.Points(P3_Idx, :);
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
    
    function EdgeIdx = GetEdgeIdx(TriangleIdx, PointA_Idx, PointB_Idx)
        TriangleEdges = TriMesh.TriangleEgdes(TriangleIdx, :);
        for EdgeIdxInTriangle = 1:3
            TempEdgeIdx = TriangleEdges(EdgeIdxInTriangle);
            if all(sort(Edges.PointIdxs(TempEdgeIdx,:)) == sort([PointA_Idx, PointB_Idx]))
                EdgeIdx = TempEdgeIdx;
                return
            end
        end
        error('Unexpected error')
    end

    function NewTriangles = SplitTriangleToThree(ABC_Idx, P_Coor)
        P_Idx = size(TriMesh.Points, 1) + 1;
        TriMesh.Points(P_Idx,:) = P_Coor; 
        
        ABP_Idx = ABC_Idx;
        BCP_Idx = size(TriMesh.ConList, 1) + 1;
        CAP_Idx = size(TriMesh.ConList, 1) + 2;
        A_Idx = TriMesh.ConList(ABC_Idx, 1);
        B_Idx = TriMesh.ConList(ABC_Idx, 2);
        C_Idx = TriMesh.ConList(ABC_Idx, 3);
       
        AB_Idx = GetEdgeIdx(ABC_Idx, A_Idx, B_Idx);
        BC_Idx = GetEdgeIdx(ABC_Idx, B_Idx, C_Idx);
        CA_Idx = GetEdgeIdx(ABC_Idx, A_Idx, C_Idx);

        TriMesh.ConList(ABP_Idx, 3) = P_Idx;
        TriMesh.ConList(BCP_Idx, :) = [B_Idx, C_Idx, P_Idx];
        TriMesh.ConList(CAP_Idx, :) = [C_Idx, A_Idx, P_Idx];

        Edges.TriangleIdxs(AB_Idx, Edges.TriangleIdxs(AB_Idx,:)==ABC_Idx) = ABP_Idx;
        Edges.TriangleIdxs(BC_Idx, Edges.TriangleIdxs(BC_Idx,:)==ABC_Idx) = BCP_Idx;
        Edges.TriangleIdxs(CA_Idx, Edges.TriangleIdxs(CA_Idx,:)==ABC_Idx) = CAP_Idx;
        
        EdgeAP = size(Edges.PointIdxs,1) + 1;
        EdgeBP = size(Edges.PointIdxs,1) + 2;
        EdgeCP = size(Edges.PointIdxs,1) + 3;
        Edges.TriangleIdxs(EdgeAP, :) = [ABP_Idx, ABC_Idx];
        Edges.TriangleIdxs(EdgeBP, :) = [BCP_Idx, ABC_Idx];
        Edges.TriangleIdxs(EdgeCP, :) = [CAP_Idx, ABC_Idx];
        Edges.PointIdxs(EdgeAP, :) = [A_Idx, P_Idx];
        Edges.PointIdxs(EdgeBP, :) = [B_Idx, P_Idx];
        Edges.PointIdxs(EdgeCP, :) = [C_Idx, P_Idx];
        LegalizeEdge(AB_Idx)
        LegalizeEdge(BC_Idx)
        LegalizeEdge(CA_Idx)
        NewTriangles = [ABP_Idx, BCP_Idx, CAP_Idx];
    end

    function NewTriangles = SplitEdgeToTwo(EgdeIdx, NewPoint)
        % NewPoint = [x,y]
        %           A                     A
        %          /||\                  /||\
        %         / || \                / || \
        %        /  ||  \              /  ||  \
        %       / T1||T2 \    Split   / T1||T2 \
        %     P1\   ||   /P2   =>   P1\---B----/P2
        %        \  ||  /              \T3||T4/
        %     T1N \ || / T2N        T1N \ || / T2N
        %          \||/                  \||/
        %           B                     C
        
        AB_Idx = EgdeIdx;
        T1_Idx = Edges.TriangleIdxs(AB_Idx, 1);
        T2_Idx = Edges.TriangleIdxs(AB_Idx, 2);
        A_Idx = Edges.PointIdxs(AB_Idx, 1);
        B_Idx = Edges.PointIdxs(AB_Idx, 2);
        
        % Create a new point (C) where B is now
        TriMesh.Points(end+1,:) = TriMesh.Points(B_Idx,:);
        C_Idx = size(TriMesh.Points,1);
        
        % Move B
        TriMesh.Points(B_Idx,:) = NewPoint;
        
        % Create the new edge BC
        Edges.PointIdxs(end+1,:) = [B_Idx, C_Idx];
        BC_Idx = size(Edges.PointIdxs,1);
        
        if ~isnan(T1_Idx)
            % Find P1, B-P1 Edge, T1N
            P1_Idx = setdiff(TriMesh.ConList(T1_Idx,:), [A_Idx, B_Idx]);
            AP1_Idx = GetEdgeIdx(T1, A_Idx, P1_Idx);
            BP1_Idx = GetEdgeIdx(T1, B_Idx, P1_Idx);
            T1N_Idx = setdiff(Edges.TriangleIdxs(BP1_Idx,:), T1);
            
            % Create new triangle
            TriMesh.ConList(end+1,:) = [P1_Idx, B_Idx, C_Idx];
            T3_Idx = size(TriMesh.ConList,1);
           
            % Update B-P1 Edge
            Edges.TriangleIdxs(BP1_Idx, Edges.TriangleIdxs(BP1_Idx,:)==T1N_Idx) = T3_Idx;
            
            % Create C-P1 Edge
            Edges.TriangleIdxs(end+1, :) = [T3_Idx, T1N_Idx];
            CP1_Idx = size(Edges.TriangleIdxs,1);
            Edges.Points(CP1_Idx, :) = [C_Idx, P1_Idx];
            
            % Update TriangleEgdes
            TriMesh.TriangleEgdes(T1N_Idx, TriMesh.TriangleEgdes(T1N_Idx, :)==BP1_Idx) = CP1_Idx;
            TriMesh.TriangleEgdes(T3_Idx, :) = [CP1_Idx, BP1_Idx, BC_Idx];
        else
            T3_Idx = nan;
        end
        if ~isnan(T2_Idx)
            % Find P2, B-P2 Edge, T2N
            P2_Idx = setdiff(TriMesh.ConList(T2_Idx,:), [A_Idx, B_Idx]);
            AP2_Idx = GetEdgeIdx(T2, A_Idx, P2_Idx);
            BP2_Idx = GetEdgeIdx(T2, B_Idx, P2_Idx);
            T2N_Idx = setdiff(Edges.TriangleIdxs(BP2_Idx,:), T1);
            
            % Create new triangle
            TriMesh.ConList(end+1,:) = [P2_Idx, B_Idx, C_Idx];
            T4_Idx = size(TriMesh.ConList,1);
           
            % Update B-P2 Edge
            Edges.TriangleIdxs(BP2_Idx, Edges.TriangleIdxs(BP2_Idx,:)==T2N_Idx) = T4_Idx;
            
            % Create C-P2 Edge
            Edges.TriangleIdxs(end+1, :) = [T4_Idx, T2N_Idx];
            CP2_Idx = size(Edges.TriangleIdxs,1);
            Edges.Points(CP2_Idx, :) = [C_Idx, P2_Idx];
            
            % Update TriangleEgdes
            TriMesh.TriangleEgdes(T2N_Idx, TriMesh.TriangleEgdes(T2N_Idx, :)==BP2_Idx) = CP2_Idx;
            TriMesh.TriangleEgdes(T4_Idx, :) = [CP2_Idx, BP2_Idx, BC_Idx];
        else
            T4_Idx = nan;
        end
        
        % Update B-C Edge
        Edges.TriangleIdxs(BC_Idx, :) = [T3_Idx, T4_Idx];
        NewTriangles = [T1_Idx, T2_Idx, T3_Idx, T4_Idx];
        NewTriangles(isnan(NewTriangles)) = [];
        % Todo: legalize all edges
        if ~isnan(T1_Idx)
            LegalizeEdge(AP1_Idx)
            LegalizeEdge(CP1_Idx)
        end
        if ~isnan(T2_Idx)
            LegalizeEdge(AP2_Idx)
            LegalizeEdge(CP2_Idx)
        end
    end

    function LegalizeEdge(EdgeIdx)
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
bc = b(:) - c(:);
ac = a(:) - c(:);
res = all(cross([bc;0], [ac;0]) == 0);
end