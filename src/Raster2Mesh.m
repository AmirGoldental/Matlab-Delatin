function TriMesh3D = Raster2Mesh(Raster, MaxError, ToPlot)
%TriMesh3D = Raster2Mesh(Raster, MaxError, ToPlot) creates a 3D triangulation from the
%   Raster

if ~exist('Raster', 'var')
    Raster = peaks(1000);
    Raster = Raster(1:2:end,:);
end
if ~exist('MaxError', 'var')
    MaxError = std(Raster,[],'all')/1000;
end
if ~exist('ToPlot', 'var')
    ToPlot = true;
end

TriMesh.ConnectivityList = []; % PointIdxs = TriMesh.ConnectivityList(TriangleIdx, 1:3)
TriMesh.Points = []; % PointCoordinates = TriMesh.Points(PointIdx, 1:2)
Meta.Candidates = []; % PointCoordinates = Meta.Candidates(TriangleIdx, 1:2)
Meta.Queue = PriorityQueue();
Edges = EdgesClass();

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
MaxErrors = [];
while GetMaxError() > MaxError
    StepIdx = StepIdx + 1;
    MaxErrors(end+1) = GetMaxError();
    disp(StepIdx)
    Step();
    if ToPlot
        if ~exist('Hndl', 'var')
            Hndl = figure;
            Hndl.Position = [159 121.4000 970.4000 569.6000];
            subplot(2,2,3)
            surf(Raster(1:10:end,1:10:end))
            title('Original Model, 0.5 Mil Points')
        end
        if ~ishandle(Hndl)
            break
        end
        subplot(2,2,1)
        plot(MaxErrors)
        ylabel('Max Error')
        xlabel('Step Number')
        subplot(2,2,2)
        plot(100*(1-(1:StepIdx)/numel(Raster)))
        ylabel('Compression Efficiency %')
        xlabel('Step Number')
        subplot(2,2,4)
        TriMesh3D = Create3DMesh();
        trisurf(TriMesh3D);
        title(['Mesh, ' num2str(StepIdx) ' Points'])
        drawnow;
    end
end
TriMesh3D = Create3DMesh();

   

    function TriMesh3D = Create3DMesh()
        Points = TriMesh.Points;
        Points(:,end+1) = Raster(sub2ind(size(Raster), TriMesh.Points(:,2), TriMesh.Points(:,1)));
        TriMesh3D = triangulation(TriMesh.ConnectivityList, Points);
    end

    function TriMesh2D = Create2DMesh()
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
    end

    function res = GetMaxError()
        res = -Meta.Queue.peek();
        res = res(1);
    end

% rasterize a triangle, find its max error, and queue it for processing
    function AddTriangleToQ(TriangleIdx)
        PointIdxs = TriMesh.ConnectivityList(TriangleIdx, :);
        PointCoordinates = TriMesh.Points(PointIdxs, :);
        [Xs, Ys] = meshgrid(...
            floor(min(PointCoordinates(:,1))):ceil(max(PointCoordinates(:,1))),...
            floor(min(PointCoordinates(:,2))):ceil(max(PointCoordinates(:,2))));
        B = cartesianToBarycentric(Create2DMesh(), repmat(TriangleIdx, numel(Xs), 1), [Xs(:), Ys(:)]);
        Outside = any(B<0,2);
        Xs(Outside) = [];
        Ys(Outside) = [];
        B(Outside,:) = [];
        C = barycentricToCartesian(Create3DMesh(), repmat(TriangleIdx, numel(Xs), 1), B);
        H = Raster(sub2ind(size(Raster), Ys(:), Xs(:)));
        [maxError, MaxIdx] = max(abs(H-C(:,3)));
        Y = Ys(MaxIdx);
        X = Xs(MaxIdx);
        rms = sum((H-C(:,3)).^2);
        
        % update triangle metadata
        Meta.Candidates(TriangleIdx,:) = [X,Y];
        Meta.rms(TriangleIdx) = rms;
        
        % add triangle to priority queue
        Meta.Queue.insert([-maxError, TriangleIdx, rms]);
    end


    function Step()
        % pop triangle with highest error from priority queue and split it
        T = Meta.Queue.remove();
        TriangleIdx = T(2);
        
        P1 = TriMesh.ConnectivityList(TriangleIdx, 1);
        P2 = TriMesh.ConnectivityList(TriangleIdx, 2);
        P3 = TriMesh.ConnectivityList(TriangleIdx, 3);
        
        P1_Coor = TriMesh.Points(P1, :);
        P2_Coor = TriMesh.Points(P2, :);
        P3_Coor = TriMesh.Points(P3, :);
        P_Coor = Meta.Candidates(TriangleIdx, 1:2);
        
        if IsColinear(P1_Coor, P2_Coor, P_Coor)
            Meta.Queue.remove(setdiff(Edges.GetEdge(P1, P2), TriangleIdx));
            NewTriangles = SplitEdgeToTwo(P1, P2, P_Coor);
        elseif IsColinear(P2_Coor, P3_Coor, P_Coor)
            Meta.Queue.remove(setdiff(Edges.GetEdge(P3, P2), TriangleIdx));
            NewTriangles = SplitEdgeToTwo(P3, P2, P_Coor);           
        elseif IsColinear(P3_Coor, P1_Coor, P_Coor)
            Meta.Queue.remove(setdiff(Edges.GetEdge(P1, P3), TriangleIdx));
            NewTriangles = SplitEdgeToTwo(P1, P3, P_Coor);
        else
            NewTriangles = SplitTriangleToThree(TriangleIdx, P_Coor);
        end
        for TriangleIdx = NewTriangles(:)'
            AddTriangleToQ(TriangleIdx);
        end
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
        %          P1                    P1
        %          /|\                  /  \
        %      TUL/ | \TUR          TUL/    \TUR
        %        /  |  \              /      \
        %       /   |   \    flip    /   T1   \
        %      L\ T1|T2 /R    =>    L\--------/R 
        %        \  |  /              \  T2  /
        %      TLD\ | /TRD          TLD\    /TRD
        %          \|/                  \  /
        %          P2                    P2 
        T_Idxs = Edges.GetEdge(P1, P2);
        if sum(isnan(T_Idxs))>0
            return
        end
        T1 = T_Idxs(1);
        T2 = T_Idxs(2);
        L = setdiff(TriMesh.ConnectivityList(T1,:),[P1 P2]);
        R = setdiff(TriMesh.ConnectivityList(T2,:),[P1 P2]);
        
        if IsInCircle(TriMesh.Points(L,:), TriMesh.Points(P1,:), ...
                TriMesh.Points(P2,:), TriMesh.Points(R,:))
            return
        end
        TriMesh.ConnectivityList(T1,:) = [L, P1, R];
        TriMesh.ConnectivityList(T2,:) = [L, P2, R];
        Edges.ClearEdge(P1,P2)
        
        Edges.SetEdge(L,R,[T1,T2])
               
        T_Idxs = Edges.GetEdge(L,P2);
        T_Idxs(T_Idxs==T1) = T2;
        Edges.SetEdge(L,P2,T_Idxs);
        
        T_Idxs = Edges.GetEdge(R,P1);
        T_Idxs(T_Idxs==T2) = T1;
        Edges.SetEdge(R,P1,T_Idxs);
        
        Meta.Queue.remove(T1);
        Meta.Queue.remove(T2);
        
        AddTriangleToQ(T1);
        AddTriangleToQ(T2);

        LegalizeEdge(L, P2);
        LegalizeEdge(R, P1);
    end
end

function res = IsInCircle(L, U, D, R)
LU = L-U;
LD = L-D;
RU = R-U;
RD = R-D;
Alpha = acosd(LU(:)'*LD(:)/((norm(LU)*norm(LD))));
if Alpha<0
    Alpha = Alpha + 180;
end
Gamma = acosd(RU(:)'*RD(:)/((norm(RU)*norm(RD))));
if Gamma<0
    Gamma = Gamma + 180;
end
res = Alpha+Gamma<180;
end
function res = IsColinear(a, b, c)
ab = b(:) - a(:);
ab_hat = ab./norm(ab);
ac = c(:) - a(:);
h = norm((eye(2)-ab_hat*ab_hat')*ac);
res = h<1;
end

