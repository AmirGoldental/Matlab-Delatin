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
this.data = Raster;
this.width = size(Raster,2);
this.height = size(Raster,1);
this.coords = []; % vertex coordinates (x, y)
this.triangles = []; % mesh triangle indices

% additional triangle data
this.halfedges = [];
this.candidates = [];
this.queueIndices = [];
this.queue = PriorityQueue();
this.pending = []; % triangles pending addition to queue
this.pendingLen = 0;

this.rmsSum = 0;

x1 = this.width - 1;
y1 = this.height - 1;
point0 = addPoint(0, 0);
point1 = addPoint(x1, 0);
point2 = addPoint(0, y1);
point3 = addPoint(x1, y1);

% add initial two triangles
Triangle0 = addTriangle(point3, point0, point2, -1, -1, -1);
addTriangle(point0, point3, point1, Triangle0, -1, -1);
flush();

% refine the mesh until its maximum error gets below the given one
StepIdx = 0;
if ToPlot
    Hndl = figure;
    Hndl.Position = [560.6000 145 804.8000 439.2000]
    subplot(1,2,1)
    surf(this.data(1:10:end,1:10:end))
    subplot(1,2,2)    
end
while getMaxError() > MaxError
    StepIdx = StepIdx + 1;
    disp(StepIdx)
    refine();
    if ToPlot
        TriMesh = create3DMesh();
        trisurf(TriMesh);
        drawnow;
    end
end
TriMesh = create3DMesh();
    function TriMesh = create3DMesh()
        Points = 1 + reshape(this.coords,2,[])';
        ConMat = 1+reshape(this.triangles,3,[])';
        Points(:,end+1) = this.data(sub2ind(size(this.data), Points(:,2), Points(:,1)));
        TriMesh = triangulation(ConMat, Points);
    end
    function TriMesh = create2DMesh()
        Points = 1 + reshape(this.coords,2,[])';
        ConMat = 1 + reshape(this.triangles,3,[])';
        TriMesh = triangulation(ConMat, Points);
    end
    function res = getMaxError()
        res = -this.queue.peek();
        res = res(1);
    end
        

% refine the mesh with a single point
    function refine()
        step();
        flush();
    end

% height value at a given position
    function res = heightAt(x, y)
        res = this.data(y+1,x+1);
    end

% rasterize and queue all triangles that got added or updated in step
    function flush()
        coords = this.coords;
        for i = 1:length(this.pending)
            t = this.pending(i);
            % rasterize triangle to find maximum pixel error
            a = 2 * this.triangles(1 + t * 3 + 0);
            b = 2 * this.triangles(1 + t * 3 + 1);
            c = 2 * this.triangles(1 + t * 3 + 2);
            findCandidate(coords(1 + a), coords(1 + a + 1),...
                coords(1 + b), coords(1 + b + 1),...
                coords(1 + c), coords(1 + c + 1), t);
        end
        this.pendingLen = 0;
    end

% rasterize a triangle, find its max error, and queue it for processing
    function findCandidateML(p0x, p0y, p1x, p1y, p2x, p2y, t)
        minX = min([p0x, p1x, p2x]);
        minY = min([p0y, p1y, p2y]);
        maxX = max([p0x, p1x, p2x]);
        maxY = max([p0y, p1y, p2y]);
        [Xs, Ys] = meshgrid(minX:maxX,minY:maxY);
        B = cartesianToBarycentric(create2DMesh(), repmat(t+1,numel(Xs),1), [Xs(:), Ys(:)]);
        C = barycentricToCartesian(create3DMesh(), repmat(t+1,numel(Xs),1), B);
        H = this.data(sub2ind(size(this.data), Ys(:)+1, Xs(:)+1));
        [maxError, MaxIdx] = max(abs(H-C(:,3)));
        [Y,X] = ind2sub(size(this.data), MaxIdx);
        Y = Y - 1;
        X = X - 1;
        rms = sum((H-C(:,3)).^2);
        % update triangle metadata
        this.candidates(1 + 2 * t) = X;
        this.candidates(1 + 2 * t + 1) = Y;
        this.rms(1 + t) = rms;
        
        % add triangle to priority queue
        this.queue.insert([-maxError, t, rms]);
    end

function findCandidate(p0x, p0y, p1x, p1y, p2x, p2y, t)
        % triangle bounding box
        minX = min([p0x, p1x, p2x]);
        minY = min([p0y, p1y, p2y]);
        maxX = max([p0x, p1x, p2x]);
        maxY = max([p0y, p1y, p2y]);
        
        % forward differencing variables
        w00 = orient(p1x, p1y, p2x, p2y, minX, minY);
        w01 = orient(p2x, p2y, p0x, p0y, minX, minY);
        w02 = orient(p0x, p0y, p1x, p1y, minX, minY);
        a01 = p1y - p0y;
        b01 = p0x - p1x;
        a12 = p2y - p1y;
        b12 = p1x - p2x;
        a20 = p0y - p2y;
        b20 = p2x - p0x;
        
        % pre-multiplied z values at vertices
        a = orient(p0x, p0y, p1x, p1y, p2x, p2y);
        z0 = heightAt(p0x, p0y) / a;
        z1 = heightAt(p1x, p1y) / a;
        z2 = heightAt(p2x, p2y) / a;
        
        % iterate over pixels in bounding box
        maxError = 0;
        mx = 0;
        my = 0;
        rms = 0;
        for y = minY:maxY
            % compute starting offset
            dx = 0;
            if (w00 < 0 && a12 ~= 0)
                dx = max([dx, floor(-w00 / a12)]);
            end
            if (w01 < 0 && a20 ~= 0)
                dx = max([dx, floor(-w01 / a20)]);
            end
            if (w02 < 0 && a01 ~= 0)
                dx = max([dx, floor(-w02 / a01)]);
            end
            
            w0 = w00 + a12 * dx;
            w1 = w01 + a20 * dx;
            w2 = w02 + a01 * dx;
            
            wasInside = false;
            
            for x = (minX + dx):maxX
                % check if inside triangle
                if (w0 >= 0 && w1 >= 0 && w2 >= 0)
                    wasInside = true;
                    
                    % compute z using barycentric coordinates
                    z = z0 * w0 + z1 * w1 + z2 * w2;
                    dz = abs(z - heightAt(x, y));
                    rms = rms + dz * dz;
                    if (dz > maxError)
                        maxError = dz;
                        mx = x;
                        my = y;
                    end
                elseif (wasInside)
                    break;
                end
                
                w0 = w0 + a12;
                w1 = w1 + a20;
                w2 = w2 + a01;
            end
            
            w00 = w00 + b12;
            w01 = w01 + b20;
            w02 = w02 + b01;
        end
        
        if (mx == p0x && my == p0y || mx == p1x && my == p1y || mx == p2x && my == p2y)
            maxError = 0;
        end
        
        % update triangle metadata
        this.candidates(1 + 2 * t) = mx;
        this.candidates(1 + 2 * t + 1) = my;
        this.rms(1 + t) = rms;
        
        % add triangle to priority queue
        this.queue.insert([-maxError, t, rms]);
    end

% process the next triangle in the queue, splitting it with a new point
    function step()
        % pop triangle with highest error from priority queue
        t = this.queue.remove();
        t = t(2);
        
        e0 = t * 3 + 0;
        e1 = t * 3 + 1;
        e2 = t * 3 + 2;
        
        p0 = this.triangles(1 + e0);
        p1 = this.triangles(1 + e1);
        p2 = this.triangles(1 + e2);
        
        ax = this.coords(1 + 2 * p0);
        ay = this.coords(1 + 2 * p0 + 1);
        bx = this.coords(1 + 2 * p1);
        by = this.coords(1 + 2 * p1 + 1);
        cx = this.coords(1 + 2 * p2);
        cy = this.coords(1 + 2 * p2 + 1);
        px = this.candidates(1 + 2 * t);
        py = this.candidates(1 + 2 * t + 1);
        
        pn = addPoint(px, py);
        
        if (orient(ax, ay, bx, by, px, py) == 0)
            handleCollinear(pn, e0);
            
        elseif (orient(bx, by, cx, cy, px, py) == 0)
            handleCollinear(pn, e1);
            
        elseif (orient(cx, cy, ax, ay, px, py) == 0)
            handleCollinear(pn, e2);
            
        else
            h0 = this.halfedges(1 + e0);
            h1 = this.halfedges(1 + e1);
            h2 = this.halfedges(1 + e2);
            
            t0 = addTriangle(p0, p1, pn, h0, -1, -1, e0);
            t1 = addTriangle(p1, p2, pn, h1, -1, t0 + 1);
            t2 = addTriangle(p2, p0, pn, h2, t0 + 2, t1 + 1);
            
            legalize(t0);
            legalize(t1);
            legalize(t2);
        end
    end

% add coordinates for a new vertex
    function i = addPoint(x, y)
        i = round(length(this.coords)/2);
        this.coords(end+(1:2)) = [x, y];
    end

% add or update a triangle in the mesh
    function e = addTriangle(a, b, c, ab, bc, ca, e)
        if ~exist('e', 'var')
            e = length(this.triangles);
        end
        t = e/3; % new triangle index
        
        % add triangle vertices
        this.triangles(1 + e + 0) = a;
        this.triangles(1 + e + 1) = b;
        this.triangles(1 + e + 2) = c;
        
        % add triangle halfedges
        this.halfedges(1 + e + 0) = ab;
        this.halfedges(1 + e + 1) = bc;
        this.halfedges(1 + e + 2) = ca;
        
        % link neighboring halfedges
        if (ab >= 0)
            this.halfedges(1 + ab) = e + 0;
        end
        if (bc >= 0)
            this.halfedges(1 + bc) = e + 1;
        end
        if (ca >= 0)
            this.halfedges(1 + ca) = e + 2;
        end
        
        % init triangle metadata
        this.candidates(1 + 2 * t + 0) = 0;
        this.candidates(1 + 2 * t + 1) = 0;
        this.queueIndices(1 + t) = -1;
        this.rms(1 + t) = 0;
        
        % add triangle to pending queue for later rasterization
        this.pending(1 + this.pendingLen) = t;
        this.pendingLen = this.pendingLen + 1;
        
    end

    function legalize(a)
        % if the pair of triangles doesn't satisfy the Delaunay condition
        % (p1 is inside the circumcircle of (p0, pl, pr)), flip them,
        % then do the same check/flip recursively for the new pair of triangles
        %
        %           pl                    pl
        %          /||\                  /  \
        %       al/ || \bl            al/    \a
        %        /  ||  \              /      \
        %       /  a||b  \    flip    /ar\
        %     p0\   ||   /p1   =>   p0\---bl---/p1
        %        \  ||  /              \      /
        %       ar\ || /br             b\    /br
        %          \||/                  \  /
        %           pr                    pr
        
        b = this.halfedges(1 + a);
        
        if (b < 0)
            return;
        end
        
        a0 = a - mod(a, 3);
        b0 = b - mod(b, 3);
        al = a0 + mod(a + 1, 3);
        ar = a0 + mod(a + 2, 3);
        bl = b0 + mod(b + 2, 3);
        br = b0 + mod(b + 1, 3);
        p0 = this.triangles(1 + ar);
        pr = this.triangles(1 + a);
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
        
        legalize(t0 + 1);
        legalize(t1 + 2);
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
            legalize(t0 + 1);
            legalize(t1 + 2);
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
        
        legalize(t0);
        legalize(t1);
        legalize(t2);
        legalize(t3);
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