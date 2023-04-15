% Marcus Oates 
% z5257541
% File contains the entirety of Project1 code (all parts ABCDEF)
% usage >>Part_ABCDEF('DataUsr_006b')

function main(file)
    load(file); 
    extract(data, file);
end

% ----------------------------------------

function extract(data, file)
    [h, hh] = initPlots(data);
    
    X = data.pose0; % [meters; meters; radians]
    X_buf = zeros(3, data.n, 'single');
    subsample_index = []; % lidar scan indexes
    vw = [0; 0]; % [meters/sec; rad/sec]
    events  = data.table;
    event0 = events(:,1);
    t_last = 0.0001 * double(event0(1));         
        
    Lidar1Cfg=data.LidarsCfg.Lidar1;  
    Lidar2Cfg=data.LidarsCfg.Lidar2;

    lidar1Times = [];
    lidar2Times = [];

    centresLidar1CF = {};
    centresGlobal = {};
    prevPlotRefs = [];
    landmarks = data.Context.Landmarks;
    disp('Begin sampling events');

    totalDistanceTravelled = 0;
    gyroBias = 0;
    lidar2Align = 0; 

    lastLidarPose = X(1:2);
    estX_buf = zeros(3, data.n, 'single');

    for i = 1:data.n
        X_buf(:,i) = X;        
        event = events(:,i);                          

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));
       
        X_last = X;
        X = kinematicModel(X, vw, dt);    
        totalDistanceTravelled = totalDistanceTravelled + distanceBetweenCartesian(X_last, X);
        fprintf('Total distance travelled: %.3f m\n', totalDistanceTravelled);

        index = event(2);
        sensorID = event(3);

        switch sensorID   
            case 1 % Lidar scan
                fprintf('LiDAR: dt=%.1f ms @ t=%.3f\n', dt*1000, t_curr);             
                scan1 = data.scans(:, index);  
                scan2 = data.scans2(:, index);

                [ranges1, intensity_idx1] = getScan(scan1);
                [ranges2, intensity_idx2] = getScan(scan2);

                [t1, localCentres, rangeToCentres, angleToCentres, alreadyFound] = processLiDAR(hh(1:3), X, ranges1, intensity_idx1, Lidar1Cfg, false);
                centresLidar1CF{end+1} = localCentres;
                processLiDAR(hh(4:6), X, ranges2, intensity_idx2, Lidar2Cfg, true);
                lidar1Times(end+1) = t1;

                subsample_index(end + 1) = i;
                
                [centresGlobal, scan1Global, scan2Global, alreadyFoundGlobal] = lidarToGlobalCF(X, localCentres, ranges1, ranges2, data, lidar2Align, alreadyFound);
                
                [pairs, rangeToPairs, angleToPairs] = dataAssociation(centresGlobal, landmarks, localCentres, rangeToCentres, angleToCentres);
            
                [estX, estY, estHeading, lastLidarPose] = estimatePose(hh(11), pairs, data, rangeToPairs, angleToPairs, lastLidarPose);
                estX_buf(:,i) = [estX, estY, estHeading];

                prevPlotRefs = plotData(prevPlotRefs, h, X, centresGlobal, scan1Global, scan2Global, pairs, alreadyFoundGlobal);
                plotOOILocal(hh(7:10), ranges1, localCentres, intensity_idx1);
                pause(0.05);
                continue;            
            case 2 % speed and gyros
                fprintf('DR: dt=%.1f ms, v=%.2f m/s, w=%.2f deg/sec\n', dt*1000, vw.*[1;180/pi]);
                vw = data.vw(:, index);
                vw(2) = vw(2) + gyroBias;

                continue;                   
            otherwise  
                continue;
        end
    end

    disp('End sampling events');
    plotDifferences(data, X_buf, subsample_index, estX_buf);
    plotLidarTimes(lidar1Times, length(subsample_index));
end

% --------------------------------------------------------------------------------

function X = kinematicModel(X, vw, dt)
    dXdt = [vw(1) * cos(X(3)); vw(1) * sin(X(3)); vw(2)];
    X = X + dXdt * dt;
end

% ---------------------------------------------------------------------------------

function [t, cartesianCentres, rangeToCentres, angleToCentres, alreadyFound] = processLiDAR(h, X, ranges, intensity_idx, cfg, isLidar2)
    tempRanges = ranges;
    for i = 1:length(ranges)
        if tempRanges(i) == 0
            tempRanges(i) = nan;
        end
    end
    set(h(1),'ydata', tempRanges);
    
    if isLidar2
        t = 0;
        cartesianCentres = [];
        return
    end
    
    fov = [-75:0.5:75]';
    set(h(2), 'xdata', fov(intensity_idx), 'ydata', ranges(intensity_idx));
    
    tic();
    cartesianCentres = {};
    rangeToCentres = [];
    angleToCentres = [];
    polarCentres = {};
    alreadyFound = [];
    for i = 1:length(intensity_idx)
        % invalid range reading
        if ranges(i) == 0
            continue
        end
        % previously searched
        if ismember(intensity_idx(i), alreadyFound)
            continue
        end
        curr_m = intensity_idx(i); curr_n = intensity_idx(i);
        next_m = curr_m;
        next_n = curr_n;
        largeSegment = false;
        % search for segment from point
        while 1
            if curr_m - 1 >= 1
                next_m = curr_m - 1;
            end
            if curr_n + 1 <= length(ranges)
                next_n = curr_n + 1;
            end            

            if distanceBetween(next_m, curr_n, ranges) > 0.2 || distanceBetween(next_m, curr_m, ranges) > 0.1
                next_m = curr_m;
            end
            if distanceBetween(next_n, next_m, ranges) > 0.2 || distanceBetween(next_n, curr_n, ranges) > 0.1
                next_n = curr_n;
            end

            if curr_m == next_m && curr_n == next_n
                if distanceBetween(next_m, curr_m, ranges) > 0.1 || distanceBetween(next_n, curr_n, ranges) > 0.1
                    largeSegment = true;
                end
                break
            else
                curr_m = next_m;
                curr_n = next_n;
            end
        end
         if ~largeSegment         
            rc = (ranges(curr_m) + ranges(curr_n)) / 2;
            idc = (curr_m + curr_n) / 2;
            polarCentres{end+1} = [rc, idc];
            [xc, yc] = indexRadiusToCartesian(idc, rc);
            cartesianCentres{end+1} = [xc, yc];
            rangeToCentres(end+1) = rc;
            angleToCentres(end+1) = indexToAngle(idc);
            indexRange = curr_m : 1 : curr_n;
            for j = 1:length(indexRange)
                if ~ismember(indexRange(j), alreadyFound)
                    alreadyFound(end+1) = indexRange(j);
                end
            end
        end
    end
    t = toc()*1000;
    plotOOI(h, polarCentres);
end

function d = distanceBetween(a, b, ranges)
    ra = ranges(a);
    rb = ranges(b);
    ang = angleBetween(a, b);
    d = sqrt(ra^2 + rb^2 - 2*ra*rb*cosd(ang));
end

function d = distanceBetweenCartesian(a, b)
    d = sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2);
end

function ang = angleBetweenAbs(a, b)
    ang = abs(a - b) * 0.5;
end

function ang = angleBetween(a, b)
    ang = (a - b) * 0.5;
end

function [x, y] = indexRangeToCartesian(i, ranges)
    zero = 151;
    ang = angleBetween(i, zero);
    x = ranges(i) * -sind(ang);
    y = ranges(i) * cosd(ang);
end

function [x, y] = indexRangeAlignToCartesian(i, ranges, align)
    zero = 151;
    ang = angleBetween(i, zero) + align;
    x = ranges(i) * -sind(ang);
    y = ranges(i) * cosd(ang);
end

function [x, y] = indexRadiusToCartesian(i, radius)
    zero = 151;
    ang = angleBetween(i, zero);
    x = radius * -sind(ang);
    y = radius * cosd(ang);
end

function [r, angIndex] = cartesianToPolarIndex(x, y)
    r = sqrt(x^2 + y^2);
    ang = -atand(x/y);
    angIndex = angleToIndex(ang);
end

function ang = indexToAngle(i)
    ang = (i - 151) * 0.5;
end

function i = angleToIndex(ang)
    i = round(2 * ang + 151); 
end

function [centresGF, global1, global2, alreadyFoundGlobal] = lidarToGlobalCF(X, centres, ranges1, ranges2, data, lidar2Align, alreadyFound)
    centresGF = {};
    L1x = data.LidarsCfg.Lidar1.Ly;
    L1y = data.LidarsCfg.Lidar1.Lx;
    beta1 = data.LidarsCfg.Lidar1.Alpha;

    L2x = data.LidarsCfg.Lidar2.Ly;
    L2y = data.LidarsCfg.Lidar2.Lx;
    beta2 = data.LidarsCfg.Lidar2.Alpha;
    
    alpha = X(3) - pi/2;
    for i = 1:length(centres)
        centre = centres{i};
        xl = centre(1);
        yl = centre(2);
        pg = lidarToGlobal([xl; yl], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        centresGF{end+1} = pg;
    end

    global1 = zeros(2, length(ranges1));
    global2 = zeros(2, length(ranges2));
    for i = 1:length(ranges1)
        [xl1, yl1] = indexRangeToCartesian(i, ranges1);
        [xl2, yl2] = indexRangeAlignToCartesian(i, ranges2, lidar2Align);
        pg1 = lidarToGlobal([xl1; yl1], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        pg2 = lidarToGlobal([xl2; yl2], alpha, beta2, [X(1); X(2)], [L2x; L2y]);
        global1(:,i) = [pg1(1) pg1(2)];
        global2(:,i) = [pg2(1) pg2(2)];
    end

    alreadyFoundGlobal = zeros(2, length(alreadyFound));
    for i = 1:length(alreadyFound)
        [xl, yl] = indexRangeToCartesian(alreadyFound(i), ranges1);
        pg = lidarToGlobal([xl; yl], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        alreadyFoundGlobal(:,i) = [pg(1) pg(2)];
    end
end

function r = rotation(a)
    r = [cos(a) -sin(a); sin(a) cos(a)];
end

function pg = lidarToGlobal(pl, alpha, beta, T1, T2)
    pg = rotation(alpha) * (rotation(beta) * pl + T2) + T1;
end

function [ranges, intensity_idx] = getScan(scan)
    mask1 = 16383;
    ranges = bitand(scan, mask1);
    ranges = single(ranges) * 0.01;

    mask2 = 49152;
    intensity = bitand(scan, mask2);
    intensity_idx = find(intensity > 0);
end

function [pairs, rangeToPairs, angleToPairs] = dataAssociation(centresGlobal, landmarks, localCentres, rangeToCentres, angleToCentres)
    pairs = {};
    rangeToPairs = [];
    angleToPairs = [];
    thresh = 0.3; % in m
    centres = zeros(2, size(centresGlobal, 2));
    for i = 1:size(centres, 2)
        centres(:, i) = centresGlobal{i};
    end
    centresLocal = zeros(2, size(localCentres, 2));
    for i = 1:size(localCentres, 2)
        centresLocal(:, i) = localCentres{i};
    end

    for i = 1:size(centres, 2)
        distances = [];
        for j = 1:size(landmarks, 2)
            distances(end+1) = distanceBetweenCartesian(centres(:, i), landmarks(:, j));
        end
        [minDist, minDistIdx] = min(distances);
        if minDist > thresh
            continue
        end
        pair = {centres(:, i), landmarks(:, minDistIdx), centresLocal(:, i)};
        pairs{end+1} = pair;
        rangeToPairs(end+1) = rangeToCentres(i);
        angleToPairs(end+1) = angleToCentres(i);
    end
end

function [estX, estY, estHeading, lastLidarPose] = estimatePose(h, pairs, data, rangeToPairs, angleToPairs, lastLidarPose)
    estX = nan; estY = nan; estHeading = nan;
    if size(pairs, 2) < 2
        set(h, 'xdata', nan,'ydata', nan);
        return
    end
    T2x = data.LidarsCfg.Lidar1.Ly;
    T2y = data.LidarsCfg.Lidar1.Lx;

    localCentres = zeros(2, size(pairs, 2));
    globalCentres = zeros(2, size(pairs, 2));
    for i = 1:size(pairs, 2)
        pair = pairs{i};
        localCentres(:, i) = pair{3};
        globalCentres(:, i) = pair{1};
    end
        
    estimates = zeros(3, nchoosek(size(localCentres, 2), 2));
    idx = 1;
    for i = 1:size(localCentres, 2)
        l1 = localCentres(:, i);
        g1 = globalCentres(:, i);
        x1 = g1(1);
        y1 = g1(2);
        r1 = rangeToPairs(i);
        for j = 1:size(localCentres, 2)
            if i == j
                continue
            end
            g2 = globalCentres(:, j);
            r2 = rangeToPairs(j);
            x2 = g2(1);
            y2 = g2(2);
           
            K = (x2^2-x1^2+y2^2-y1^2+r1^2-r2^2)/2;
            a = 1 + (x2-x1)^2/((y2-y1)^2);
            b = -2*x1-2*K*(x2-x1)/((y2-y1)^2)+2*y1*(x2-x1)/(y2-y1);
            c = x1^2+K^2/((y2-y1)^2)+y1^2-2*K*y1/(y2-y1)-r1^2;
            xa = (-b-sqrt(b^2-4*a*c))/(2*a);
            xb = (-b+sqrt(b^2-4*a*c))/(2*a);
            ya = K/(y2-y1)-(x2-x1)/(y2-y1)*xa;
            yb = K/(y2-y1)-(x2-x1)/(y2-y1)*xb;
            
            dista = distanceBetweenCartesian(lastLidarPose, [xa ya]);
            distb = distanceBetweenCartesian(lastLidarPose, [xb yb]);

            if dista < distb
                x = xa;
                y = ya;
            else
                x = xb;
                y = yb;
            end

            heading = atan2(y1-y,x1-x) - angleToPairs(i)*pi/180;
            pp = [1 0; 0 1] * l1 + [T2x; T2y];
            T1 = g1 - [cos(heading-pi/2) -sin(heading-pi/2); sin(heading-pi/2) cos(heading-pi/2)] * pp; 
            estimates(:, idx) = [T1(1), T1(2), heading];            
            idx = idx + 1;
        end
    end
    estX = sum(estimates(1, :)) / size(estimates, 2);
    estY = sum(estimates(2, :)) / size(estimates, 2);
    estHeading = sum(estimates(3, :)) / size(estimates, 2);
    
    lastLidarPose = [estX; estY];

    set(h, 'xdata', estX,'ydata', estY);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function plotRef = plotEstimatedX(X)
    figure(11)
    hold on;
    legend('AutoUpdate','off')
    plotRef = plot(X(1), X(2), 'b.');
end

function plotDifferences(data, X_buf, subsample_index, estX_buf)
    X_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(X_subsample)
        X_subsample(:,i) = X_buf(:, subsample_index(i));
    end
    
    estX_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(estX_subsample)
        estX_subsample(:,i) = estX_buf(:, subsample_index(i));
    end

    ground = data.verify.poseL;
    x_diff = abs(ground(1,:) - X_subsample(1,:));
    y_diff = abs(ground(2,:) - X_subsample(2,:));

    estX_diff = abs(ground(1,:) - estX_subsample(1,:));
    estY_diff = abs(ground(2,:) - estX_subsample(2,:));
    
    groundPose = sqrt(ground(1,:).^2 + ground(2,:).^2);
    pose = sqrt(X_subsample(1,:).^2 + X_subsample(2,:).^2);
    pose_diff = abs(groundPose - pose);

    estPose = sqrt(estX_subsample(1,:).^2 + estX_subsample(2,:).^2);
    estPoseDiff = abs(groundPose - estPose);

    for i = 1:length(estX_diff)
        if isnan(estX_subsample(1))
            estPoseDiff(i) = 0;
        end
    end

    heading_diff = (X_subsample(3,:) - ground(3,:))*180/pi;
    for i = 1:length(heading_diff)
        if abs(heading_diff(i)) > 180
            if ground(3,i) > X_subsample(3,i)
                heading_diff(i) = (2*pi - ground(3,i)) + X_subsample(3,i);
            else
                heading_diff(i) = ((2*pi - X_subsample(3,i) + ground(3,i))) * -1 ;
            end
            heading_diff(i) = heading_diff(i)*180/pi;
        end
    end

    estHeading_diff = (estX_subsample(3,:) - ground(3,:))*180/pi;
    for i = 1:length(estHeading_diff)
        if abs(estHeading_diff(i)) > 180
            if ground(3,i) > estX_subsample(3,i)
                estHeading_diff(i) = (2*pi - ground(3,i)) + estX_subsample(3,i);
            else
                estHeading_diff(i) = ((2*pi - estX_subsample(3,i) + ground(3,i))) * -1 ;
            end
            estHeading_diff(i) = estHeading_diff(i)*180/pi;
        end
    end

    for i = 1:length(estHeading_diff)
        if isnan(estX_subsample(1))
            estHeading_diff(i) = 0;
        end
    end

    t = 1 : length(subsample_index);

    figure(12)
    subplot(211)
    plot(t, pose_diff * 100);
    title("Pose Difference");
    xlabel('LiDAR event');
    ylabel('difference (cm)');

    subplot(212)
    plot(t, heading_diff);
    title("Heading Difference");
    xlabel('LiDAR event');
    ylabel('difference (deg)');

    figure(15)
    subplot(211)
    plot(t, estPoseDiff * 100);
    title("Estimated Pose Difference");
    xlabel('LiDAR event');
    ylabel('difference (cm)');

    subplot(212)
    plot(t, estHeading_diff);
    title("Estimated Heading Difference");
    xlabel('LiDAR event');
    ylabel('difference (deg)');
    
    maxPoseDiff = max(pose_diff) * 100;
    fprintf("Maximum pose difference: %.2f m\n", maxPoseDiff);

    maxHeadingDiff1 = max(heading_diff);
    maxHeadingDiff2 = abs(min(heading_diff));
    if maxHeadingDiff1 > maxHeadingDiff2
        maxHeadingDiff = maxHeadingDiff1;
        direction = "CCW";
    else 
        maxHeadingDiff = maxHeadingDiff2;
        direction = "CW";
    end

    fprintf("Max heading difference: %.2f deg ", maxHeadingDiff);
    fprintf(direction);
    fprintf("\n");

    avgEstPose = 0;
    count = 0;
    for i = 1:length(estPoseDiff)
        if ~isnan(estPoseDiff(i))
            avgEstPose = avgEstPose + estPoseDiff(i);
            count = count + 1;
        end
    end
    avgEstPose = avgEstPose/count;
    
    avgEstHeading = 0;
    count = 0;
    for i = 1:length(estHeading_diff)
        if ~isnan(estHeading_diff(i))
            avgEstHeading = avgEstHeading + estHeading_diff(i);
            count = count + 1;
        end
    end
    avgEstHeading = avgEstHeading/count;

    fprintf("Average estimated pose difference: %.2f m\n", avgEstPose);
    fprintf("Average estimated heading difference: %.2f deg\n", avgEstHeading);
end

function plotLidarTimes(times1, nLidarEvents)
    figure(13); clf;
    t = 1 : nLidarEvents;
    %subplot(211)
    plot(t, times1);
    title("LiDAR#1 processing time")
    xlabel('LiDAR#1 event')
    ylabel('time (ms)')
    avg1 = sum(times1) / length(times1);
    hold on;
    yline(avg1, '--', 'color', 'r');
    legend({'processing time', 'average time'})

    avg1 = sum(times1) / length(times1);

    fprintf("Average processing time LiDAR#1: %.2f ms\n", avg1);
end

function plotOOI(h, polarCentres)
    ranges = [];
    angleIndex = [];
    for i = 1:length(polarCentres)
        idx = polarCentres{i};
        ranges(end+1) = idx(1);
        angleIndex(end+1) = indexToAngle(idx(2));
    end
    set(h(3), 'xdata', angleIndex, 'ydata', ranges);
end

function prev = plotData(prev, h, X, centresGlobal, scan1Global, scan2Global, pairs, alreadyFoundGlobal)
    figure(11);
    hold on;
    legend('AutoUpdate','off')
    for i = 1:length(prev)
        delete(prev(i));
    end
    prev = [];
    set(h(1), 'xdata', X(1), 'ydata', X(2))
    set(h(3), 'xdata', scan1Global(1,:), 'ydata', scan1Global(2,:))
    set(h(4), 'xdata', scan2Global(1,:), 'ydata', scan2Global(2,:))

    centres = zeros(2, size(centresGlobal, 2));
    for i = 1:length(centresGlobal)
        centres(:, i) = centresGlobal{i};
    end
    set(h(2), 'xdata', centres(1,:), 'ydata', centres(2,:))

    figure(11); hold on;
    for i = 1:length(pairs)
        pair = pairs{i};
        centre = pair{1};
        landmark = pair{2};
        prev(end+1) = plot([centre(1) landmark(1)], [centre(2) landmark(2)], 'm');
    end

    set(h(6), 'xdata', alreadyFoundGlobal(1,:), 'ydata', alreadyFoundGlobal(2,:));
end

function plotOOILocal(h, ranges, localCentres, intensity_idx)
    cartesian = zeros(2, length(ranges));
    for i = 1:length(ranges)
        [x, y] = indexRangeToCartesian(i, ranges);
        cartesian(:, i) = [x; y];
    end
    set(h(2), 'xdata', cartesian(1,:), 'ydata', cartesian(2,:));
    intensityCartesian = zeros(2, length(intensity_idx));
    for i = 1:length(intensity_idx)
        intensityCartesian(:, i) = cartesian(:, intensity_idx(i));
    end
    set(h(3), 'xdata', intensityCartesian(1,:), 'ydata', intensityCartesian(2,:));
    
    centres = zeros(2, length(localCentres));
    for i = 1:length(localCentres)
        centres(:,i) = localCentres{i};
    end
    set(h(4), 'xdata', centres(1,:), 'ydata', centres(2,:));
end

function [h, hh] = initPlots(data)
    figure(11); clf();
    
    landmarks = data.Context.Landmarks;
    plot(landmarks(1,:), landmarks(2,:), 'ko');
    title('Global CF (Displaying Data)');
    xlabel('x (m)'); 
    ylabel('y (m)');

    hold on;
    walls = data.Context.Walls;
    plot(walls(1,:), walls(2,:), 'color', [0,1,0]*0.7, 'linewidth', 3);    
    
    p0=data.pose0;
    plot(p0(1),p0(2),'r*','markersize',10);

    h1a = plot(nan, nan,'b.');
    h1c = plot(nan, nan, 'c.');
    h1d = plot(nan, nan, 'y.');
    h1e = plot(nan, nan, 'm');
    h1f = plot(nan, nan, 'k.');
    h1b = plot(nan, nan,'r+');

    legend({'landmarks','walls (middle planes)','initial position', 'estimated position', 'LiDAR#1 scan', 'LiDAR#2 scan', 'DA links', 'OOI part', 'OOI centre estimate'});
    hold off;

    figure(32); clf();
    hold on;
    landmarks = data.Context.Landmarks;
    plot(landmarks(1,:), landmarks(2,:), 'ko')
    walls = data.Context.Walls;
    plot(walls(1,:), walls(2,:), 'color', [0,1,0]*0.7, 'linewidth', 3);
    p0=data.pose0;
    plot(p0(1),p0(2),'r*','markersize',10);
    plotRef = plot(p0(1),p0(2),'b.');
    title('Global CF (Estimated Pose)');
    xlabel('x (m)'); 
    ylabel('y (m)');
    legend({'landmarks','walls (middle planes)','initial position','estimated position'});

    hh = initPolarPlots();
    hh(end+1) = plotRef;
    h = [h1a, h1b, h1c, h1d, h1e, h1f];
end

function hh = initPolarPlots()
    figure(10); clf();
    fov = [-75:0.5:75];
    r=fov*0;
    
    subplot(211);  
    h1 = plot(fov, r, '.b');
    title('LiDAR1 Polar Scan');  
    xlabel('angle (deg)');  
    ylabel('range (m)'); 
    axis([-75, 75, 0, 20]); 
    grid on;
    hold on;  
    h1b = plot(nan, nan, 'g*');
    hold on;
    h1c = plot(nan, nan, 'r+');
    legend({'opaque pixels','intensity pixels', 'OOI centers'});
        
    
    subplot(212);  
    h2 = plot(fov, r, '.b'); 
    title('LiDAR2 Polar Scan');  
    xlabel('angle (deg)');  
    ylabel('range (m)'); 
    axis([-75, 75, 0, 20]); 
    grid on;
    hold on;  
    h2b = plot(nan, nan, 'g*');
    h2c = plot(nan, nan, 'r+');

    figure(31); clf();
    hold on;
    h3 = plot(0, 0, 'c-s');
    h3b = plot(nan, nan, 'b.');
    h3c = plot(nan, nan, 'g*');
    h3d = plot(nan, nan, 'r+');
    title('LiDAR1 Local CF');
    xlabel('x (m)');
    ylabel('y (m)');
    legend({'LiDAR scanner', 'opaque pixels', 'intensity pixels', 'OOI centres'});
    xlim([-10 10]);
    ylim([0 20]);
    hh = [h1, h1b, h1c, h2, h2b, h2c, h3, h3b, h3c, h3d];
end