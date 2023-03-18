function main()
    file = '.\DataUsr_007b.mat';
    load(file); 
    extract(data);
end

% ----------------------------------------

function extract(data)
    hh = initPlots(data);
    
    X = data.pose0; % [meters; meters; radians]
    X_buf = zeros(3, data.n, 'single');
    subsample_index = [];
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
    
    % '.\DataUsr_p020.mat' values:
    %gyroBias = -0.3984375 * pi/180; % rad / sec
    %lidar2Align = 5; % deg

    % '.\DataUsr_p021.mat' values:
    %gyroBias = 0.3984375 * pi/180; % rad / sec
    %lidar2Align = -5.703125; % deg

    for i = 1:data.n
        X_buf(:,i) = X;        
        event = events(:,i);                          

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));
       
        X_last = X;
        X = kinematicModel(X, vw, dt);    
        totalDistanceTravelled = totalDistanceTravelled + distanceBetweenCartesian(X_last, X);
        %disp(totalDistanceTravelled);
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

                [t1, localCentres] = processLiDAR(hh(1:3), X, ranges1, intensity_idx1, Lidar1Cfg, false);
                centresLidar1CF{end+1} = localCentres;
                [t2, ~] = processLiDAR(hh(4:6), X, ranges2, intensity_idx2, Lidar2Cfg, true);
                lidar1Times(end+1) = t1;
                lidar2Times(end+1) = t2;

                subsample_index(end + 1) = i;
                
                [centresGlobal, scan1Global, scan2Global] = lidarToGlobalCF(X, localCentres, ranges1, ranges2, data, lidar2Align);
                
                pairs = dataAssociation(centresGlobal, landmarks);
                prevPlotRefs = plotData(prevPlotRefs, X, centresGlobal, scan1Global, scan2Global, pairs);
                plotOOILocal(hh(7:10), ranges1, localCentres, intensity_idx1);
                %pause(0.05);
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
    plotDifferences(data, X_buf, subsample_index);
    plotLidarTimes(lidar1Times, lidar2Times, length(subsample_index));
end

% --------------------------------------------------------------------------------

function X = kinematicModel(X, vw, dt)
    dXdt = [vw(1) * cos(X(3)); vw(1) * sin(X(3)); vw(2)];
    X = X + dXdt * dt;
end

% ---------------------------------------------------------------------------------

function [t, cartesianCentres] = processLiDAR(h, X, ranges, intensity_idx, cfg, isLidar2)
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
    polarCentres = {};
    alreadyFound = [];
    for i = 1:length(intensity_idx)
        if ranges(i) == 0
            continue
        end
        if ismember(intensity_idx(i), alreadyFound)
            continue
        end
        curr_m = intensity_idx(i); curr_n = intensity_idx(i);
        next_m = curr_m;
        next_n = curr_n;
        largeSegment = false;
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
            indexRange = curr_m : 0.5 : curr_n;
            for j = 1:length(indexRange)
                if ~ismember(indexRange(j), alreadyFound)
                    alreadyFound(end+1) = indexRange(j);
                end
            end
        end
    end
    t = toc()*1000;
    plotOOI(h, polarCentres);
    %plotOOILocal(polarCentres);
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

function [centresGF, global1, global2] = lidarToGlobalCF(X, centres, ranges1, ranges2, data, lidar2Align)
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

function pairs = dataAssociation(centresGlobal, landmarks)
    pairs = {};
    thresh = 0.3; % in m
    centres = zeros(2, size(centresGlobal, 2));
    for i = 1:size(centres, 2)
        centres(:, i) = centresGlobal{i};
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
        pair = {centres(:, i), landmarks(:, minDistIdx)};
        pairs{end+1} = pair;
    end
end

% -------------------------------------------------------------------------

function plotRef = plotEstimatedX(X)
    figure(11)
    hold on;
    legend('AutoUpdate','off')
    plotRef = plot(X(1), X(2), 'b.');
end

function plotDifferences(data, X_buf, subsample_index)
    X_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(X_subsample)
        X_subsample(:,i) = X_buf(:, subsample_index(i));
    end
    
    ground = data.verify.poseL;
    x_diff = abs(ground(1,:) - X_subsample(1,:));
    y_diff = abs(ground(2,:) - X_subsample(2,:));
    
    groundPose = sqrt(ground(1,:).^2 + ground(2,:).^2);
    pose = sqrt(X_subsample(1,:).^2 + X_subsample(2,:).^2);
    pose_diff = abs(groundPose - pose);

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

    t = 1 : length(subsample_index);

    figure(12)
    subplot(211)
    plot(t, pose_diff * 100);
    title("Pose difference");
    xlabel('LiDAR event');
    ylabel('difference (cm)');

    subplot(212)
    plot(t, heading_diff);
    title("Heading difference");
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
    sprintf(direction, "\n");
end

function plotLidarTimes(times1, times2, nLidarEvents)
    figure(13); clf;
    t = 1 : nLidarEvents;
    subplot(211)
    plot(t, times1);
    title("LiDAR#1 processing time")
    xlabel('LiDAR#1 event')
    ylabel('time (ms)')
    avg1 = sum(times1) / length(times1);
    hold on;
    yline(avg1, '--', 'color', 'r');
    legend({'processing time', 'average time'})

    subplot(212)
    plot(t, times2);
    title("LiDAR#2 processing time")
    xlabel('LiDAR#2 event')
    ylabel('time (ms)')
    avg2 = sum(times2) / length(times2);
    hold on;
    yline(avg2, '--', 'color', 'r');
    legend({'processing time', 'average time'})

    avg1 = sum(times1) / length(times1);
    avg2 = sum(times2) / length(times2);

    fprintf("Average processing time LiDAR#1: %.2f ms\n", avg1);
    fprintf("Average processing time LiDAR#2: %.2f ms\n", avg2);
end

function plotOOI(h, polarCentres)
    ranges = [];
    angleIndex = [];
    %disp(polarCentres);
    for i = 1:length(polarCentres)
        idx = polarCentres{i};
        ranges(end+1) = idx(1);
        angleIndex(end+1) = indexToAngle(idx(2));
    end
    %disp(ranges);
    %disp(angleIndex);
    set(h(3), 'xdata', angleIndex, 'ydata', ranges);
end

function prev = plotData(prev, X, centresGlobal, scan1Global, scan2Global, pairs)
    figure(11);
    hold on;
    legend('AutoUpdate','off')
    for i = 1:length(prev)
        delete(prev(i));
    end
    prev = [];
    %plotRef1 = plotEstimatedX(X);
    prev(end+1) = plot(X(1), X(2), 'b.');
    
    for i = 1:length(scan1Global)
        if scan1Global(1,i) < -3 || scan1Global(1,i) > 19 || scan1Global(2,i) < -4 || scan1Global(2,i) > 18
            %scan1Global(:,i) = nan;
        end
        if scan2Global(1,i) < -3 || scan2Global(1,i) > 19 || scan2Global(2,i) < -4 || scan2Global(2,i) > 18
            %scan2Global(:,i) = nan;
        end
    end

    prev(end+1) = plot(scan1Global(1,:), scan1Global(2,:), 'c.');
    prev(end+1) = plot(scan2Global(1,:), scan2Global(2,:), 'y.');

    for i = 1:length(centresGlobal)
        centre = centresGlobal{i};
        prev(end+1) = plot(centre(1), centre(2), 'r+');
    end

    for i = 1:length(pairs)
        pair = pairs{i};
        centre = pair{1};
        landmark = pair{2};
        prev(end+1) = plot([centre(1) landmark(1)], [centre(2) landmark(2)], 'm');
    end
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

function hh=initPlots(data)
    figure(11); clf();
    
    landmarks = data.Context.Landmarks;
    plot(landmarks(1,:), landmarks(2,:), 'ko')
    title('Global CF');
    xlabel('x (m)'); 
    ylabel('y (m)');

    hold on;
    walls = data.Context.Walls;
    plot(walls(1,:), walls(2,:), 'color', [0,1,0]*0.7, 'linewidth', 3);    
    
    p0=data.pose0;
    plot(p0(1),p0(2),'r*','markersize',10);
    
    p = data.verify.poseL;
    plot(p(1,:), p(2,:), 'c.');

    plot(nan, nan,'b.');
    plot(nan, nan,'r+');
    plot(nan, nan, 'c.');
    plot(nan, nan, 'y.');
    plot(nan, nan, 'm');

    %legend({'landmarks','walls (middle planes)','initial position', 'estimated position', 'OOI centre estimate', 'LiDAR#1 scan', 'LiDAR#2 scan', 'DA links'});
    legend({'landmarks','walls (middle planes)','initial position', 'ground truth', 'estimated position', 'OOI centre estimate', 'LiDAR#1 scan', 'LiDAR#2 scan', 'DA links'});
    hold off;

    hh = initPolarPlots();
end

function h = initPolarPlots()
    figure(10); clf();
    fov = [-75:0.5:75];
    r=fov*0;
    
    subplot(211);  
    h1 = plot(fov, r, '.b');
    title('LiDAR1 polar scan');  
    xlabel('angle (degrees)');  
    ylabel('range (m)'); 
    axis([-75, 75, 0, 20]); 
    grid on;
    hold on;  
    h1b = plot(nan, nan, 'g*');
    hold on;
    h1c = plot(nan, nan, 'r+');
    legend({'opaque pixels','brilliant pixels', 'OOI centers'});
        
    
    subplot(212);  
    h2 = plot(fov, r, '.b'); 
    title('LiDAR2 polar scan');  
    xlabel('angle (degrees)');  
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
    legend({'LiDAR scanner', 'opaque pixels', 'brilliant pixels', 'OOI centres'});
    xlim([-10 10]);
    ylim([0 20]);
    h = [h1, h1b, h1c, h2, h2b, h2c, h3, h3b, h3c, h3d];
end