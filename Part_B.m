% Marcus Oates 
% z5257541
% File contains the entirety of Project1 code (all parts ABCDEF)
% usage >>Part_B('DataUsr_006k')

function main(file)
    load(file); 
    extract(data, file);
end

% ----------------------------------------

function extract(data, file)
    disp('Begin optimizing');
    gyroBias = optimize(data);
    disp('End optimizing');

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

    plot(nan, nan,'g.');
    plot(nan, nan, 'r.');
    plot(nan, nan, 'b.');

    legend({'landmarks','walls (middle planes)','initial position', 'prediction (bias=0)', 'prediction (bias)', 'ground truth'});
    hold off;

    X1 = data.pose0; % [meters; meters; radians] with estimated bias
    vw1 = [0; 0]; % [meters/sec; rad/sec] with estimated bias
    X2 = data.pose0; % [meters; meters; radians] assuming bias = 0
    vw2 = [0; 0]; % [meters/sec; rad/sec] assuming bias = 0
    events  = data.table;
    event0 = events(:,1);
    t_last = 0.0001 * double(event0(1));    
    ground = data.verify.poseL;
        
    disp('Begin sampling events'); 
    subsampleIdx = 1;

    for i = 1:data.n
        event = events(:,i);
        index = event(2);
        sensorID = event(3);

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));

        X1 = kinematicModel(X1, vw1, dt);
        X2 = kinematicModel(X2, vw2, dt);

        switch sensorID
            case 1
                figure(11); hold on; legend('AutoUpdate','off');
                plot(X2(1), X2(2), 'g.');
                plot(ground(1,subsampleIdx), ground(2,subsampleIdx), 'b.');
                plot(X1(1), X1(2), 'r.');
                hold off;
                subsampleIdx = subsampleIdx + 1;
            case 2
                vw1 = data.vw(:, index);
                vw1(2) = vw1(2) + gyroBias;
                vw2 = data.vw(:, index);
            otherwise
                continue
        end
    end

    disp('End sampling events');
end

% --------------------------------------------------------------------------------

function X = kinematicModel(X, vw, dt)
    dXdt = [vw(1) * cos(X(3)); vw(1) * sin(X(3)); vw(2)];
    X = X + dXdt * dt;
end

function bias = optimize(data)
    X = data.pose0; % [meters; meters; radians]
    X_buf = zeros(3, data.n, 'single');
    X2 = X;

    L1x = data.LidarsCfg.Lidar1.Ly;
    L1y = data.LidarsCfg.Lidar1.Lx;

    L2x = data.LidarsCfg.Lidar2.Ly;
    L2y = data.LidarsCfg.Lidar2.Lx;

    subsample_index = []; % lidar scan indexes
    vw = [0; 0]; % [meters/sec; rad/sec]
    vw2 = [0; 0];
    events  = data.table;
    event0 = events(:,1);
    t_last = 0.0001 * double(event0(1));         
        
    Lidar1Cfg=data.LidarsCfg.Lidar1;  
    Lidar2Cfg=data.LidarsCfg.Lidar2;

    prevPlotRefs = [];
    landmarks = data.Context.Landmarks;
    ground = data.verify.poseL;

    gyroBias = 0;
    lidar2Align = 0;
    biases = [];

    biasTime = 0;
    for i = 1:data.n
        X_buf(:,i) = X;
        event = events(:,i);                          

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));
        biasTime = biasTime + dt;
       
        X = kinematicModel(X, vw, dt);
        X2 = kinematicModel(X, vw, dt);

        index = event(2);
        sensorID = event(3);

        switch sensorID   
            case 1 % Lidar scan
                %fprintf('LiDAR: dt=%.1f ms @ t=%.3f\n', dt*1000, t_curr);             
                scan1 = data.scans(:, index);  
                scan2 = data.scans2(:, index);

                [ranges1, intensity_idx1] = getScan(scan1);
                [ranges2, intensity_idx2] = getScan(scan2);

                [~, localCentres, rangeToCentres, angleToCentres, alreadyFound] = processLiDAR(X, ranges1, intensity_idx1, Lidar1Cfg, false);
                [~, localCentres2, rangeToCentres2, angleToCentres2, ~] = processLiDAR(X, ranges2, intensity_idx2, Lidar2Cfg, false);
                tempLocal2 = localCentres2;
                localCentres2 = {};
                for q = 1:size(tempLocal2, 2)
                    tempPoint = tempLocal2{q};
                    localCentres2{end+1} = [tempPoint(1), tempPoint(2)];
                end

                subsample_index(end + 1) = i;
                
                [centresGlobal, scan1Global, scan2Global, alreadyFoundGlobal] = lidarToGlobalCF(X, localCentres, ranges1, ranges2, data, lidar2Align, alreadyFound);
                globalCentres2 = zeros(2, size(localCentres2, 2));
                for q = 1:size(globalCentres2, 2)
                    localPoint = localCentres2{q};
                    centreGCF = lidarToGlobal(localPoint', (X(3)-pi/2), pi, X(1:2), [L2x; L2y]);
                    globalCentres2(:, q) = [centreGCF(1); centreGCF(2)];                    
                end

                for q = 1:size(globalCentres2, 2)
                    centresGlobal{end+1} = globalCentres2(:, q);
                    localCentres{end+1} = localCentres2{q};
                    rangeToCentres(end+1) = rangeToCentres2(q);
                    angleToCentres(end+1) = angleToCentres2(q);
                end
                
                [pairs, ~, ~] = dataAssociation(centresGlobal, landmarks, localCentres, rangeToCentres, angleToCentres);

                pairsCCF = {};
                pairsGT = {};
                for q = 1:length(pairs)
                    pair = pairs{q};
                    pairGCF = pair{1};
                    pairsCCF{end+1} = inv(rotation(X(3)-pi/2))*(pairGCF-X(1:2));
                    %pairsGT{end+1} = pair{2}; 
                    pairsGT{end+1} = inv(rotation(ground(3,length(subsample_index))-pi/2))*(pairGCF-X(1:2));
                end

                options = optimset('TolFun', 0.01);
                disp('n')
                disp(X);
                lastX = X;
                [X, ~, ~] = fminsearch(@(pose) CostFunction(pose, pairsCCF, pairsGT), double(X), options);
                disp(X);
                biases(end+1) = (ground(3,length(subsample_index))-lastX(3))/biasTime;
                %disp(biases)
                biasTime = 0;
                %return
                %prevPlotRefs = plotData(prevPlotRefs, h, X, centresGlobal, scan1Global, scan2Global, pairs, alreadyFoundGlobal);
                %figure(11); hold on;

                continue;            
            case 2 % speed and gyros
                %fprintf('DR: dt=%.1f ms, v=%.2f m/s, w=%.2f deg/sec\n', dt*1000, vw.*[1;180/pi]);
                vw = data.vw(:, index);
                vw(2) = vw(2) + gyroBias;
                vw2 = data.vw(:, index);

                continue;                   
            otherwise  
                continue;
        end
    end
    bias = sum(biases) / length(biases);
    disp(bias);
end

function cost = CostFunction(X, ccfLandmarks, gtLandmarks) 
    estLandmarks = {};
    for i = 1:length(ccfLandmarks)
        localOOI = ccfLandmarks{i};
        estLandmarks{end+1} = rotation(X(3)-pi/2)*localOOI+X(1:2);
    end
    dist = [];
    for i = 1:length(estLandmarks)
        estLandmark = estLandmarks{i};
        gtLandmark = gtLandmarks{i};
        dist = sqrt((estLandmark(1)-gtLandmark(1))^2+(estLandmark(2)-gtLandmark(2))^2);
    end
    cost = sum(dist);
end

% ---------------------------------------------------------------------------------

function [t, cartesianCentres, rangeToCentres, angleToCentres, alreadyFound] = processLiDAR(X, ranges, intensity_idx, cfg, isLidar2)
    tempRanges = ranges;
    for i = 1:length(ranges)
        if tempRanges(i) == 0
            tempRanges(i) = nan;
        end
    end
    
    if isLidar2
        t = 0;
        cartesianCentres = [];
        return
    end
    
    fov = [-75:0.5:75]';
    
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
    thresh = 1; % in m
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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function plotConsistency(data, X_buf, Px_buf, X50_buf, Px50_buf, X25_buf, Px25_buf, subsample_index)
    X_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(X_subsample)
        X_subsample(:,i) = X_buf(:, subsample_index(i));
    end

    X50_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(X50_subsample)
        X50_subsample(:,i) = X50_buf(:, subsample_index(i));
    end

    X25_subsample = zeros(3, length(subsample_index), 'single');
    for i = 1:length(X25_subsample)
        X25_subsample(:,i) = X25_buf(:, subsample_index(i));
    end

    ground = data.verify.poseL;
    x_diff = abs(ground(1,:) - X_subsample(1,:));
    y_diff = abs(ground(2,:) - X_subsample(2,:));
    x50_diff = abs(ground(1,:) - X50_subsample(1,:));
    y50_diff = abs(ground(2,:) - X50_subsample(2,:));
    x25_diff = abs(ground(1,:) - X25_subsample(1,:));
    y25_diff = abs(ground(2,:) - X25_subsample(2,:));

    %groundPose = sqrt(ground(1,:).^2 + ground(2,:).^2);
    %pose = sqrt(X_subsample(1,:).^2 + X_subsample(2,:).^2);
    %pose_diff = abs(groundPose - pose);
    %pose50 = sqrt(X50_subsample(1,:).^2 + X50_subsample(2,:).^2);
    %pose_diff50 = abs(groundPose - pose50);
   % pose25 = sqrt(X25_subsample(1,:).^2 + X25_subsample(2,:).^2);
   % pose_diff25 = abs(groundPose - pose25);

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

    heading_diff50 = (X50_subsample(3,:) - ground(3,:))*180/pi;
    for i = 1:length(heading_diff50)
        if abs(heading_diff50(i)) > 180
            if ground(3,i) > X50_subsample(3,i)
                heading_diff50(i) = (2*pi - ground(3,i)) + X50_subsample(3,i);
            else
                heading_diff50(i) = ((2*pi - X50_subsample(3,i) + ground(3,i))) * -1 ;
            end
            heading_diff50(i) = heading_diff50(i)*180/pi;
        end
    end

    heading_diff25 = (X25_subsample(3,:) - ground(3,:))*180/pi;
    for i = 1:length(heading_diff25)
        if abs(heading_diff25(i)) > 180
            if ground(3,i) > X25_subsample(3,i)
                heading_diff25(i) = (2*pi - ground(3,i)) + X25_subsample(3,i);
            else
                heading_diff25(i) = ((2*pi - X25_subsample(3,i) + ground(3,i))) * -1 ;
            end
            heading_diff25(i) = heading_diff25(i)*180/pi;
        end
    end

    t = 1 : length(subsample_index);
    
    margin = zeros(3, length(subsample_index));
    for i = 1:length(margin)
        Px = Px_buf{i};
        margin(1, i) = 2*sqrt(Px(1,1));
        margin(2, i) = 2*sqrt(Px(2,2));
        margin(3, i) = 2*sqrt(Px(3,3));
    end

    margin50 = zeros(3, length(subsample_index));
    for i = 1:length(margin50)
        Px50 = Px50_buf{i};
        margin50(1, i) = 2*sqrt(Px50(1,1));
        margin50(2, i) = 2*sqrt(Px50(2,2));
        margin50(3, i) = 2*sqrt(Px50(3,3));
    end

    margin25 = zeros(3, length(subsample_index));
    for i = 1:length(margin25)
        Px25 = Px25_buf{i};
        margin25(1, i) = 2*sqrt(Px25(1,1));
        margin25(2, i) = 2*sqrt(Px25(2,2));
        margin25(3, i) = 2*sqrt(Px25(3,3));
    end

    figure(500); clf();
    subplot(3,1,1); hold on;
    plot(t, x_diff * 100);
    plot(t, margin(1,:) * 100);
    plot(t, -margin(1,:) * 100);
    title("X Difference (100% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,2); hold on;
    plot(t, y_diff * 100);
    plot(t, margin(2,:) * 100);
    plot(t, -margin(2,:) * 100);
    title("Y Difference (100% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,3); hold on;
    plot(t, heading_diff);
    plot(t, margin(3,:)*180/pi);
    plot(t, -margin(3,:)*180/pi);
    title("Heading Difference (100% Density)");
    xlabel('LiDAR event');
    ylabel('difference (deg)');
    hold off;

    figure(501); clf();
    subplot(3,1,1); hold on;
    plot(t, x50_diff * 100);
    plot(t, margin50(1,:) * 100);
    plot(t, -margin50(1,:) * 100);
    title("X Difference (50% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,2); hold on;
    plot(t, y50_diff * 100);
    plot(t, margin50(2,:) * 100);
    plot(t, -margin50(2,:) * 100);
    title("Y Difference (50% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,3); hold on;
    plot(t, heading_diff50);
    plot(t, margin50(3,:)*180/pi);
    plot(t, -margin50(3,:)*180/pi);
    title("Heading Difference (50% Density)");
    xlabel('LiDAR event');
    ylabel('difference (deg)');
    hold off;

    figure(502); clf();
    subplot(3,1,1); hold on;
    plot(t, x25_diff * 100);
    plot(t, margin25(1,:) * 100);
    plot(t, -margin25(1,:) * 100);
    title("X Difference (25% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,2); hold on;
    plot(t, y_diff * 100);
    plot(t, margin25(2,:) * 100);
    plot(t, -margin25(2,:) * 100);
    title("Y Difference (25% Density)");
    xlabel('LiDAR event');
    ylabel('difference (cm)');
    hold off;

    subplot(3,1,3); hold on;
    plot(t, heading_diff25);
    plot(t, margin25(3,:)*180/pi);
    plot(t, -margin25(3,:)*180/pi);
    title("Heading Difference (25% Density)");
    xlabel('LiDAR event');
    ylabel('difference (deg)');
    hold off;
    
end

function prev = plotData(prev, h, X, centresGlobal, scan1Global, scan2Global, pairs, alreadyFoundGlobal)
    figure(11);
    hold on;
    legend('AutoUpdate','off')
    for i = 1:length(prev)
        delete(prev(i));
    end
    prev = [];
    %set(h(1), 'xdata', X(1), 'ydata', X(2))
    %set(h(3), 'xdata', scan1Global(1,:), 'ydata', scan1Global(2,:))
    %set(h(4), 'xdata', scan2Global(1,:), 'ydata', scan2Global(2,:))

    centres = zeros(2, size(centresGlobal, 2));
    for i = 1:length(centresGlobal)
        centres(:, i) = centresGlobal{i};
    end
    %set(h(2), 'xdata', centres(1,:), 'ydata', centres(2,:))

    figure(11); hold on;
    for i = 1:length(pairs)
        pair = pairs{i};
        centre = pair{1};
        landmark = pair{2};
        prev(end+1) = plot([centre(1) landmark(1)], [centre(2) landmark(2)], 'm');
    end

    %set(h(6), 'xdata', alreadyFoundGlobal(1,:), 'ydata', alreadyFoundGlobal(2,:));
end

function [h] = initPlots(data)
    figure(11); clf();
    landmarks = data.Context.Landmarks;
    plot(landmarks(1,:), landmarks(2,:), 'ko');
    title('Global CF (100% Density)');
    xlabel('x (m)'); 
    ylabel('y (m)');

    hold on;
    walls = data.Context.Walls;
    plot(walls(1,:), walls(2,:), 'color', [0,1,0]*0.7, 'linewidth', 3);    
    
    p0=data.pose0;
    plot(p0(1),p0(2),'r*','markersize',10);

    h1a = plot(nan, nan,'g.');
    h1b = plot(nan, nan, 'r.');
    h1c = plot(nan, nan, 'b.');

    legend({'landmarks','walls (middle planes)','initial position', 'position (bias=0)', 'position (bias=est)', 'position (gt)'});
    hold off;

    h = [h1a, h1b, h1c];
end