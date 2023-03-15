function main()
    file = '.\DataUsr_002.mat';
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
    centres = {};

    disp('Begin sampling events');
    for i = 1:data.n
        X_buf(:,i) = X;        
        event = events(:,i);                          

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));
       
        X = kinematicModel(X, vw, dt);    

        index = event(2);
        sensorID = event(3);

        switch sensorID   
            case 1 % Lidar scan
                fprintf('LiDAR: dt=%.1f ms @ t=%.3f\n', dt*1000, t_curr);             
                scan1 = data.scans(:, index);  
                scan2 = data.scans2(:, index);            
                
                [t1, localCentres] = processLiDAR(hh(1:3), X, scan1, Lidar1Cfg);
                centresLidar1CF{end+1} = localCentres;
                [t2, ~] = processLiDAR(hh(4:6), X, scan2, Lidar2Cfg);
                lidar1Times(end+1) = t1;
                lidar2Times(end+1) = t2;

                subsample_index(end + 1) = i;
                plotEstimatedX(X);
                
                lidarToGlobalCF(X, localCentres, scan1, scan2, data);

                pause(0.05);
                continue;            
            case 2 % speed and gyros
                fprintf('DR: dt=%.1f ms, v=%.2f m/s, w=%.2f deg/sec\n', dt*1000, vw.*[1;180/pi]);
                vw = data.vw(:, index);
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

function [t, cartesianCentres] = processLiDAR(h, X, scan, cfg)
    mask1 = 16383;
    ranges = bitand(scan, mask1);
    ranges = single(ranges) * 0.01 ;
    set(h(1),'ydata', ranges);

    mask2 = 49152;
    intensity = bitand(scan, mask2);
    intensity_idx = find(intensity > 0);   
    fov = [-75:0.5:75]';
    set(h(2), 'xdata', fov(intensity_idx), 'ydata', ranges(intensity_idx));
    
    tic();
    cartesianCentres = {};
    polarCentres = {};
    alreadyFound = [];
    for i = 1:length(intensity_idx)
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
            if curr_m + 2 <= curr_n
                j = floor((curr_m + curr_n) / 2);
                [xm, ym] = indexRangeToCartesian(curr_m, ranges);
                [xn, yn] = indexRangeToCartesian(curr_n, ranges);
                [xj, yj] = indexRangeToCartesian(j, ranges);

                xc = ((xm^2+ym^2-xj^2-yj^2)*(yn-ym)-(xm^2+ym^2-xn^2-yn^2)*(yj-ym))/(2*((xn-xm)*(yj-ym)-(xj-xm)*(yn-ym)));
                yc = -1/(2*(yn-ym))*(2*xc*(xn-xm)+xm^2+ym^2-xn^2-yn^2);
                
                cartesianCentres{end+1} = [xc, yc];
                [rc, idc] = cartesianToPolarIndex(xc, yc);
                polarCentres{end+1} = [rc, idc];
            else
                rc = (ranges(curr_m) + ranges(curr_n)) / 2;
                idc = (curr_m + curr_n) / 2;
                
                polarCentres{end+1} = [rc, idc];
                [xc, yc] = indexRadiusToCartesian(idc, rc);
                cartesianCentres{end+1} = [xc, yc];
            end
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
    for i = 1:length(cartesianCentres)
        p = cartesianCentres{i};
        %cartesianCentres{i} = [-1 * p(1), p(2)];
    end
end

function d = distanceBetween(a, b, ranges)
    ra = ranges(a);
    rb = ranges(b);
    ang = angleBetween(a, b);
    d = sqrt(ra^2 + rb^2 - 2*ra*rb*cosd(ang));
end

function ang = angleBetweenAbs(a, b)
    ang = abs(a - b) * 0.5;
end

function ang = angleBetween(a, b)
    ang = (a - b) * 0.5;
    %ang = (b - a) * 0.5;
end

function [x, y] = indexRangeToCartesian(i, ranges)
    zero = 151;
    ang = angleBetween(i, zero);
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
    %if x > 0
    %    ang = -ang;
    %end
    angIndex = angleToIndex(ang);
end

function ang = indexToAngle(i)
    ang = (i - 151) * 0.5;
end

function i = angleToIndex(ang)
    i = round(2 * ang + 151); 
end

function [centresGF] = lidarToGlobalCF(X, centres, scan1, scan2, data)
    centresGF = {};
    L1x = data.LidarsCfg.Lidar1.Lx;
    L1y = data.LidarsCfg.Lidar1.Ly;
    beta1 = data.LidarsCfg.Lidar1.Alpha;

    L2x = data.LidarsCfg.Lidar2.Lx;
    L2y = data.LidarsCfg.Lidar2.Ly;
    beta2 = data.LidarsCfg.Lidar2.Alpha;
    
    mask = 16383;
    ranges1 = bitand(scan1, mask);
    ranges1 = single(ranges1) * 0.01;

    ranges2 = bitand(scan2, mask);
    ranges2 = single(ranges2) * 0.01;
    
    alpha = X(3) - pi/2;
    %disp(alpha);
    for i = 1:length(centres)
        centre = centres{i};
        %disp(centre);
        xl = centre(1);
        yl = centre(2);
        %pg = lidarToGlobal([xl; yl], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        %pp = rotation(beta1) * [xl; yl] + [L1x; L1y]; 
        %disp(pp);
        %pg = rotation(alpha) * pp + [X(1); X(2)];
        pg = lidarToGlobal([xl; yl], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        %disp(pg);
        centresGF{end+1} = pg;

        %figure(11)
        %hold on;
        %legend('AutoUpdate','off')
        %plot(pg(1), pg(2), 'r+')
    end

    global1 = zeros(2, length(ranges1));
    for i = 1:length(ranges1)
        [xl, yl] = indexRangeToCartesian(i, ranges1);
        pg = lidarToGlobal([xl; yl], alpha, beta1, [X(1); X(2)], [L1x; L1y]);
        global1(:,i) = [pg(1) pg(2)];
        %global1(:,i) = [xl yl];
    end

    figure(11)
    hold on;
    legend('AutoUpdate','off')
    plot(global1(1,:), global1(2,:), 'c.')
        
end

function r = rotation(a)
    r = [cos(a) -sin(a); sin(a) cos(a)];
end

function pg = lidarToGlobal(pl, alpha, beta, T1, T2)
    pg = rotation(alpha) * (rotation(beta) * pl + T2) + T1;
end

% -------------------------------------------------------------------------

function plotEstimatedX(X)
    figure(11)
    hold on;
    legend('AutoUpdate','off')
    plot(X(1), X(2), 'b.')
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
    plot(t, pose_diff*100);
    title("Pose difference");
    xlabel('LiDAR event');
    ylabel('difference (cm)');

    subplot(212)
    plot(t, heading_diff);
    title("Heading difference");
    xlabel('LiDAR event');
    ylabel('difference (deg)');
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
    plot(p(1,:), p(2,:), 'r.');
    
    plot(p0(1), p0(2),'b.');
    plot(nan, nan,'r+');
    legend({'landmarks','walls (middle planes)','initial position', 'ground truth', 'estimated position', 'OOI centre estimate'});
    
    hh = initPolarPlots();
end

function h = initPolarPlots()
    figure(10); clf();
    fov = [-75:0.5:75];
    r=fov*0;
    
    subplot(211);  
    h1 = plot(fov, r, '.b');
    title('LiDAR1(shown in Polar)');  
    xlabel('angle (degrees)');  
    ylabel('range (m)'); 
    axis([-75, 75, 0, 20]); 
    grid on;
    hold on;  
    h1b = plot(0, 0, 'g*');
    hold on;
    h1c = plot(0, 0, 'r+');
    legend({'opaque pixels','brilliant pixels', 'OOI centers'});
        
    
    subplot(212);  
    h2 = plot(fov, r, '.b'); 
    title('LiDAR2(shown in Polar)');  
    xlabel('angle (degrees)');  
    ylabel('range (m)'); 
    axis([-75, 75, 0, 20]); 
    grid on;
    hold on;  
    h2b = plot(0, 0, 'r+');
    h2c = plot(0, 0, 'g*');
    h = [h1, h1b, h1c, h2, h2b, h2c];
end