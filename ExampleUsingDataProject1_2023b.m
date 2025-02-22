
% Example, version 2
% Example, which shows how to read the data used in in Project1/parts A,B,C, etc.

% This example also shows a way to plot, dynamically, LiDARs' scans, etc,
% using handles of graphic objects in Matlab, 

% Read this program in detail, to understand well, how to use the data (which is critical for you). 

% You may modify this program, if you want, for implementing your solution to Project 1.


% MTRN4010.T1.2023
% If you have questions, ask the lecturer, via Moodle or email(j.guivant@unsw.edu.au)

% ---------------------------------------------------------------------------------
function Main()

% Load data, to be played back. This is just one of the may datasets you may use.
file='.\DataUsr_002.mat';   %one of the datasets we can use (simulated data, noise free, for Project 1).
% As it is "noise-free", you can debug your implementation more easily.

% load the data, now.
load(file); 

% It will load a variable named data (it is a structure)  

% Use the data.
ExploreData(data);
end

% ----------------------------------------

function ExploreData(data)
% Here, sensors' measurements (which are provided ordered, chronologically by timestamp.)
% We can easily read them in the same order in which they had occurred, chronologically, 
% when the data was collected/generated (in a real experiment o in simulation).
% Each row of the table of events contains an entry which describes a particular event (e.g. a measurement):
% --> The sampling time of the measurement, the sensor type, and an index to a record in the list of the recorded measurements of that
% particular sensor. So, we use that index to read to the specific measurement.


% You may initialize your program, before iterating through the list of events.
% You may initialize some figures and graphic objects, of which you may
% keep handles for dynamically updating them during play-back session.
hh=InitCertainPartOfMyProgram(data);


X=data.pose0;   %platform's initial pose; [x0;y0;heading0]   [meters;meters;radians]
% initial pose is necessary, for our task.

% I "copy" variables, for easy access (btw: Matlab copies vars "by reference", if used for reading)
% so that you are not wasting memory (in case you have that concern)
ne     = data.n;                % how many events?
EventsList  = data.table;            % table of events.
event0 = EventsList(:,1);            % first event.

t0=event0(1) ; t0=0.0001*double(t0); % initial time (the time of event0).
% "event" is of integer type, that is why we convert t0 to double ( i.e. to
% "real") 
% the original timestamps are expressed as unsigned integers, in which 1 unit = 0.1 millsecond.
% but in our program we express the time in seconds, using "real" (e.g.
% double precision floating point; single precision would be OK as well)


vw=[0;0];  % Program variable to keep last [speed (aka "longitudinal velocify"),heading rate] measurement.

XX=zeros(3,ne,'single');     % A buffer for recording my results (estimated poses).  size=3 x ne.
% Because I intend to record the estimated pose at the time of each event.
% You may intend to record other things, so you would define proper buffers
% for those purposes (it is your program, so you decide)


%................
% LiDARs' parameters (position and oriantation in car's frame, FoV, etc)
Lidar1Cfg=data.LidarsCfg.Lidar1;  %Info about LiDAR installation (position and orientation, ..
% .. in platform's coordinate frame.). 
Lidar2Cfg=data.LidarsCfg.Lidar2;
% these may change in other datasets, as we can install LiDARs in different
% ways (although I do not expect to do it, so that I may keep the current configuration).

% info: 'LiDAR's pose in UGV's CF=[Lx;Ly;Alpha], in [m,m,RADIANS]'
% It needs to be considered in your calculations.
%................

% Loop: read entries, one by one, for sequential processing.
for i=1:ne,
    
    XX(:,i)=X;  %record current pose; so, we can analyze the estimated trajectory after the loop ends.
    
    % get event description.
    event   = EventsList(:,i);        %event #i                    
    sensorID=event(3);                          % source (i.e. which sensor generated this event?)
    tNow=0.0001*double(event(1));               % when was this measurement taken?
    % Time in tNow is expressed in seconds.
    dt=tNow-t0;    % dt since last event ( "dt" is needed for prediction steps).
    t0=tNow ;      % remember current time, so we can calculate dt in next iteration.            
    
    % Perform prediction X(t+dt) = f(X,vw,dt) ; vw model's inputs (speed and gyroZ) 
    % apply Kinematic model here, using the last values of inputs ( in vw), 
    % and the last value of the system's state, X
        
    X=MyKinematicModel(X,vw,dt);   
    %X(t+dt) = discrete time model,  F(X(t),u(t))
    % at this line, the program variable X  contains the predicted state at the current time (in
    % variable tNow)
    
    
    index = event(2);  % where to read the actual measurement, from that sensor recorder.
        
    switch sensorID    % measurement is from which sensor?
        
        case 1  %  it is a scan from  LiDAR#1, and scan from LiDAR#2! (both, as they are synchronized)
        fprintf('Event (LiDAR), dt=[%.1fms], at t=[%.3f]\n',dt*1000,tNow);   % print some info, for debugging our code.
        
        % both LiDARs are synchronized.
        scan1 = data.scans(:,index);  
        scan2 = data.scans2(:,index);  
        
        % do something with the data.
        
        processLiDAR(hh(1:2),X,scan1,Lidar1Cfg);  % e.g. for showing scan in "global CF" ( IF you had X(t) ) 
        processLiDAR(hh(3:4),X,scan2,Lidar2Cfg);
        %etc.
        
        pause(0.05);
        continue;  %"done, next event!"
        
        %...............................
        case 2  %  It is speed encoder + gyro  (they are packed together, synchonized readings)
        vw=data.vw(:,index);    % speed and gyroZ, last updated copy.
        % I keep the variable "vw" updated to the last reading of those sensors.

        fprintf('Event (DR),dt=[%.1f ms]; v=[%.2f]m/s,w=[%.2f]deg/sec\n',dt*1000,vw.*[1;180/pi]);
        continue;  %"next!" (go to process next even, inmmediately)
        
        otherwise  % It may happen, if the dataset contains measurements from sensors 
                     %which you had not expected to process.
        %fprintf('unknown sensor, type[%d], at t=[%d]\n',sensorID, t);         
        % because I do not know what to do with this measurement, I ignore it.
        continue;
    end;
end;       % end loop, reading chronologically sensors' events.
% .....................................................

%plot( data.verify.poseL(1,:), data.verify.poseL(2,:),'m+');
disp('Loop of events ends.');

% you may, before ending your program, show some results in an "off-line fashion"
disp('Showing ground truth (your estimated trajectory should be close).)');
ShowVerification1(data);  

end
% --------------------------------------------------------------------------------

function hh=InitCertainPartOfMyProgram(data)
        
% you may initialize things you need.
% E.g.: context for some dynamic plots,etc.
    
% for local CF, in a figure.
% and for global representation, in other one (you decide how to show
% your results)

% create some figure for your stuff.
figure(11); clf();    % global CF.

% Show the map landmarks and, if it is of interest to verify your solution, the
% walls/infrastructure present there.
% (All of them are provided in Global CF)

Landmarks=data.Context.Landmarks;
% plot centres of landmarks. 
plot(Landmarks(1,:),Landmarks(2,:),'ko')
% later, during play back session, some LiDAR's pixels will appear close to some of these leandmarks. It means he LiDAR scan is
% detecting those poles. 


% plot interior of walls (they have ~20cm thickness; but the provided info just includes the ideal center of the walls
% your LiDAR scans will appear at ~ 10cm from some of those lines.    
% Wall transversal section:  :  wall left side [<--10cm-->|<--10cm-->] the other  wall side. 
hold on;
Walls = data.Context.Walls;
plot(Walls(1,:),Walls(2,:),'color',[0,1,0]*0.7,'linewidth',3);
legend({'Centers of landmarks','Walls (middle planes) '});

title('Global CF (you should show some results here)');
xlabel('X (m)'); 
ylabel('Y (m)');
p0=data.pose0;
plot(p0(1),p0(2),'r*','markersize',10);
legend({'Landmarks','Walls (middle planes)','initial position'});


hh = CreateFigureToShowScansInPolar();
%hh=[];  % array of handles you may want to use in other parts of the
%program, for updating certain graphic objects (for animations)

end

function ShowVerification1(data)

% plot some provided verification points (of platfom's pose).
% Those are the "ground truth".
% Do not expect your solution path to perfectly match those points, as those are
% the real positions, and yours are approximate ones, based on
% predictions, and using sampled inputs. 
% When using a "noise-free"  dataset, the discrepancy should be just fraction of cm, as the inputs are not
% polluted by noise, and the simulated model is the nominal analog model.
% The errors are mostly due to time discretization and sampled inputs.
% Inputs were sampled ~ @100Hz (10ms) (you can infer that from "dt".
figure(11)
hold on;
p=data.verify.poseL;
plot(p(1,:),p(2,:),'r.');
legend({'Landmarks','Walls (middle planes)','Initial pose','Ground truth (subsampled)'});
end


% ---------------------------------------------------------------------------------

function processLiDAR(h, X, scan, cfg)
    % process LiDAR scan.
    % It is your implementation, you decide what to include here.
    % You decide the input and ouput variables, etc.
    % You may introduce a delay (using pause(...) ), here; if we are showing some dynamic plots at lowe rate.
    % pause() can be used, alternatively, in the main loop (that is a better place, for that delay).
    
    
    % showing LiDAR scan, in native polar (in its Local Coordinate frame)
    % here I am just showing the ranges and also showing, in red, the brilliant pixels
    % (which are usually reflections associated to poles)
    
    % this model of LiDAR, in the configuration being use, provide 301 ranges
    % (for FoV = [-75,+75], angular resolution =0.5 degrees.)
    % each range is a 16 bits unsigned integer (class "uint16", in Matlab)
    % its first 14 bist is the actual meausured distance, in cm.
    % the highest 2 bits indicates the "intensity" of the reflection,
    % poles were painted in highly relective paint, so thet usually result
    %  in "high intensity" readings.
    
    mask1 = 16383;   %  0xFFF7 % extract lowest 14 bits (actual ranges, in cm).
    mask2 = 49152;   %  0xC000 % extract highest 2 bits ("intensity")
    ranges = bitand(scan, mask1);        % in cms, uint16
    ranges = single(ranges)*0.01 ;       % now in meters, single precision.
    % extract intensity (I call it "color")
    color = bitand(scan,mask2);
    
    
    % update some graphic object, for animation.
    set(h(1),'ydata',ranges);
    
    % what ranges are "colored"?
    ii = find(color>0);
    
    angles = [-75:0.5:75]'; 
    % refresh those "colored" ones, in the figure.
    set(h(2),'xdata',angles(ii),'ydata',ranges(ii));
end

% ---------------------------------------------------------------------------------

function h = CreateFigureToShowScansInPolar()
    
    %I create some figures in which to show some animations (e.g. LiDAR
    %scans in native polar representation)

    figure(10); clf();
    aa = [-75: 0.5 : 75];  % LiDAR's FoV  ( -75 to +75 degrees), discretized by angular resolution (0.5 deg).
    r=aa*0;   % same number of ranges  ( i.e. 301 )
    
    % create figures, and graphic objects for showing, later, data dynamically.
    % I decided to use subfigures ("subplot") (just cosmetic details, you
    % decide your style.)
    subplot(211);  h1 = plot(aa,r,'.b');
    title('LiDAR1(shown in Polar)');  xlabel('angle (degrees)');  ylabel('range (m)'); axis([-75,75,0,20]); grid on;
    hold on;  h1b = plot(0,0,'r+');
    legend({'opaque pixels','brilliant pixels'});
        
    
    subplot(212);  
    h2 = plot(aa,r,'.b'); 
    title('LiDAR2(shown in Polar)');  xlabel('angle (degrees)');  ylabel('range (m)'); axis([-75,75,0,20]); grid on;
    hold on;  h2b = plot(0,0,'r+');
    h = [h1,h1b,h2,h2b];
end    
    

function X=MyKinematicModel(X,vw,dt)
        % Here, you implement your discrete time model; e.g. using Euler's,     
        % approach, etc.
        % 
end   
% ---------------------------------------------------------------------------------
% If you have questions, ask the lecturer, via Moodle or email(
% j.guivant@unsw.edu.au)
% or ask the demonstrators.

% If you find a bug in this code, you may ask the lecturer (I will
% apreciate it, and the rest to the students, as well!)
% and if some commnent is not clear, or posted in the wrong section of the
% code, please let's us know.

% ---------------------------------------------------------------------------------
