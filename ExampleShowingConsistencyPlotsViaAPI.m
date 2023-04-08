
function MyMain()


% In Part A1 of Project 2, at the end of a trip we plot some figures to
% inspect the acuracy and consistency of the results of the estimation
% process.
 
% We provide an small API to simplify your job, and, also, to mitigate certain
% mistakes which may confuse the marker who evaluates your solution.



file = 'DataUsr_008k.mat';
temp=load(file); data=temp.data; clear temp;
% Initialize simple API for recording expected values and covariances, so we plot error,
% and standard deviations of marginal PDFs at the end of the run.
AA=API_4010_verifyEKF(data);



P = diag([36,36,4]); 
Xe = [0;0;0];

% supose this s the events loop/


for i=1:400,
    
    % at each LiDAR event we we records current X_expected and covariance matrix.
    
    % here I just pretend these Xe and P are expected values and covariance matrix I got from EKF
    Xe=data.verify.poseL(:,i)*0.8;     
    AA.Rec(Xe,P);   % record 
    
    
end;    

% Here, later, I show the basic analisys in figure 10 (it will not make sense, because those dummy Xe and P that were recorded.;
AA.Show(10,'My Nice title');

    
% End of example using the API.
% ----------------------------------    

% note: you may implement your way of validating performance, but this one
% is necessary for marking.
