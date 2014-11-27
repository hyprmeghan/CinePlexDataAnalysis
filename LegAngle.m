function [timeMat, thetaMat] = LegAngle(filepath, filenames, frameMat)
%CinePlex Data Analysis Function
%Parses text files from Cineplex to determine time and angle
%Removes any data points where any coordinate was zero
%assuming that in those cases the dots were out of the camera's
%field of view. Finds theta, checks for NAN results, determines
%beginning and end of signal by checking for stability. 
%Plots theta vs. time.

%Inputs: filepath (string) excluding file name
%        filenames (cell array of strings) file names separated by commas
%        frameMat (matrix) number of frames for each stimulation pattern
%                          each file gets a row, doesn't include patterns
%                          1 & 6

%Meghan Jimenez
%SINAPSE
%2 October 2014

%Clear anything in the figure
clf

%Create theta matrix
thetaMat = [];
timeMat = [];

for c = 1:length(filenames)
%Make sure nothing carries over from previous runs except what we want
clearvars -except filepath filenames c thetaMat timeMat frameMat
    
%Construct file name
filename = cell2mat(strcat(filepath,'/',filenames(c),'.txt'));

%opens the data file
fid = fopen(filename, 'r');

%scans the data file for comma delimited strings
data = textscan(fid, '%s', 'delimiter', ',');
data = data{1}; %something about cell arrays and formatting

%Calculates the length each variable will be from the length
%of the total data string since there are 6 variables
Length = length(data)/6;

%Iniialize matrix for holding data and vector for tracking
%when a zero has been found in the data
datMat(1:Length,6) = 0;
foundZero = [];

%Variable to keep track of the line we are on in data
k = 1;

%Loops over the empty datMat adding in numbers from data,
%checks for zeros, and adds the row number for zeros found
%to the list of rows to delete (foundZero)
for i = 1:Length
    for j = 1:6
        newDat = str2num(data{k});
        if newDat == 0 || isnan(newDat)
            foundZero(end + 1) = i;
        end        
        datMat(i,j) = newDat;
        k = k + 1;
    end
end

%datMat

%Remove repeats in foundZero
foundZero = unique(foundZero);

%Loop through datMat and remove any row that had a zero
for i = 1:length(foundZero)
    ZRow = foundZero(i) - (i-1);
    datMat = datMat([1:ZRow - 1, ZRow + 1:length(datMat)],1:6);
end

%Update Length
Length = length(datMat);

%Separate variables into vectors
time = datMat(1:length(datMat), 2);
xr = datMat(1:length(datMat), 3);
yr = datMat(1:length(datMat), 4);
xb = datMat(1:length(datMat), 5);
yb = datMat(1:length(datMat), 6);

%Find the change in x and y
dx = xb-xr;
dy = yb-yr;

%Find the angle (theta) between the two points
theta = atan(dx./dy);

%Check for NAN since we're using atan and NAN should be pi/2
for i = 1:Length
   if isnan(theta(i))
       theta(i) = pi/2;
   end
end

if (c == 2)
    plot(time, theta, '-*')
    figure
    foundZero
    find(time < 2.6)
end


%DS tracks the point where the angle destabilizes
DS = 1;

%Destabilization threshold in radians
thresh = 0.05;

%Loop through datMat to find where theta de-stabilizes
%Destabilization is determined based on the threshold and the average of
%all of the points that have been looked at already
for i = 2:Length
    diff = abs(theta(i) - mean2(theta(1:i-1)));
    if diff > thresh
        DS = i;
        break
    end 
end

DS

foundZero - DS;

for i = 1:length(foundZero)
    for j = 1:length(frameMat(c,:))
            if foundZero(i) - DS <= frameMat(c,j) && foundZero(i) - DS > 1
                frameMat(c,j) = frameMat(c,j) + 1;
            end
    end
end

frameMat;

%Remove all data before the destabilization
time = time(DS-2:Length); %Yes. the 2 is arbitrary.
theta = theta(DS-2:Length);

%Update Length
Length = length(time);

%Reset DS to cut off the end (now it stands for data stabilizes)
DS = Length;
%Reset thresh
thresh = 0.05;

%Reverse theta so we can look from the end forward
rTheta = flipud(theta);

%Loop through the reverse of theta to find where theta stabilizes
%This averages the last 10 points rather than all previous points
for i = 2:Length
    if i < 11
        diff = abs(rTheta(i) - mean2(rTheta(1:i-1)));
    else
        diff = abs(rTheta(i) - mean2(rTheta(i-10:i-1)));
    end
    if diff > thresh
        DS = i;
        break
    end 
end

%Check that the endpoint isn't longer than the length
endPoint = Length - DS + 2; %The two is still arbitrary
if endPoint > Length
    endPoint = Length;
end

%Cut off the stable stuff
time = time(1:endPoint); 
theta = theta(1:endPoint);

%Update Length
Length = length(time);

%Subtract off the difference in time so everything starts at zero
timeOff = time(1);

for i = 1:Length
    time(i) = time(i) - timeOff;
end

%Create output time and theta matrices
for i = 1:Length
    thetaMat(i, c) = theta(i);
    timeMat(i, c) = time(i);
end

if c == length(filenames)
    plot(timeMat, thetaMat, '-*')
    figure
end

end

Length = length(thetaMat);
timeMat;

m = 1;

%Loop to make sure that all times are aligned by adding zeros to ignore
%It's pretty gross but it works. If you want to make it better, find a
%work around where you keep a vector of the rows with 0's from the start
%and reconstruct using that
while m < Length
    m;
    [~, ~, val] = find(timeMat(m,:));
    minCol = find(val == min(val));

    if length(val) == 3
         if aequal(val(1), val(2)) && aequal(val(1), val(3))
             'doop';
             m = m + 1;
         else
             timeMat(end + 1, :) = [0 0 0];
             thetaMat(end + 1, :) = [0 0 0];
             Length = length(timeMat);
            if aequal(val(1), val(2))
                if val(3) > val(1)
                    'case 1a';
                    timeMat(:,3) = timeMat([1:m-1, Length, m:Length-1], 3);
                    
                    thetaMat(:,3) = thetaMat([1:m-1, Length, m:Length-1], 3);
                else
                    'case 1b';
                    timeMat(:,1) = timeMat([1:m-1, Length, m:Length-1], 1);
                    timeMat(:,2) = timeMat([1:m-1, Length, m:Length-1], 2);
                    
                    thetaMat(:,1) = thetaMat([1:m-1, Length, m:Length-1], 1);
                    thetaMat(:,2) = thetaMat([1:m-1, Length, m:Length-1], 2);
                end
            elseif aequal(val(2), val(3))
                if val(1) > val(3)
                    'case 2a';
                    timeMat(:,1) = timeMat([1:m-1, Length, m:Length-1], 1);
                    
                    thetaMat(:,1) = thetaMat([1:m-1, Length, m:Length-1], 1);
                else
                    'case 2b';
                    timeMat(:,2) = timeMat([1:m-1, Length, m:Length-1], 1);
                    timeMat(:,3) = timeMat([1:m-1, Length, m:Length-1], 2);
                    
                    thetaMat(:,2) = thetaMat([1:m-1, Length, m:Length-1], 1);
                    thetaMat(:,3) = thetaMat([1:m-1, Length, m:Length-1], 2);
                end 
            elseif aequal(val(1), val(3))
                'case 3';
                if val(2) > val(1)
                    'case 3a';
                    timeMat(:,2) = timeMat([1:m-1, Length, m:Length-1], 2);
                    
                    thetaMat(:,2) = thetaMat([1:m-1, Length, m:Length-1], 2);
                else
                    'case 3b';
                    thetaMat(:,1) = thetaMat([1:m-1, Length, m:Length-1], 1);
                    thetaMat(:,3) = thetaMat([1:m-1, Length, m:Length-1], 3);
                end
            else
                'case 4';
                thetaMat(:,minCol) = timeMat([1:m-1, Length, m:Length-1], minCol); 
            end
         end
    
    elseif length(val) == 2 && ~aequal(val(1), val(2))
        timeMat(end + 1, :) = [0 0 0];
        thetaMat(end + 1, :) = [0 0 0];
        Length = length(timeMat);
        'case 5';
        timeMat(:,minCol) = timeMat([1:m-1, Length, m:Length-1], minCol);
        thetaMat(:,minCol) = thetaMat([1:m-1, Length, m:Length-1], minCol);
        m = m + 1;
    else
        m = m + 1;
    end
    Length = length(timeMat);
end

%Gets rid of trailing zeros on the matrices
endOfZeros = Length;
rTimeMat = flipud(timeMat);

for i = 1:Length
   if rTimeMat(i,:) == [0, 0, 0]
       endOfZeros = endOfZeros - 1;
   else
       break
   end
end

timeMat = timeMat(1:endOfZeros,:);
thetaMat = thetaMat(1:endOfZeros,:);
Length = length(timeMat);

%plot(timeMat(:,1), thetaMat(:,1), '--*')

x = 0;
y = 0;

ColOrd = get(gca,'ColorOrder');

hold on
for j = 1:3
    for i = 1:6
        point = thetaMat(frameMat(j, i),j);

        if point == 0
            point = thetaMat(frameMat(j,i) + 1,j);
        end
        
        point;

        y = 46*cos(point);
        x = 46*sin(point);
        
        Col = ColOrd(j,:);

        subplot(1,6,i)
        line([50,50 + x],[50,50 - y], 'Color', Col, 'LineWidth', 2)

        axis([0,100,0,100])
        axis square
    end
end
hold off
