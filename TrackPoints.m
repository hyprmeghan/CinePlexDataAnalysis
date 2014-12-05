function [timeMat, thetaMat] = TrackPoints(filepath, filenames, numVars)
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
%27 November 2014

%Clear anything in the figure
clf

%Create theta matrix
thetaMat = [];
timeMat = [];

for c = 1:length(filenames)
%Make sure nothing carries over from previous runs except what we want
clearvars -except filepath filenames c thetaMat timeMat frameMat numVars
    
%Construct file name
filename = cell2mat(strcat(filepath,'/',filenames(c),'.txt'));

%opens the data file
fid = fopen(filename, 'r');

%scans the data file for comma delimited strings
data = textscan(fid, '%s', 'delimiter', ',');
data = data{1}; %something about cell arrays and formatting

%Calculates the number of variables (columns) in the data file
%based on the user input number of variables
numVars = 2 + (2*numVars);

%Calculates the length each variable will be from the length
%of the total data string based on the number of variables
Length = length(data)/numVars;

%Iniialize matrix for holding data and vector for tracking
%when a zero has been found in the data
dataMat(1:Length,numVars) = 0;

%Variable to keep track of the line we are on in data
k = 1;

%Loops over the empty datMat adding in numbers from data,
%checks for zeros, and adds the row number for zeros found
%to the list of rows to delete (foundZero)
for i = 1:Length
    for j = 1:numVars
        newDat = str2num(data{k});
              
        dataMat(i,j) = newDat;
        k = k + 1;
    end
end

%Set mins and maxes to the first of each x and y
xmin = dataMat(1,3);
xmax = dataMat(1,3);
ymin = dataMat(2,3);
ymax = dataMat(2,3);

%Find max and min for x and y

for j = 3:numVars
    cmax = max(dataMat(:,j));
    cmin = min(dataMat(:,j));
    
    if mod(j,2) == 0
        if cmax > ymax
            ymax = cmax;
        end
        if cmin < ymin
            ymin = cmin;
        end
    
    else
        if cmax > ymax
            xmax = cmax;
        end
        if cmin < ymin
            xmin = cmin;
        end
    end      
        
end

ColOrd = get(gca,'ColorOrder');
for i = 1:Length
    for j = 3:numVars - 2
        if mod(j,2) ~= 0
            hold on
            px = dataMat(i,j);
            py = dataMat(i,j + 1);

            x = dataMat(1,j + 2);
            y = dataMat(1,j + 3);

            if px == 0 || py == 0 || x == 0 || y == 0
                break
            end

            dx = px - x;
            dy = py - y;

            t1 = tan(dx/dy);
            t2 = tan(dy/dx);
            t3 = 180 - t1 - t2;

            line([px, x],[py, y], 'Color', ColOrd(j,:))
            axis([xmin,xmax,ymin,ymax])
        end

    end
    hold off
    waitforbuttonpress()
    clf
    %saveas(gcf,strcat('/Users/meghan/Desktop/Testing/Testing Photos/boop', int2str(i), '.png'))
    %i
    
end

%I = imread(strcat('/Users/meghan/Desktop/Testing/Testing Photos/boop', int2str(i),'.png'));
%imshow(I);
%implay(I);

end
