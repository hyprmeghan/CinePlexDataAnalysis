function [dataMat] = TrackPoints(filepath, filenames, numVars)
%CinePlex Data Analysis Function
%Parses text files from Cineplex to determine time and angle
%Plots each tracked point connected to the others with lines
%Tracked points should be in columns in the order that you would like
%those points connected in the figures

%Inputs: filepath (string) excluding file name
%        filenames (cell array of strings) file names separated by commas
%        numVars (int) number of points tracked

%Meghan Jimenez
%SINAPSE
%8 December 2014

%%%*******************************Fix then number of angles you get
%or check it at least. it might be right.

%Clear anything in the figure
clf

for c = 1:length(filenames)
    
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

%Iniialize matrix for holding data
dataMat(1:Length,numVars) = 0;

%Initialize matrix for holding theta values
thetaMat(1:Length, (numVars - 2)/2 - 2) = 0;

%Variable to keep track of the line we are on in data
k = 1;

%Loops over the empty datMat adding in numbers from data
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

%Find max and min for x and y by checking through each tracked point's 
%x and y values and comparing
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

%Plotting and display
i = 1;
ColOrd = get(gca,'ColorOrder');

%Loop over all of the variables in each from of the video
%Skips time and frame columns, looks at variables in x,y pairs and so
%skips the even (y) columns when looping
%Tracks the previous (px, py) and current (x, y) point to draw lines
%between them
lineLen = 0;
px = 0;
py = 0;
for i = 1:Length
    for j = 3:numVars - 2
        if mod(j,2) ~= 0
            ppx = px; %2-back x,y
            ppy = py;
            
            px = dataMat(i,j); %Previous x, y
            py = dataMat(i,j + 1);

            x = dataMat(i,j + 2); %Current x, y
            y = dataMat(i,j + 3);
            
            %{
            %Skip the point if anything is zero
            %This probably doesn't work as intended.
            if px == 0 || py == 0 || x == 0 || y == 0
                break
            end
            %}
            
            %Store the previous line length
            pLineLen = lineLen;
            
            %Calculate the distance between points
            dx = x - px;
            dy = y - py;
            lineLen = sqrt(dx^2 + dy^2);
            
            %Calculate the angle
            if j > 3
                hypDx = ppx - x;
                hypDy = ppy - y;
                hyp = sqrt(hypDx^2 + hypDy^2);
                theta = acosd((pLineLen^2 + lineLen^2 - hyp^2)...
                   /(2*pLineLen*lineLen));
                thetaMat(i,(j-1)/2 - 1) = theta;
            end            
        end
    end
end

i = 1;
while i <= Length
    for j = 3:numVars - 2
        if mod(j,2) ~= 0
            hold on
            
            px = dataMat(i,j); %Previous x, y
            py = dataMat(i,j + 1);

            x = dataMat(i,j + 2); %Current x, y
            y = dataMat(i,j + 3);
            
            %Plot lines
            line([px, x],[py, y], 'Color', ColOrd(j,:))
            axis([xmin,xmax,ymin,ymax])
            title(strcat('Stimulated Rat Limb Motion @ t = ', ...
                num2str( dataMat(i,2)), ', f = ', num2str(dataMat(i,1))))
            %{
                slider = uicontrol('Style','slider', 'Max', Length, ...
                'Min', 1, 'SliderStep', 1)
                %}
            if j > 3
                text(px, py + (ymax - ymin)/25, ...
                strcat(num2str(thetaMat(i,(j-1)/2 - 1)), '°'))
            end
            
            hold off

        end
        
    end
    if j > 3
        w = waitforbuttonpress();
        
        if w
            CH = get(gcf,'CurrentCharacter');
        end
        
        
        %Keyboard press navigates forward, mouse click navigates backward
        if CH == 98
            if i ~= 1
                i = i - 1;
            end
        end
        if CH == 110
            i = i + 1;
        end
        if CH == 32
            break
        end
        if CH == 102
            frame = input('Frame:');
            if frame > Length
                frame = Length;
            end
            i = frame;
        end
        
        if CH == 116
            time = input('Time (s):');
            tFrame = find(abs((dataMat(:,2) - time)) < 0.017);
            i = tFrame(1);
        end
        
        clf
        
    end
end

%figure
plot(dataMat(:,2), thetaMat)
title('Rat Leg Angular Motion')
xlabel('Time (s)')
ylabel('Angle (degrees)')

end
