filepath = '/Users/meghan/Desktop/CinePlexbased videos/Rat Videos and CinePlex Data/Side View Data Text Files';
filenames = {'20140723002_2','20140723009_2', '20140723010_2'};
frameMat = [1 3 7 8 9 10; 1 3 6 9 11 13; 1  3 5 9 11 21];

%frameMat = [2 4 8 10 11 19; 2 4 7 12 14 17; 2 5 7 12 14 27];
[aTime, aTheta] = LegAngle(filepath, filenames,frameMat);
aTheta(1,:)
%aTheta(10,:)
%aTheta(19,:)

'Next'
