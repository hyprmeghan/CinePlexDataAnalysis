function [ result ] = aequal( v1, v2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if abs(v1 - v2) < 0.1/6
    result = true;
else
    result = false;
end

end

