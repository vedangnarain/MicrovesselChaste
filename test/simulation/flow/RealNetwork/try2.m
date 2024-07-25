clear all
clc
close all

folder = 'home/jakub/Desktop';  % You specify this!
fullMatFileName = fullfile(folder,  'LargestComponentCoordinates.mat');
s = load(fullMatFileName)
C= struct2cell(s)
%A= cell2mat(C);
%dlmwrite('myFile.txt', s, 'delimiter','\t')
dlmwrite('NodeCoordinates.txt', C, 'delimiter','\t')
