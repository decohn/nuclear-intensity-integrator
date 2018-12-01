%% Nuclear Intensity Integrator
% See Github repository for README and more details regarding the project.

% It looks like the main issue with this code is probably going to be
% getting the thresholding of the nuclei correct. I think that using a max
% projection has to be the correct way to do that, but I don't know if
% there are better algorithms than graythresh and adaptthresh available to
% me. I could always threshold a brightfield image, but the
% issue is that the cells aren't really all in the same focal plane, and
% there's that ugly line running through the middle.

% Replace all my 17s with the stackSize variable.

clc
clear variables

%% Configuration Variables
intensityPerFPMolecule = 3500;
stackSize = 17;
nucleusSizeTolerance = 2.5;
imageNames = {'PCNA_1s_100%_1', 'PCNA_1s_100%_2', 'PCNA_1s_100%_3', 'PCNA_1s_100%_4'};
numImages = size(imageNames, 2);
allCopyNumbers = cell(numImages,1);

for k=1:numImages
    %% Image Read

    % Read a max-projection
    Imax = imread([pwd, '/MAX_', imageNames{k}, '.tif']);

    % Read the image plane-by-plane
    I = cell(stackSize,1);

    for i=1:stackSize
        I{i} = imread([pwd, '/', imageNames{k}, '.tif'], i);
    end

    %% Segment Image into Individual Nuclei
    % Multiplying graythresh by 1.2 seems plausible, as does using adaptthresh
    % with a sensitivity of 0.15. I think the latter (for the most part) looks
    % better.
    threshold = adaptthresh(Imax,0.15);
    MAXbinary = imbinarize(Imax, threshold);

    % Uncomment when testing thresholding.
    % figure(2*k-1);
    % imshow(MAXbinary,'InitialMagnification','Fit');

    % remove any objects that contact the border of the image
    MAXbinary = imclearborder(MAXbinary, 4);

    % fill in any holes in objects
    MAXbinary = imfill(MAXbinary, 'holes');

    % make measurements for all objects remaining in the image
    cc = bwconncomp(MAXbinary);
    stats = regionprops(cc, 'Area','Eccentricity');
    data = table2array(struct2table(stats));

    % remove any objects that are too big, too small, or not round enough
    allAreas = sort(data(:,1));
    allAreas = allAreas(allAreas > 10);
    medianArea = median(allAreas);
    idx = find([stats.Area] > (medianArea / nucleusSizeTolerance) & [stats.Area] < (medianArea * nucleusSizeTolerance) & [stats.Eccentricity] < 0.8); 
    MAXbinary = ismember(labelmatrix(cc), idx);

    
    % Uncomment when testing thresholding.
    % figure(2*k);
    % imshow(MAXbinary,'InitialMagnification','Fit');

    %% Compute Background Intensity Value for each Plane
    % The values that I'm getting are 900-1000, which seem on the high side.
    % Just by inspection, 'true background' seems to be in the vicinity of 425.
    % It's really hard for me to decide whether I should be aiming for that
    % true background, or whether it's ok that I have some bleed from the
    % nucleus included in it. In reality, the ideal would be to know the
    % background autofluorescence in cells that didn't have PCNA (this would 
    % have to be a fair bit higher than 425, so maybe my value is ok-ish).

    backgroundInt = zeros(stackSize,1);

    for i=1:stackSize
        backgroundRegions = uint16(imcomplement(MAXbinary));
        numBackgroundPixels = sum(sum(backgroundRegions));
        backgroundInt(i) = sum(sum(backgroundRegions .* I{i})) / numBackgroundPixels;
    end

    %% Compute the Integrated Intensity in each Plane for each Nucleus
    intensityStructs = cell(stackSize,1);
    pixelValueArrays = cell(stackSize,1);

    for i=1:stackSize
        intensityStructs{i} = regionprops(MAXbinary, I{i}, 'PixelValues');
    end

    for i=1:stackSize
        pixelValueArrays{i} = table2array(struct2table(intensityStructs{i}));
    end

    % for each object, go through and calculate the intensity in each plane

    copyNumbers = zeros(size(pixelValueArrays{1},1), 1);

    for i=1:size(pixelValueArrays{1},1)
        highestIntensity = 0;

        for j=1:stackSize
            pixels = pixelValueArrays{j}{i};
            intensityInThisPlane = sum(pixels) - (backgroundInt(j) * size(pixels,1));
            highestIntensity = max([highestIntensity, intensityInThisPlane]);
        end

        copyNumbers(i) = highestIntensity / intensityPerFPMolecule;
    end
    
    allCopyNumbers{k} = copyNumbers;
end

%% Output Results
% pretty these up SIGNIFICANTLY, and find good ways to display the mean or
% median. consider doing the 1D scatter plot thingy that papers seem to
% love doing

figure(2*numImages)

for i=1:numImages
    pointsPerImage(i) = size(allCopyNumbers{i}, 1);
end

xvals = [];
yvals = [];
means = zeros(4,1);
stdevs = zeros(4,1);

for i=1:numImages
    xvals = [xvals; i*ones(pointsPerImage(i),1) - 0.125 + 0.25 * rand(pointsPerImage(i),1)];
    yvals = [yvals; allCopyNumbers{i}];
    means(i) = mean(allCopyNumbers{i});
    stdevs(i) = std(allCopyNumbers{i}); 
end

scatter(xvals, yvals, 20, 'filled');
hold on

% Plot standard error bars
xlim([0.5, numImages + 0.5]);
xlabel('Image Number');
ylabel('PCNA Copy Number');
xticks(1:numImages);
title('PCNA Copy Number per Nucleus');

errorbar(1:numImages, means, stdevs, '+m', 'LineWidth', 3, 'CapSize', 20);

%% Debug Thresholding
%BWoutline = bwperim(MAXbinary); 
%Imax(BWoutline) = 0; 
%figure(3) 
%imshow(imadjust(Imax),'InitialMagnification','Fit');
