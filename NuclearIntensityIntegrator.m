%% Nuclear Intensity Integrator
% See Github repository for README and more details regarding the project.

% It looks like the main issue with this code is probably going to be
% getting the thresholding of the nuclei correct. I think that using a max
% projection has to be the correct way to do that, but I don't know if
% there are better algorithms than graythresh and adaptthresh available to
% me. 

% Change background intensity calculation to function on a cell-by-cell basis.
% Try segmenting the image into individual cells first? Otherwise, I'm not
% sure how you could easily pair nuclei to cells. Compare YDC05 results
% using CFP as the nuclear thresholder to those using GFP. Produce a
% comparison table with all the different sensitivity values, etc.      

clc
clear variables

%% Configuration Variables
figureNumber = 1;
intensityPerFPMolecule = 155;
stackSize = 17;
nucleusSizeTolerance = 2.5;
imageNames = {'YDC05_1s_100%_1', 'YDC05_1s_100%_2', 'YDC05_1s_100%_3', 'YDC05_1s_100%_4'};
proteinName = 'YDC05';
useCFPChannel = 0;
numImages = size(imageNames, 2);
numThresholds = 20;
allCopyNumbers = zeros(numImages,numThresholds+1);
copyNumbersByObject = cell(numThresholds+1, numImages);

% This is the main control loop, which will perform analysis for each
% separate .tif file.
for k=1:numImages
    %% Image Read
    
    % Read a max-projection in the channel that will be used to segment
    % nuclei.
    Imax = readMax(proteinName, imageNames, k, useCFPChannel);
   
    % Now, read the GFP channel plane-by-plane.
    I_GFP = readStack(proteinName, imageNames, k, stackSize);
    
    % attempting to test whether averaging the planes will make
    % thresholding easier. I'm not sure yet, leave this option open.
    Imax = I_GFP{1};
    for i=2:17
        Imax = Imax + I_GFP{i};
    end
    Imax = Imax ./ 17;
    
    %% Segment Image into Individual Nuclei
    % For PCNA: 0.15. For Mcm4: 0.25. For the CFP
    % channel: typically around 0.5, but this varies significantly from
    % image to image, which is a huge issue.
    MAXbinary = cell(numThresholds + 1, 1);
    allBackground = cell(numThresholds + 1, 1);
    
    % Creates a bunch of binary masks, using different sensitivity values.
    for i=0:1:numThresholds
        MAXbinary{i+1} = segmentNuclei(Imax, i/numThresholds);
    end
    
    % Displays the result of the initial thresholding. Make a smartImShow
    % method that does this automatically.
    figure(figureNumber);
    imshow(MAXbinary{10},'InitialMagnification','Fit');
    figureNumber = figureNumber + 1;
    
    % Filters each binary mask and collects background pixels.
    for i=0:1:numThresholds
        [MAXbinary{i+1}, allBackground{i+1}] = filterNuclei(MAXbinary{i+1}, nucleusSizeTolerance);
    end
    
    % Uncomment when testing filtering.
    figure(figureNumber);
    imshow(MAXbinary{10},'InitialMagnification','Fit');
    figureNumber = figureNumber + 1;
     
    %% Compute Background Intensity Value for each Plane
    backgroundInt = zeros(stackSize, numThresholds + 1);
    
    for i=1:stackSize
        for j=1:(numThresholds+1)
            if (i == 9 && j == 10)
                boolDisplay = 1;
            else
                boolDisplay = 0;
            end
            
            [backgroundInt(i, j), figureNumber] = getBkgrndInt(allBackground{j}, I_GFP{i}, 0.35, boolDisplay, figureNumber);          
        end
    end
    
    %% Compute the Integrated Intensity in each Plane for each Nucleus
    for j = 1:(numThresholds+1)
        totalIntensity = 0;
        
        % The purpose of this block is just to determine how many objects
        % are in the binary image.
        intensityStruct = regionprops(MAXbinary{j}, I_GFP{1}, 'PixelValues');    
        pixelValues = table2array(struct2table(intensityStruct));    
        numberOfObjects = size(pixelValues, 1);
        
        integratedIntensitiesPerPlane = zeros(stackSize, numberOfObjects);
        
        for i = 1:stackSize
            integratedIntensitiesPerPlane(i, :) = integrateIntensityPerObject(MAXbinary{j}, I_GFP{i}, backgroundInt(i,j));
        end
        
        objectMaxIntensities = zeros(numberOfObjects, 1);
        
        for i = 1:numberOfObjects
            objectIntensity = max(integratedIntensitiesPerPlane(:, i));
            totalIntensity = totalIntensity + objectIntensity;
            objectMaxIntensities(i, 1) = objectIntensity;
        end
        
        copyNumbersByObject{j, k} = objectMaxIntensities ./ intensityPerFPMolecule;
        
        allCopyNumbers(k, j) = totalIntensity / numberOfObjects / intensityPerFPMolecule;
    end
end
% The above is the end of the main control loop.
    
%% Output Results
for i = 1:(numThresholds+1)
    figure(figureNumber);
    figureNumber = figureNumber + 1;
    
    sens = (i-1) ./ numThresholds;
    makePlotsWithSeparateColumnsForEachImage(copyNumbersByObject(i, :), numImages, proteinName, sens, useCFPChannel);
end

%     
%     % Scatter plot with all points in one column
%     %figure(4*numImages+2)
%     scatter(ones(sum(pointsPerImage), 1) - 0.125 + 0.25 * rand(sum(pointsPerImage), 1), yvals, 20, 'filled');
%     hold on
%     xlim([0.8,1.2]);
%     ylabel([proteinName ,' Copy Number']);
%     xticks([]);
%     title([proteinName, ' Copy Number per Nucleus']);
%     errorbar(1, mean(yvals), std(yvals) / sqrt(sum(pointsPerImage)), '+m', 'LineWidth', 3, 'CapSize', 20);

%% Functions
function maxProjection = readMax(protein, imageNames, imageNum, boolCFP)
    if(boolCFP == 0)
        prefix = '/MAX_';
    else
        prefix = '/MAX_CFP_';
    end
    
    maxProjection = imread([pwd, '/', protein, prefix, imageNames{imageNum}, '.tif']);   
end
function imagePlanes = readStack(protein, imageNames, imageNum, numPlanes)
    imagePlanes = cell(numPlanes,1);

    for i=1:numPlanes
        imagePlanes{i} = imread([pwd, '/', protein, '/', imageNames{imageNum}, '.tif'], i);
    end
end
function maxBinaryMask = segmentNuclei(maxProjection, sensitivity)
    
    if(sensitivity ~= -1)
        threshold = adaptthresh(maxProjection, sensitivity);
    else
        threshold = adaptthresh(maxProjection);
    end
        
    maxBinaryMask = imbinarize(maxProjection, threshold);
end
function [filteredBinaryMask, backgroundPixels] = filterNuclei(binaryMask, sizeTolerance)
    % fill in any holes in objects
    filteredBinaryMask = imfill(binaryMask, 'holes');

    % collect all background pixels by obtaining a copy of the original
    % image where the nuclei have their pixel intensity set to zero
    backgroundPixels = uint16(imcomplement(filteredBinaryMask));

    % remove any objects that contact the border of the image
    filteredBinaryMask = imclearborder(filteredBinaryMask, 4);
    
    % If there are no objects left, end the filtering immediately.
    if(sum(sum(filteredBinaryMask)) == 0)
        return
    end

    % make measurements for all objects remaining in the image
    cc = bwconncomp(filteredBinaryMask);
    stats = regionprops(cc, 'Area','Eccentricity');
    data = table2array(struct2table(stats));

    % remove any objects that are too big, too small, or not round enough
    allAreas = sort(data(:,1));
    allAreas = allAreas(allAreas > 20);
    medianArea = median(allAreas);
    idx = find([stats.Area] > (medianArea / sizeTolerance) & [stats.Area] < (medianArea * sizeTolerance) & [stats.Eccentricity] < 0.8); 
    filteredBinaryMask = ismember(labelmatrix(cc), idx); 
end
function [meanBackgroundIntensity, newFigNum] = getBkgrndInt(pixelLocations, greyscaleImage, sensitivity, boolDisplay, figNum)
    allBackgroundPixels = pixelLocations .* greyscaleImage;

    % set this sensitivity as high as possible without getting a bunch
    % of background pixels included in the foreground. Some room for
    % improvement here: is there a way to crank up the sensitivity to
    % get more of the cell body, and then discard the isolated pixels
    % that lie in the background? 
    threshold = adaptthresh(allBackgroundPixels, sensitivity);
    cellularBackgrounds = uint16(imbinarize(allBackgroundPixels, threshold));
    cellularBackgroundsGray = cellularBackgrounds .* greyscaleImage;

    cellularBackgroundsGray = cellularBackgroundsGray(cellularBackgroundsGray > 0);

    meanBackgroundIntensity = mean(mean(cellularBackgroundsGray));
    
    if (boolDisplay == 1)
        figure(figNum)
        imshow(imadjust(cellularBackgrounds), 'InitialMagnification', 'Fit');
        newFigNum = figNum + 1;
    else
        newFigNum = figNum;
    end
end
function totalIntensities = integrateIntensityPerObject(nucleiLocations, greyScaleImage, backgroundLevel)
    intensityStruct = regionprops(nucleiLocations, greyScaleImage, 'PixelValues');    
    pixelValues = table2array(struct2table(intensityStruct));
    
    if isa(pixelValues, 'cell')
        numberOfObjects = size(pixelValues, 1);
    else
        numberOfObjects = 1;
    end
    
    totalIntensities = zeros(numberOfObjects, 1);
    
    if ~isa(pixelValues, 'double')   
        for i = 1:numberOfObjects
            if numberOfObjects ~= 1
                totalIntensities(i) = sum(pixelValues{i}) - (backgroundLevel * size(pixelValues{i}, 1));
            else
                totalIntensities(i) = sum(pixelValues);
            end
            
            if(totalIntensities(i) < 0)
                totalIntensities(i) = 0;
            end
        end
    end
end
function noReturn = makePlotsWithSeparateColumnsForEachImage(copyNumbers, numberOfImages, protein, sensitivity, boolCFP)
    noReturn = -1;
    
    pointsPerImage = zeros(numberOfImages, 1);
    
    for i=1:numberOfImages
         pointsPerImage(i) = size(copyNumbers{1, i}, 1);
    end

    xvals = [];
    yvals = [];
    means = zeros(numberOfImages,1);
    stdevs = zeros(numberOfImages,1);

    for i=1:numberOfImages
        xvals = [xvals; i*ones(pointsPerImage(i),1) - 0.125 + 0.25 * rand(pointsPerImage(i),1)];
        yvals = [yvals; copyNumbers{1,i}];
        means(i) = mean(copyNumbers{1,i});
        stdevs(i) = std(copyNumbers{1,i}) / sqrt(pointsPerImage(i)); 
    end
 
    scatter(xvals, yvals, 20, 'filled');
    hold on
    
    if(boolCFP == 1)
        channel = 'CFP';
    else
        channel = 'GFP';
    end
    
    % Plot standard error bars
    xlim([0.5, numberOfImages + 0.5]);
    xlabel('Image Number');
    ylabel([protein ,' Copy Number']);
    xticks(1:numberOfImages);
    % Can't concatenate doubles to a string?? is a cast needed?
    title([protein, ' Copy # per Nucleus with ', channel, ' sens=', num2str(sensitivity)]);
    errorbar(1:numberOfImages, means, stdevs, '+m', 'LineWidth', 3, 'CapSize', 20); 
end
%% Results
% Median copy # for PCNA with 0.15 sensitivity and 155 intensity per FP
% molecule, with smart adjustments for the cellular autofluorescence is 
% 9140. The mean copy # is 11384.
% Median copy # for Mcm4 with 0.25 sensitivity and 155 intensity per FP
% molecule, with smart adjustments for the cellular autofluorescence is 
% 662. The mean copy # is 696.
% Median copy # for YDC05 with automated sensitivity and 155 intensity per FP
% molecule, with smart adjustments for the cellular autofluorescence, and
% using the CFP channel for the nuclear thresholding, is 418. Mean is 476. 