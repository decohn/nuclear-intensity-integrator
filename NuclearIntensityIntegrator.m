%% Nuclear Intensity Integrator
% See Github repository for README and more details regarding the project.

% Change background intensity calculation to function on a cell-by-cell basis.
% Try segmenting the image into individual cells first? Otherwise, I'm not
% sure how you could easily pair nuclei to cells.      

clc
clear variables

%% Configuration Variables
figureNumber = 1;
intensityPerFPMolecule = 155;
stackSize = 17;
nucleusSizeTolerance = 2.5;
imageNames = {'POL32_1s_100%_exp2_1', 'POL32_1s_100%_exp2_2', 'POL32_1s_100%_exp2_3', 'POL32_1s_100%_exp2_4'};
proteinName = 'POL32';
useCFPChannel = 1;
useMaxProjection = 0;
planesToIgnoreInProjections = 4;
numImages = size(imageNames, 2);
numThresholds = 100;

% The below value, in pixels, is used in order to determine which
% thresholding sensitivity value is optimal.
desiredNuclearArea = 595;

allCopyNumbers = zeros(numImages,1);
copyNumbersByObject = cell(numThresholds+1, numImages);

% This is the main control loop, which will perform analysis for each
% separate .tif file.
for k=1:numImages
    %% Image Read
    
    % Now, read the GFP channel plane-by-plane.
    I_GFP = uint16(readStack(proteinName, imageNames, k, stackSize, 0));
    
    if(useMaxProjection == 0)
        % Finds average projection of the channel that will be used for
        % thresholding.
        if(useCFPChannel == 0)
            I_GFP_temp = I_GFP(:,:,(1+planesToIgnoreInProjections):(stackSize-planesToIgnoreInProjections));
            Ithresh = sum(I_GFP_temp, 3);
        else
            I_CFP = readStack(proteinName, imageNames, k, stackSize, 1);
            I_CFP_temp = I_CFP(:,:,(1+planesToIgnoreInProjections):(stackSize-planesToIgnoreInProjections));
            Ithresh = sum(I_CFP_temp, 3);
        end
    
        Ithresh = uint16(Ithresh ./ (stackSize-(2*planesToIgnoreInProjections)));
    else
        % Finds maximum projection of the channel that will be used for
        % thresholding.
        if(useCFPChannel == 0)
            I_GFP_temp = I_GFP(:,:,(1+planesToIgnoreInProjections):(stackSize-planesToIgnoreInProjections));
            Ithresh = uint16(max(I_GFP_temp, [], 3));
        else
            I_CFP = uint16(readStack(proteinName, imageNames, k, stackSize, 1));
            I_CFP_temp = I_CFP(:,:,(1+planesToIgnoreInProjections):(stackSize-planesToIgnoreInProjections));
            Ithresh = uint16(max(I_CFP_temp, [], 3));
        end
    end
        
    %% Segment Image into Individual Nuclei
    MAXbinary = cell(numThresholds + 1, 1);
    allBackground = cell(numThresholds + 1, 1);
    
    % Creates a bunch of binary masks, using different sensitivity values.
    for i=0:1:numThresholds
        MAXbinary{i+1} = segmentNuclei(Ithresh, i/numThresholds);
    end
    
    % Displays the result of the initial thresholding. Make a smartImShow
    % method that does this automatically.
    % figure(figureNumber);
    % imshow(MAXbinary{10},'InitialMagnification','Fit');
    % figureNumber = figureNumber + 1;
    
    % Filters each binary mask, collects background pixels, and finds
    % the median area of every object that is left for each filtered
    % image. Selects (and continues with) the filtered image whose median
    % object area is closest to the desired area. This is how the program
    % decides what sensitivity value to ultimately use for the
    % thresholding.
    medianArea = zeros(i,1);
    
    for i=0:1:numThresholds
        [MAXbinary{i+1}, allBackground{i+1}, medianArea(i+1)] = filterNuclei(MAXbinary{i+1}, nucleusSizeTolerance);
    end
    
    % finds the index of MAXbinary that is best thresholded
    [M,I] = min(abs(medianArea - desiredNuclearArea));
    I = I(1);
    
    % Uncomment when testing filtering.
    % figure(figureNumber);
    % imshow(MAXbinary{I},'InitialMagnification','Fit');
    % figureNumber = figureNumber + 1;
     
    %% Compute Background Intensity Value for each Plane
    backgroundInt = zeros(i,1);
    
    for i=1:stackSize     
        [backgroundInt(i), figureNumber] = getBkgrndInt(allBackground{I}, I_GFP(:,:,i), 0.35, figureNumber);          
    end
    
    %% Compute the Integrated Intensity in each Plane for each Nucleus
    totalIntensity = 0;
        
    % The purpose of this block is just to determine how many objects
    % are in the binary image.
    intensityStruct = regionprops(MAXbinary{I}, I_GFP(:,:,1), 'PixelValues');    
    pixelValues = table2array(struct2table(intensityStruct));    
    numberOfObjects = size(pixelValues, 1);
        
    integratedIntensitiesPerPlane = zeros(stackSize, numberOfObjects);

    for i = 1:stackSize
        integratedIntensitiesPerPlane(i, :) = integrateIntensityPerObject(MAXbinary{I}, I_GFP(:,:,i), backgroundInt(i));
    end

    objectMaxIntensities = zeros(numberOfObjects, 1);

    for i = 1:numberOfObjects
        objectIntensity = max(integratedIntensitiesPerPlane(:, i));
        totalIntensity = totalIntensity + objectIntensity;
        objectMaxIntensities(i, 1) = objectIntensity;
    end

    copyNumbersByObject{k} = objectMaxIntensities ./ intensityPerFPMolecule;

    allCopyNumbers(k) = totalIntensity / numberOfObjects / intensityPerFPMolecule;
end
% The above is the end of the main control loop.
    
%% Output Results

% for quick and dirty output
disp(mean(allCopyNumbers, 1));

figure(figureNumber);
figureNumber = figureNumber + 1;
    
sens = (I-1) ./ numThresholds;
makePlotsWithSeparateColumnsForEachImage(copyNumbersByObject, numImages, proteinName, sens, useCFPChannel);

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
% maxProjection currently not used
function maxProjection = readMax(protein, imageNames, imageNum, boolCFP)
    if(boolCFP == 0)
        prefix = '/MAX_';
    else
        prefix = '/MAX_CFP_';
    end
    
    maxProjection = imread([pwd, '/', protein, prefix, imageNames{imageNum}, '.tif']);   
end
function imagePlanes = readStack(protein, imageNames, imageNum, numPlanes, boolCFP)
    if(boolCFP == 1)
        suffix = '_CFP.tif';
    else
        suffix = '.tif';
    end
    
    testPlane = imread([pwd, '/', protein, '/', imageNames{imageNum}, suffix], 1);
    
    imagePlanes = zeros([size(testPlane), numPlanes]);
    imagePlanes(:,:,1) = testPlane;

    for i=2:numPlanes
        imagePlanes(:,:,i) = imread([pwd, '/', protein, '/', imageNames{imageNum}, suffix], i);
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
function [filteredBinaryMask, backgroundPixels, medianArea] = filterNuclei(binaryMask, sizeTolerance)
    medianArea = -1;

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
    medianArea = median(allAreas);
end
function [meanBackgroundIntensity, newFigNum] = getBkgrndInt(pixelLocations, greyscaleImage, sensitivity, figNum)
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
    
    %figure(figNum)
    %imshow(imadjust(cellularBackgrounds), 'InitialMagnification', 'Fit');
    newFigNum = figNum;
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
         pointsPerImage(i) = size(copyNumbers{i}, 1);
    end

    xvals = [];
    yvals = [];
    means = zeros(numberOfImages,1);
    stdevs = zeros(numberOfImages,1);

    for i=1:numberOfImages
        xvals = [xvals; i*ones(pointsPerImage(i),1) - 0.125 + 0.25 * rand(pointsPerImage(i),1)];
        yvals = [yvals; copyNumbers{i}];
        means(i) = mean(copyNumbers{i});
        stdevs(i) = std(copyNumbers{i}) / sqrt(pointsPerImage(i)); 
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
    title([protein, ' Copy # per Nucleus using ', channel, ' channel.']);
    errorbar(1:numberOfImages, means, stdevs, '+m', 'LineWidth', 3, 'CapSize', 20); 
end