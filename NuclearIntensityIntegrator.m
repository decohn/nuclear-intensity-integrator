%% Nuclear Intensity Integrator
% See Github repository for README and more details regarding the project.

% It looks like the main issue with this code is probably going to be
% getting the thresholding of the nuclei correct. I think that using a max
% projection has to be the correct way to do that, but I don't know if
% there are better algorithms than graythresh and adaptthresh available to
% me. I could always threshold a brightfield image, but the
% issue is that the cells aren't really all in the same focal plane, and
% there's that ugly line running through the middle.

% Figure out whether or not standard deviation or standard error is the
% appropriate thing to use for my error bars.

% Change background intensity calculation to function on a cell-by-cell basis.

clc
clear variables

%% Configuration Variables
intensityPerFPMolecule = 155;
stackSize = 17;
nucleusSizeTolerance = 2.5;
imageNames = {'PCNA_1s_100%_1','PCNA_1s_100%_2','PCNA_1s_100%_3','PCNA_1s_100%_4'};
proteinName = 'PCNA';
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
    % Optimal sensitivity value seems to vary by protein. Any way to do
    % this programmatically? For PCNA: 0.15. For Mcm4: 0.25. PCNA
    % sensitivity value honestly might be a little on the large side? But
    % how on Earth could you possibly quantify that sort of thing? This is
    % why I don't like having a sensitivity parameter in the first place.
    % Other effective ways to do this thresholding?
    threshold = adaptthresh(Imax,0.25);
    MAXbinary = imbinarize(Imax, threshold);

    % Uncomment when testing thresholding.
    figure(4*k-3);
    imshow(MAXbinary,'InitialMagnification','Fit');
    
    % fill in any holes in objects
    MAXbinary = imfill(MAXbinary, 'holes');
    
    % collect all background pixels by obtaining a copy of the original
    % image where the nuclei have their pixel intensity set to zero
    allBackground = uint16(imcomplement(MAXbinary));

    % remove any objects that contact the border of the image
    MAXbinary = imclearborder(MAXbinary, 4);

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
    figure(4*k-2);
    imshow(MAXbinary,'InitialMagnification','Fit');
    
    %% Compute Background Intensity Value for each Plane

    backgroundInt = zeros(stackSize,1);
    
    for i=1:stackSize
        allBackgroundPixels = allBackground .* I{i};
        
        % set this sensitivity as high as possible without getting a bunch
        % of background pixels included in the foreground. Some room for
        % improvement here: is there a way to crank up the sensitivity to
        % get more of the cell body, and then discard the isolated pixels
        % that lie in the background? Note that unfortunately the regions
        % found in cellularBackgrounds cannot be filled, because then they
        % will actually fill in the nuclei in their centres! Looking at the
        % thresholding results, I don't even think that filling would help
        % (even if the nuclei weren't an issue). The thresholding genuinely
        % seems extremely good.
        threshold = adaptthresh(allBackgroundPixels, 0.35);
        cellularBackgrounds = uint16(imbinarize(allBackgroundPixels, threshold));
        cellularBackgroundsGray = cellularBackgrounds .* I{i};
        
        cellularBackgroundsGray = cellularBackgroundsGray(cellularBackgroundsGray > 0);
        
        backgroundInt(i) = mean(mean(cellularBackgroundsGray));
        
        if(i == 9)
            figure(4*k-1)
            imshow(imadjust(cellularBackgrounds), 'InitialMagnification', 'Fit');
            BWoutline = bwperim(cellularBackgrounds);
            allBackgroundPixels(BWoutline) = 0;
            %imshow(imadjust(allBackgroundPixels),'InitialMagnification','Fit');
        end
        
        % This is the old background computation method, which just
        % averaged together the values of all pixels that were not in the
        % nucleus. 
        %backgroundRegions = uint16(imcomplement(MAXbinary));
        %numBackgroundPixels = sum(sum(backgroundRegions));
        %backgroundInt(i) = sum(sum(backgroundRegions .* I{i})) / numBackgroundPixels;
    end

    BWoutline = bwperim(MAXbinary); 
    Imax(BWoutline) = 0; 
    figure(4*k) 
    imshow(imadjust(Imax),'InitialMagnification','Fit');
    
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

figure(4*numImages + 1)

for i=1:numImages
    pointsPerImage(i) = size(allCopyNumbers{i}, 1);
end

xvals = [];
yvals = [];
means = zeros(numImages,1);
stdevs = zeros(numImages,1);

for i=1:numImages
    xvals = [xvals; i*ones(pointsPerImage(i),1) - 0.125 + 0.25 * rand(pointsPerImage(i),1)];
    yvals = [yvals; allCopyNumbers{i}];
    means(i) = mean(allCopyNumbers{i});
    stdevs(i) = std(allCopyNumbers{i}) / sqrt(pointsPerImage(i)); 
end

scatter(xvals, yvals, 20, 'filled');
hold on

% Plot standard error bars
xlim([0.5, numImages + 0.5]);
xlabel('Image Number');
ylabel([proteinName ,' Copy Number']);
xticks(1:numImages);
title([proteinName, ' Copy Number per Nucleus']);
errorbar(1:numImages, means, stdevs, '+m', 'LineWidth', 3, 'CapSize', 20);

% Scatter plot with all points in one column
figure(4*numImages+2)
scatter(ones(sum(pointsPerImage), 1) - 0.125 + 0.25 * rand(sum(pointsPerImage), 1), yvals, 20, 'filled');
hold on
xlim([0.8,1.2]);
ylabel([proteinName ,' Copy Number']);
xticks([]);
title([proteinName, ' Copy Number per Nucleus']);
errorbar(1, mean(yvals), std(yvals) / sqrt(sum(pointsPerImage)), '+m', 'LineWidth', 3, 'CapSize', 20);

%% Results
% Median copy # for PCNA with 0.15 sensitivity and 155 intensity per FP
% molecule, with smart adjustments for the cellular autofluorescence is 
% 9140. The mean copy # is 11384.
% Median copy # for Mcm4 with 0.25 sensitivity and 155 intensity per FP
% molecule, with smart adjustments for the cellular autofluorescence is 
% 662. The mean copy # is 696.