function placeFigure()
    % Get handles of all open figures
    figHandles = findall(groot, 'Type', 'figure');
    
    % Define the default figure size (width, height)
    defaultSize = [560, 400];  % Default MATLAB figure size
    
    % Define screen size and figure spacing
    screenSize = get(groot, 'ScreenSize');
    xSpacing = 25; % Horizontal spacing between figures
    ySpacing = 100; % Vertical spacing between figures
    
    % Compute maximum rows and columns that can fit on the screen
    maxCols = floor((screenSize(3) - xSpacing) / (defaultSize(1) + xSpacing));

    % Calculate the next figure position
    nextPosition = calculateNextPosition(figHandles, maxCols, defaultSize, xSpacing, ySpacing, screenSize);
    
    % Create the new figure at the computed position
    figure('Position', nextPosition);
end

% Helper function to calculate the next position
function nextPos = calculateNextPosition(figHandles, maxCols, defaultSize, xSpacing, ySpacing, screenSize)
    % Get the positions of all existing figures
    existingPositions = arrayfun(@(fh) get(fh, 'Position'), figHandles, 'UniformOutput', false);
    existingPositions = cat(1, existingPositions{:});
    
    % Calculate grid position for the next figure
    numFigures = size(existingPositions, 1);
    col = mod(numFigures, maxCols);  % Column index
    row = floor(numFigures / maxCols);  % Row index
    
    % Calculate the new position based on row and column
    left = col * (defaultSize(1) + xSpacing) + xSpacing;
    bottom = screenSize(4) - (row + 1) * (defaultSize(2) + ySpacing);
    
    % Set the position for the next figure
    nextPos = [left, bottom, defaultSize];
end
