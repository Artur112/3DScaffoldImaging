figure;
numpics = 15;
%celldata.BBox = celldata2.BBox{1,1}
for n = 1:numpics
    if(numpics<6)
        subplot(1,5,n);
    elseif(numpics<11)
        subplot(2,5,n);
    else
        subplot(3,5,n);
    end
    
    %Take bounding box region only
    blob = image(floor(celldata.BBox(n,2)):ceil(celldata.BBox(n,2)) + celldata.BBox(n,4),floor(celldata.BBox(n,1)):ceil(celldata.BBox(n,1)) + celldata.BBox(n,3));
    imshow(blob);
    hold on;
    cc = bwconncomp(blob); %Connected Components
    centroid = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Get centroid of cell in bounding box area
    minorax = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))'); %Minor Axis Length
    majorax = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))'); %Major Axis Length
    orient = cell2mat(struct2cell(regionprops(cc,'Orientation'))'); %Orientation of Major Axis

    deltax = majorax*cosd(orient);
    deltay = majorax*sind(orient);
    %xMaj = [centroid(1) - deltax/2, centroid(1) + deltax/2];
    %yMaj = [centroid(2) + deltay/2, centroid(2) - deltay/2];

    xMaj = centroid(1) + [-1 1]*(majorax/2)*cosd(orient);
    yMaj = centroid(2) + [1 -1]*(majorax/2)*sind(orient);

    deltax = minorax/(2*sind(orient));
    deltay = minorax*sind(90-orient);
    %xMin = [centroid(1) + deltax/2, centroid(1) - deltax/2];
    %yMin = [centroid(2) + deltay/2, centroid(2) - deltay/2];

    xMin= centroid(1) + [-1 1]*(minorax/2)*sind(orient);
    yMin= centroid(2) + [-1 1]*(minorax/2)*cosd(orient);

    line(xMaj,yMaj,'Color','r','LineWidth', 2);
    line(xMin,yMin,'Color','b','LineWidth', 2);

    %Plot Ellipsoid 
    a = 1/2*majorax;
    b = 1/2*minorax;
    t = linspace(0,2*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(yMaj(2)-yMaj(1),xMaj(2)-xMaj(1));
    x = (xMaj(1)+xMaj(2))/2 + X*cos(w) - Y*sin(w);
    y = (yMaj(1)+yMaj(2))/2 + X*sin(w) + Y*cos(w);
    plot(x,y,'g','LineWidth',2)
    hold off;
    title("W/H Ratio: " + num2str(majorax/minorax));
end
legend('MajorAxis','MinorAxis');
sgtitle('Finding cellular Width to Height ratios');