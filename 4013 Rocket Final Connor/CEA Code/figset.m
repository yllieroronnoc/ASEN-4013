function figset(FontSize,FigLineWidth,FontWeight,MarkerSize, MarkerLineWidth, Xresolution, Yresolution, vargin)
% Formats the currently selected matlab figure
% Travis Sippel
% tsippel@purdue.edu
% February 12, 2009

% vargin = h PlotLineWidth


% FigPosition = x y position and width height (in pixels) of figure
left = 100; %pixels
bottom = 100; %pixels
width = Xresolution; %pixels
height = Yresolution; %pixels
FigPosition = [left, bottom, width, height];



% % % % input parameters
% h is the plot handle. when plotting use: h = plot(x,y).
FS = FontSize;% Font size in points
FLW = FigLineWidth;%line width in points
FW = FontWeight;%font weight, either b -bold, n -normal l -light
% PLW = PlotLineWidth; %plot object line width in points


% Set line specification for each line on the current plot
LineHandles = get(gca,'children');
for i = 1:length(LineHandles);
    L = LineHandles(i);
    set(L, 'MarkerSize', MarkerSize);
    set(L, 'LineWidth', MarkerLineWidth);
end;

% figure line width
set(gca,'LineWidth',FLW);

% figure color white
set(gcf,'Color','white');

% legend font size
set(findobj(gcf,'Type','text'),'FontSize',FS);

% ylabel font size
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',FS);

% xlabel font size
h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',FS);

% set axes number font size
set(gca,'FontSize',FS);

% set tile font size
h_title = get(gca,'Title');
set(h_title,'FontSize',FS);

% Apply figure resizing and position
set(gcf,'Position',FigPosition)

% % % % % % % % % % % % % % % % % % % % % 
% APPLY TEXT BOLDING
% % % % % % % % % % % % % % % % % % % % % 
if FW=='b'
    set(h_xlabel,'FontWeight','bold')
    set(h_ylabel,'FontWeight','bold')
    set(h_title,'FontWeight','bold')
    set(findobj(gcf,'Type','text'),'FontWeight','bold');
    set(gca,'FontWeight','bold');
elseif FW=='n'
    set(h_xlabel,'FontWeight','normal')
    set(h_ylabel,'FontWeight','normal')
    set(h_title,'FontWeight','normal')
    set(findobj(gcf,'Type','text'),'FontWeight','normal');
    set(gca,'FontWeight','normal');
elseif FW=='l'
    set(h_xlabel,'FontWeight','light')
    set(h_ylabel,'FontWeight','light')
    set(h_title,'FontWeight','light')
    set(findobj(gcf,'Type','text'),'FontWeight','light');
    set(gca,'FontWeight','light');
end;

        
% % % % % % % % % % % % % % % % % % % % % 
