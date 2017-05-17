function [] = Plotme(Plotmeinput)

%% Define Mode.............................................................

x                   = Plotmeinput.x;
y                   = Plotmeinput.y;
xlabelname          = Plotmeinput.xlabelname;
ylabelname          = Plotmeinput.ylabelname;
filename            = Plotmeinput.filename;
nVar                = Plotmeinput.nVar;
titlename           = Plotmeinput.title;
legendoption        = Plotmeinput.legendoption;
legendnames         = Plotmeinput.legendnames;

%% Color Map...............................................................
colororder = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.75  0.75  0.00
    0.25  0.25  0.25
    0.80  0.50  0.50
    0.95  0.95  0.00
    0.25  0.25  0.75
    0.75  0.75  0.75
    0.00  1.00  0.00
    0.76  0.57  0.17
    0.54  0.63  0.22
    0.34  0.57  0.92
    1.00  0.50  0.30
    0.10  0.49  0.47
    0.94  0.80  0.23
    0.66  0.34  0.65
    0.5   0.5   0.5
    ];

%% Plot....................................................................

h_fig = figure;

for i = 1 : nVar
    
    plot(x{1,i}, y{1, i}, 'Color', colororder(i,:),...
        'Linewidth', 2);
    hold on
end

hold off
h_xlabel = xlabel(xlabelname);
h_ylabel = ylabel(ylabelname);

title(titlename);

grid on;
grid minor;

set(gca,'FontSize',16)
set(h_xlabel,'FontSize',20);
set(h_ylabel,'FontSize',20);
set(h_fig, 'Position', [100, 100, 1000, 700]);

ax = gca;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.GridAlpha = 0.7;
ax.MinorGridAlpha = 0.7;

if strcmp(legendoption, 'On') == 1
    legend(legendnames, 'Location', 'southeast','Fontsize',16);
end

print(filename,'-depsc') % Saves the file in EPS vector graphics format
print(filename,'-dpng') % Saves the file in the PNG Bitmap graphics format
print(filename,'-dpdf') % Saves the file in the PDF format
clear h_xlabel h_ylabel ax

%% Developed by Suraj Bansal, (Distributed to AE419 students on 02/07/2017)
%  Graduate Research Assistant
%  M.S., Aerospace Engineering
%  University of Illinois at Urbana-Champaign