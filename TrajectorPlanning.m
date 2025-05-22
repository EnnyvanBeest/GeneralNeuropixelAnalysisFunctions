%% Acute Neuropixels, different insertion sites

figure;
x = [-5000:10:5000-10];
y = [-5000:10:5000-10];
imagesc(x,y,ones(800,800));

% ML; AP
Center = [4557,-2973];
Radius = 500;
points = [4673,-3230;4693,-2940;4474,-2681;4390,-3189];

%draw circle
 drawcircle('Center',Center,'Radius',Radius,'FaceAlpha',0.1,'Color',[1 0 0],'LineWidth',2)
hold on
cols = lines(length(points));
legendname={};
for pid=1:length(points)
    h(pid)=plot(points(pid,1),points(pid,2),'*','color',cols(pid,:));
    legendname={legendname{:} ['Insertion ' num2str(pid)]};
    Distance(pid) = sqrt((points(pid,1)-Center(1))^2+(points(pid,2)-Center(2))^2);
end
legend([h(:)],legendname)
set(gca,'XDir','normal','YDir','normal')
%Zoom in on craniotomy
xlim([Center(1)-Radius*1.1 Center(1)+Radius*1.1])
ylim([Center(2)-Radius*1.1 Center(2)+Radius*1.1])

ylabel('AnteriorPosterior')
xlabel('MedialLateral')
box off

Distance

%% Chronic Neuropixels - distance between probes
% AP - ML (Relative to Bregma) - Given coordinates are for Shank 0
Shank0_Probe0 = [-1650 680]; % Probe 0, Shank 0 position
Shank0_Probe1Prov = [-2480 -2460]; % Probe 1, Shank 0 estimated position

Yaw = 244; % Yaw in trajectory planner coordinates
Radius = 750;

% Correct yaw to match MATLAB’s coordinate system
CorrectedYaw_Probe0 = Yaw + 90; % Adjust from 9 o'clock reference frame to standard
CorrectedYaw_Probe1 = Yaw + 270; % Mirrored yaw for Probe 1
YawRad_Probe0 = deg2rad(CorrectedYaw_Probe0); % Convert to radians
YawRad_Probe1 = deg2rad(CorrectedYaw_Probe1); % Convert to radians

% Shank separation
ShankSep = 250; % Distance between shanks
NumShanks = 4; % Number of shanks

% **Ensure the craniotomy centers align symmetrically**
ShankRecenter = -(Radius*2 - (ShankSep*(NumShanks-1))) / 2;
Craniotomy0 = Shank0_Probe0 + [ShankRecenter*cos(YawRad_Probe0) ShankRecenter*sin(YawRad_Probe0)];

% **Compute symmetrically spaced shank offsets**
ShankOffsets = linspace(-ShankSep * (NumShanks - 1) / 2, ShankSep * (NumShanks - 1) / 2, NumShanks);

% Compute AP and ML shifts relative to the craniotomy center
AP_Shift = ShankOffsets * cos(YawRad_Probe0); % AP component scales correctly
ML_Shift = ShankOffsets * sin(YawRad_Probe0); % ML component scales correctly
ShankPositions = [AP_Shift', ML_Shift']; % [AP, ML] relative to craniotomy center

% Compute rotated shank positions for Probe 0
RotatedShanks_Probe0 = ShankPositions + Craniotomy0;

% **Find Craniotomy1 close to Shank0_Probe1Prov**
% Compute vector from Craniotomy0 to estimated Probe 1 Shank 0 position
ProbeVector = RotatedShanks_Probe0(4,:)-Shank0_Probe1Prov;

% Project ProbeVector onto the yaw axis of Probe 1
ProjectedDist = norm(ProbeVector);

% Generate candidate positions in a **circle** around Craniotomy0
theta = linspace(0, 2*pi, 100); % 100 candidate points on a circle
CircleCandidates = Craniotomy0 + ProjectedDist * [cos(theta') sin(theta')];

% Find the candidate that is **closest to the estimated Probe 1 position**
[~, idx] = min(vecnorm(CircleCandidates - Shank0_Probe1Prov, 2, 2));
Craniotomy1 = CircleCandidates(idx, :); % Best candidate for Craniotomy1
% **Ensure mirroring of shanks across probes**
MirroredShankPositions = flipud(ShankPositions); % Flip shanks to maintain mirroring

% Compute rotated shank positions for Probe 1
RotatedShanks_Probe1 = MirroredShankPositions + Craniotomy1;

% Compute corrected probe distance
CorrectedDistance = norm(Craniotomy1 - Craniotomy0);

% **Plotting**
figure;
x = -5000:10:5000-10;
y = -5000:10:5000-10;
imagesc(x, y, zeros(800, 800));
set(gca, 'clim', [-1 1]);
colormap(gray);
hold on;

% Draw craniotomies
drawcircle('Center', fliplr(Craniotomy0), 'Radius', Radius, 'FaceAlpha', 0.1, 'Color', [1 0 0], 'LineWidth', 2);
drawcircle('Center', fliplr(Craniotomy1), 'Radius', Radius, 'FaceAlpha', 0.1, 'Color', [0 0 1], 'LineWidth', 2);

% Plot shanks
scatter(RotatedShanks_Probe0(:,2), RotatedShanks_Probe0(:,1), 100, 'r', 'filled');
scatter(RotatedShanks_Probe1(:,2), RotatedShanks_Probe1(:,1), 100, 'b', 'filled');

% Draw a line connecting the two probe craniotomy centers
plot([Craniotomy0(2), Craniotomy1(2)], [Craniotomy0(1), Craniotomy1(1)], 'k--', 'LineWidth', 2);

% Annotate the distance
midpoint = (Craniotomy0 + Craniotomy1) / 2;
text(midpoint(2), midpoint(1), sprintf('%.1f µm', CorrectedDistance), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

ylabel('Anterior-Posterior');
xlabel('Medial-Lateral');
box off;
axis equal;
title('Neuropixels Shank Placement with Aligned & Mirrored Shanks');
legend('Probe0 Craniotomy', 'Probe1 Craniotomy (Aligned)', 'Shanks Probe0', 'Shanks Probe1', 'Probe Distance');
set(gca, 'ydir', 'normal');

% **Output results**
disp('Chosen Probe0 Craniotomy Center (AP, ML):');
disp(Craniotomy0);

disp('Shank positions for Probe 0 (AP, ML):');
disp(RotatedShanks_Probe0);

disp('Rotation Probe 0 (degrees):');
disp(Yaw);

disp('Chosen Probe1 Craniotomy Center (AP, ML):');
disp(Craniotomy1);

disp('Shank positions for Probe 1 (AP, ML):');
disp(RotatedShanks_Probe1);

disp('Rotation Probe 1 (degrees):');
disp(mod(Yaw + 180, 360));

disp(['Corrected distance between probe craniotomies: ', num2str(CorrectedDistance), ' µm']);

disp('Add 90 degrees for Neuropixels trajectory planner to match rotation');
