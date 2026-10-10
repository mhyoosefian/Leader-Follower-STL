function animateLeaderFollower(t, loc_l, loc_f, P, opts)
% ANIMATELEADERFOLLOWER  Animation of the leader-follower solution.
%
%   animateLeaderFollower(t, loc_l, loc_f, P)        plays the animation.
%   animateLeaderFollower(t, loc_l, loc_f, P, opts)  with options:
%       opts.video = 'leader_follower.mp4'  also saves it as an MP4 file
%       opts.speed = 4                      playback speed (4 = four times real time)
%       opts.fps   = 30                     frames per second
%
%   t     : 1 x (N+1) time steps of the solution (s)
%   loc_l : 2 x (N+1) positions of the leader
%   loc_f : 2 x (N+1) positions of the follower (from the optimization)
%   P     : struct with fields r, omega_l, theta0, a, b, region = [xLeft xRight yBottom yTop],
%           window = [start finish]

if nargin < 5, opts = struct; end
if ~isfield(opts, 'video'), opts.video = ''; end
if ~isfield(opts, 'speed'), opts.speed = 4; end
if ~isfield(opts, 'fps'), opts.fps = 30; end

% Smooth motion between the time steps: the leader moves on its circle and
% the follower moves with constant velocity during each step.
tf = (t(1):opts.speed/opts.fps:t(end));
th = P.theta0 + P.omega_l*(tf - t(1));
L = [P.r*cos(th); P.r*sin(th)];
F = [interp1(t, loc_f(1, :), tf); interp1(t, loc_f(2, :), tf)];
dF = gradient(F, tf(2) - tf(1)); dL = gradient(L, tf(2) - tf(1));
R = P.region; inRegion = F(1, :) >= R(1) - 1e-6 & F(1, :) <= R(2) + 1e-6 & F(2, :) >= R(3) - 1e-6 & F(2, :) <= R(4) + 1e-6;

% ---------------- Figure ----------------
fig = figure('Color', 'w', 'Position', [100 100 1100 880], 'Name', 'Leader-follower');
ax = axes(fig, 'Position', [.07 .17 .9 .77]); hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
axis(ax, [-5.2 3 -3.6 3]); grid(ax, 'on');
set(ax, 'FontSize', 13, 'GridAlpha', .12, 'Layer', 'top');
xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)');
green = [.20 .65 .30]; blue = [.15 .35 .80]; red = [.85 .20 .20]; dark = [.15 .15 .2];

region = patch(ax, R([1 2 2 1]), R([3 3 4 4]), green, 'FaceAlpha', .25, 'EdgeColor', green, 'LineWidth', 1.5);
text(ax, mean(R(1:2)), R(4) + .25, 'goal region', 'HorizontalAlignment', 'center', 'Color', green, 'FontSize', 13);
a = linspace(0, 2*pi, 200);
plot(ax, P.r*cos(a), P.r*sin(a), '--', 'Color', [.6 .6 .65], 'LineWidth', 1);   % leader's path

% Allowed region around the leader (square of half-width b) and keep-out
% square (half-width a)
sq = @(c, h) [c(1) + h*[-1 1 1 -1]; c(2) + h*[-1 -1 1 1]];
s = sq(L(:, 1), P.b); safeBox = patch(ax, s(1, :), s(2, :), blue, 'FaceAlpha', .07, 'EdgeColor', blue, 'LineWidth', 1.2, 'LineStyle', '--');
s = sq(L(:, 1), P.a); keepOut = patch(ax, s(1, :), s(2, :), red, 'FaceAlpha', .18, 'EdgeColor', red, 'LineWidth', 1.2);

trail = plot(ax, nan, nan, '-', 'Color', [blue .6], 'LineWidth', 1.5);
leader = patch(ax, nan, nan, dark, 'EdgeColor', 'none');
follower = patch(ax, nan, nan, blue, 'EdgeColor', 'none');
legend(ax, [leader follower safeBox keepOut region], {'leader', 'follower', ...
    sprintf('follower must stay within %.2f m (\\infty-norm)', P.b), sprintf('keep-out zone (%.1f m)', P.a), 'goal region'}, ...
    'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 12, 'AutoUpdate', 'off', 'Box', 'off');
label = title(ax, '', 'FontSize', 15, 'FontWeight', 'normal');
check = text(ax, R(1), R(3) - .35, '', 'Color', green, 'FontSize', 15, 'FontWeight', 'bold');

% Robots are drawn as arrowheads pointing in their direction of motion.
arrow = [.22 -.13 -.06 -.13; 0 .13 0 -.13];
heading = @(v, prev) atan2(v(2), v(1)) .* (norm(v) > 1e-6) + prev .* (norm(v) <= 1e-6);

if ~isempty(opts.video)
    vw = VideoWriter(opts.video, 'MPEG-4'); vw.FrameRate = opts.fps; vw.Quality = 95; open(vw);
end
hl = 0; hf = 0; visited = false; held = false;
for k = 1:numel(tf)
    hl = heading(dL(:, k), hl); hf = heading(dF(:, k), hf);
    Rl = [cos(hl) -sin(hl); sin(hl) cos(hl)]; Rf = [cos(hf) -sin(hf); sin(hf) cos(hf)];
    pl = Rl*arrow + L(:, k); pf = Rf*arrow + F(:, k);
    set(leader, 'XData', pl(1, :), 'YData', pl(2, :));
    set(follower, 'XData', pf(1, :), 'YData', pf(2, :));
    s = sq(L(:, k), P.b); set(safeBox, 'XData', s(1, :), 'YData', s(2, :));
    s = sq(L(:, k), P.a); set(keepOut, 'XData', s(1, :), 'YData', s(2, :));
    set(trail, 'XData', F(1, 1:k), 'YData', F(2, 1:k));
    if inRegion(k) && tf(k) >= P.window(1) && tf(k) <= P.window(2)
        visited = true;
    end
    if visited
        set(region, 'FaceAlpha', .5); set(check, 'String', 'goal reached');
    end
    set(label, 'String', sprintf('t = %4.1f s   (reach the goal region within [%g, %g] s)', tf(k), P.window));
    drawnow;
    nRep = 1; if visited && ~held, nRep = round(.8*opts.fps); held = true; end   % pause briefly at the goal
    if isempty(opts.video)
        pause(nRep/opts.fps);
        if ~isvalid(fig), return; end
    else
        fr = getframe(fig); for j = 1:nRep, writeVideo(vw, fr); end
    end
end
if ~isempty(opts.video), close(vw); fprintf('Saved %s\n', opts.video); end
end
