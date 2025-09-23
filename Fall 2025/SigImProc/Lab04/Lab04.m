%%
% Slade Brooks
% brooksl@mail.uc.edu
% 09.17.25
% BME6013C
% Lab 04

clear variables
close all

%% Part 1
% set image url
imageURL = "https://pressbooks.pub/app/uploads/sites/3987/2022/09/elbow-1.png";

% read in image
I = webread(imageURL);

% plot image to verify
figure(); imshow(I); axis image;
xlabel("x (px)"); ylabel("y (px)");
title("original image")

%% Part 2
% take max of third dimension (color data) to get grayscale
M = max(I, [], 3);

% plot to verify
figure(); imagesc(M); colormap("gray"); colorbar; axis image;
xlabel("x (px)"); ylabel("y (px)");
title("grayscale image")

%% Part 3
% filter image w/ medfilt (experimented to find good vals)
filt = medfilt2(M, [13 10]);

% the 15, 10 gets rid of the text for the most part, only leaving a large
% watermark and some ghost text that can be filtered out further. any
% larger filter steps would get rid of the segments between the bones

% plot to verify
figure(); imagesc(filt); colormap("gray"); colorbar; axis image;
xlabel("x (px)"); ylabel("y (px)");
title("medfilt image")

% the range shown in the plot is slightly lower than the grayscale plot in
% Part 2, which makes sense as we removed the text which was perfect white
% and the image we are left with should only have grays

%% Part 4
% first get the background out
background = (filt < 30);

% now find soft tissue and bone
bone = ~background & (filt > 135);
tissue = ~background & ~bone;

% create segmented image w/ correct vals
segmented = 0.*background + 1.*tissue + 2.*bone;

% plot segmented image
figure(); imagesc(segmented); colormap("gray"); axis image;
xlabel("x (px)"); ylabel("y (px)");
title("segmented image")

% these values do a decent job of segmenting. there are some issues since
% we are using only grayscale, as some values of tissue are higher than
% shadowy areas of bone and it isn't possible with any values to perfectly
% separate them because of this

%% Part 5
% first do some erosion to get rid of the text remnants and reflections
se = strel("line", 14, 45);
final = imerode(segmented, se);

% also erode some more to smooth the bumps on the borders of the tissue and
% expand the space between the two bones a little more
se = strel("line", 12, 90);
final = imerode(final, se);

% plot segmented image
figure(); imagesc(final); colormap("gray"); axis image;
xlabel("x (px)"); ylabel("y (px)");
title("final image")

% we can't use holes because it fills in gaps between bones and some of the
% gaps need to stay. we do have some holes in the bone that are actually
% just dark spots

% the final image is fairly well segmented. as mentioned before, the
% inability to use holes without filling in important features leads to a
% couple voids in the bone where there are shadows. as mentioned in part 4,
% there is some detail missing from the left side of the xray where there
% is a shadow and since the processing is in grayscale, there is no way to
% differentiate between highlighted tissue and shadowed bone. the erosion
% steps do a good job of removing the remaining text, smoothing edges, and
% re-expanding some of the desired spaces between bones to be more
% accurate. the erosion also helps address the reflections at the top and
% bottom. compared to the original image, the final segmented image matches
% the background well besides a couple blemishes, matches the tissue shape
% very closely, and does a decent job with the bone shaping. the top,
% right, and bottom of the bone are fairly well represented. the only major
% issue is the shadow on the left side of the xray which makes the bone not
% visible there. there are also a couple spots that show tissue within the
% bone in locations where there are shadows, but are unable to be repaired
% with imfill because it fills in gaps that should exist as well.