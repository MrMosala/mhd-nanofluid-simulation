# Research Videos

Place your research videos here.

## Recommended Content

1. **Flow Animation** - Visualization of velocity field
2. **Temperature Evolution** - Thermal field development
3. **Parameter Sensitivity** - Effect of Ha, Re, Pr, Ec
4. **Numerical Method Explanation** - SQLM walkthrough
5. **Entropy Generation** - Thermodynamic analysis

## Recommended Format

- **Format:** MP4 (H.264 codec)
- **Resolution:** 1080p or 720p
- **File Size:** Keep under 50MB per video for web optimization
- **Length:** 30 seconds to 3 minutes

## Tips for Creating Videos

### From MATLAB Animations
```matlab
% In MATLAB, export animations as video:
v = VideoWriter('flow_animation.mp4', 'MPEG-4');
v.FrameRate = 30;
open(v);
% ... your animation loop ...
writeVideo(v, getframe(gcf));
% ... end loop ...
close(v);
```

### From Screen Recording
- Use OBS Studio (free): https://obsproject.com
- Or Windows Game Bar: Win+G

## After Adding Videos

Update the `renderVideos` function in `src/App.js` to display your videos.

Example:
```jsx
<video 
  controls 
  width="100%"
  poster="/images/video_thumbnail.png"
>
  <source src="/videos/flow_animation.mp4" type="video/mp4" />
  Your browser does not support the video tag.
</video>
```

## Netlify Considerations

For large videos (>10MB), consider hosting on:
- YouTube (embed using iframe)
- Vimeo (embed using iframe)
- Google Drive (share link)

This reduces your Netlify build size and improves loading speed.

### YouTube Embed Example
```jsx
<iframe
  width="100%"
  height="100%"
  src="https://www.youtube.com/embed/YOUR_VIDEO_ID"
  frameBorder="0"
  allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
  allowFullScreen
/>
```
