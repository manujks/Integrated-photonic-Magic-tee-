% Create and save warp function
function createandsave_warpfunction2(wstart,wend,zvector,z,fv)
[xxx,www]=getwidthwarp(wstart, wend, zvector, [z(:) fv(:)]);
figure; plot(xxx.',[+www.'/2 -www.'/2]);
warpfcn=[z(:) fv(:)];
save('warpfunctionxy_lamXXXXnm.txt','warpfcn','-ASCII');    % For Jason

savefimmwarpfcn('Sitaper_lamXXXXnm_warpfcn.txt',z,fv);      % For FIMMWAVE
