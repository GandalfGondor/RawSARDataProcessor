function image = multilookProcessing(slcimg,sx,sy)
[nx,ny] = size(slcimg);
nfx = floor(nx/sx);
nfy = floor(ny/sy);
image = (zeros(nfx,nfy));
for i=1:nfx
    for j = 1:nfy
        fimg=0;
        for ix = 1:sx
            for jy = 1:sy
                fimg = fimg+slcimg(((i-1)*sx)+ix,((j-1)*sy)+jy);
            end
        end
        image(i,j) = fimg/(sx*sy);
    end
end
end