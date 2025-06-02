function gausbg=FPN(dc,hsize,sigma)
[m,n,b]=size(dc);
scale=6;
gm1=cell(scale,1);
gm2=cell(scale+1,1);
gausFilter = fspecial('gaussian', [hsize,hsize], sigma);
for j=1:b                                %constructing hyperspectral feature pyramid 
    gm2{1}(:,:,j)=dc(:,:,j);
    for i=1:scale
        gm1{i}(:,:,j)=gm2{i}(:,:,j);
        gm1{i}(:,:,j)= gm1{i}(:,:,j)-imfilter(gm1{i}(:,:,j), gausFilter, 'replicate');
        gm2{i+1}(:,:,j)= imresize(gm1{i}(:,:,j), ceil(size(gm1{i}(:,:,j))/2), 'bilinear');
    end
end


gas=cell(scale,scale);                     %cross-scale subtraction
for c=1:3
    for delta = 1:3
        s = c + delta;
        for j=1:b
            gas{c,s}(:,:,j) = abs(Subtract(gm1{c}(:,:,j), gm1{s}(:,:,j)));
        end
    end
end

gausbg=zeros(m,n,b);
for c=1:3
    for delta = 1:3
        s = c + delta;
        for j=1:b
            gausbg(:,:,j) = Add(gausbg(:,:,j),gas{c,s}(:,:,j));   %cross-scale summation
        end
    end
end

function result = Subtract(im1, im2)
    im2 = imresize(im2, size(im1), 'bilinear');
    result = im1 - im2;
end

function result = Add(im1, im2)
    im1 = imresize(im1, [m,n], 'bilinear');
    im2 = imresize(im2, [m,n], 'bilinear');
    result = im1 + im2;
end

end