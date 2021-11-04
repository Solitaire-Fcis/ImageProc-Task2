% ---------- MAIN ----------
function Task2IP()
    I = imread('I.bmp'); 
    fprintf("Enter Number of Desired Filter:::\n\t1.Mean Filter\n\t2.Gaussian 1 (Size,Sigma) Filter\n\t3.Gaussian 2 (Sigma) Filter\n\t4.Laplacian Sharpness\n\t5.Sobel\n\t6.Edge Magnit\n");
    choice = input('Enter Number>>>');
    if choice == 1
        W = input('Enter Width: ');
        H = input('Enter Height: ');
        mask = MeanMask(W,H);
        res = LinearFilter(I, mask, 'none');
    end
    if choice == 2
        size = input('Enter Size: ');
        sig = input('Enter Sigma: ');
        mask = Gauss1(size,sig);
        res = LinearFilter(I, mask, 'none');
    end
    if choice == 3
        sig = input('Enter Sigma: ');
        mask = Gauss2(sig);
        res = LinearFilter(I, mask, 'none');
    end
    imshow(res)
end

% ---------- Main Linear Filteration Function ----------
function res = LinearFilter(I, Filter, Postproc)
    if Postproc == 'none'
        res = imfilter(I, Filter, 'conv');
    end
    %if Postproc == 'cutoff'
    %end
    %if Postproc == 'absolute'
    %end
    
end
% ----------   MEAN FILTER -----------
function mask = MeanMask(W,H)
    mask = ones(W,H);
    mask = (1/(W*H))*mask;
end
function mask = Gauss1(size,sig)
    halfsize = cast(size/2, 'int32')-1;
    mask = zeros(size,size);
    x = (-halfsize);
    y = (-halfsize);
    if mod(size,2) == 1
        x = x - 1;
        y = y - 1;
    end
    for i = 1:size
        tmp = y;
        for j = 1:size
            p1 = (1/(2*pi*sig*sig));
            p2 = (exp(double(-(x*x + y*y)/(2*sig*sig))));
            y = y+1;
            mask(i,j) = double(p1*p2);
        end
        y = tmp;
        x = x+1;
    end
    mask = mask/sum(mask(:));
end
% ---------- Gaussian2 Filter ---------
function mask = Gauss2(sig)
    N = cast((3.7*sig) - 0.5, 'uint32');
    maskSize = (2*N)+1;
    halfsize = cast(maskSize/2, 'int32');
    mask = zeros(maskSize,maskSize);
    x = (-halfsize);
    y = (-halfsize);
    if mod(maskSize,2) == 1
        x = x - 1;
        y = y - 1;
    end
    for i = 1:maskSize
        tmp = y;
        for j = 1:maskSize
            p1 = (1/(2*pi*sig*sig));
            p2 = (exp(double(-(x*x + y*y)/(2*sig*sig))));
            y = y+1;
            mask(i,j) = double(p1*p2);
        end
        y = tmp;
        x = x+1;
    end
    mask = mask/sum(mask(:));
end
% ---------- Laplacian Sharpnes -------
function mask = LaplacianSharp()

end

% ---------- SOBEL ----------
function masks = Sobel(mask)

end

% ---------- EDGE MAGNIT -----------
function goal = EdgeMagnit(I)

end
