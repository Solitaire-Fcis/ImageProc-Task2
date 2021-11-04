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
    if choice == 4
        mask = LaplacianSharp();
        res = LinearFilter(I, mask, 'cutoff');
    end
    if choice == 5
        fprintf('Enter Which Direction [1.Horizontal,2.Vertical]\n');
        choice = input('Enter Number>>>');
        if choice == 1
            choice = 'H';
        end
        if choice == 2
            choice = 'V';
        end
        [Sob1,Sob2] = Sobel(choice);
        res1 = LinearFilter(I, Sob1, 'absolute');
        res2 = LinearFilter(I, Sob2, 'absolute');
        res = res1+res2;
    end
    if choice == 6
        res = EdgeMagnit(I);
    end
    imshow(res)
end

% ---------- Main Linear Filteration Function ----------
function res = LinearFilter(I, Filter, Postproc)
    [rows, columns] = size(I);
    I = padarray(I,[1 1], 0);
    classOfImage = class(I);
    if isa(I,'uint8') || isa(I,'uint16') || isa(I,'int8') || isa(I,'int16') || isa(I,'logical')
        I = cast(I,'single');
        Filter = cast(Filter,'single');
    end
    res = convn(I,Filter,'same');
    if(~isa(res,classOfImage))
        res = cast(res, classOfImage);
    end
    if strcmpi(Postproc, 'cutoff')
        for i = 1:rows
            for j = 1:columns
                if res(i,j) > 255
                    res(i,j) = 255;
                end
                if res(i,j) < 0
                    res(i,j) = 0;
                end
            end
        end
    end
    if strcmpi(Postproc, 'absolute')
        for i = 1:rows
            for j = 1:columns
                if res(i,j) < 0
                    res(i,j) = abs(res(i,j));
                end
            end
        end
        for i = 1:rows
            for j = 1:columns
                if res(i,j) > 255
                    res(i,j) = 255;
                end
                if res(i,j) < 0
                    res(i,j) = 0;
                end
            end
        end
    end
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
    lap1 = [ 0,-1,0; -1,5,-1; 0,-1,0 ];
    lap2 = [ -1,-1,-1; -1,9,-1; -1,-1,-1 ];
    fprintf('Enter Which Filter To Use [1.(5)Centered,2.(9)Centered]\n');
    choice = input('Enter Number>>>');
    if choice == 1
        mask = lap1;
    end
    if choice == 2
        mask = lap2;
    end
end

% ---------- SOBEL ----------
function [Sob1, Sob2] = Sobel(maskType)
    if strcmpi(maskType, 'H')
        Sob1 = [ -1,-2,-1; 0,0,0; 1,2,1 ];
        Sob2 = [ 1,2,1; 0,0,0; -1,-2,-1 ];
    end
    if strcmpi(maskType, 'V')
        Sob1 = [ -1,0,1; -2,0,2; -1,0,1 ];
        Sob2 = [ 1,0,-1; 2,0,-2; 1,0,-1 ];
    end
end

% ---------- EDGE MAGNIT -----------
function goal = EdgeMagnit(I)
    [SobH1, SobH2] = Sobel('H');
    [SobV1, SobV2] = Sobel('V');
    resH1 = LinearFilter(I, SobH1, 'absolute');
    resH2 = LinearFilter(I, SobH2, 'absolute');
    resH = resH1+resH2;
    resV1 = LinearFilter(I, SobV1, 'absolute');
    resV2 = LinearFilter(I, SobV2, 'absolute');
    resV = resV1+resV2;
    goal = resH+resV;
end
