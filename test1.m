%reading image
name='ata.jpg';
img = double(imread(name))/255;

%setting parameters
output='m2.jpg';
sigma=5;
numberOfSectors=8;
q=5;
epsilon = 10^-4;
[size1, size2, N] = size(img);

%this part is used for creating gaussian kernel
[x,y] = meshgrid(linspace(-size2/2, size2/2, size2), linspace(size1/2, -size1/2, size1));
gaussianKernel = exp( -(x.^2+y.^2) / (2*sigma^2) );

%this part converts each channel of image to fourier domain.
for i = 1 : N
    imW(:,:,i)  = fft2(img(:,:,i));
    im2W(:,:,i) = fft2(img(:,:,i).^2);
end

%this part creates zero matrices
num = zeros(size1,size2,N);
den = zeros(size1,size2);

for i = 0 : numberOfSectors-1
    
%this part creates cutting functions and multiplies it with gaussian kernel
%this operation creates weighting functions
    G = smoothgaux(sector(size1,size2, i*2*pi/numberOfSectors, pi/numberOfSectors), 1, 2.5) .* gaussianKernel;
    G = G / sum(G(:));
    
%this part converts weighting functions to fourier domain
    G = fft2(G);
     
    S = zeros(size1,size2);
    
%this part calculates mean and standart deviation values
%this operation made in fourier domain and after that it converted back
    for k = 1 : N
        m(:,:,k) = ifft2(G .* imW(:,:,k));
        S = S +    ifft2(G .* im2W(:,:,k)) - m(:,:,k).^2;
    end

%this part applies same operation for each channel in color images
    S = (S+epsilon).^(-q/2);
    den = den + S;
    for k = 1 : N
        num(:,:,k) = num(:,:,k) + m(:,:,k).*S;
    end
   
end

%this part creates output image
for k = 1 : N
    y(:,:,k) = fftshift(num(:,:,k) ./ den);
end
figure, imshow(y);
imwrite(y, output);
