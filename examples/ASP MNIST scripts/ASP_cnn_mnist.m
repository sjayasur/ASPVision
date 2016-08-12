function [net, info] = cnn_mnist(varargin)
% CNN_MNIST  Demonstrated MatConNet on MNIST

run(fullfile(fileparts(mfilename('fullpath')),...
  '..', 'matlab', 'vl_setupnn.m')) ;
addpath(‘..’)
opts.expDir = fullfile('data','mnist-ASP') ;
[opts, varargin] = vl_argparse(opts, varargin) ;

opts.dataDir = fullfile('data','mnist') ;
opts.imdbPath = fullfile(opts.expDir, 'imdb.mat');
opts.useBnorm = false ;
opts.train.batchSize = 100 ;
opts.train.numEpochs = 20 ;
opts.train.continue = true ;
opts.train.gpus = 1 ;
opts.train.learningRate = 0.001 ;
opts.train.expDir = opts.expDir ;
opts = vl_argparse(opts, varargin) ;

% --------------------------------------------------------------------
%                                                         Prepare data
% --------------------------------------------------------------------

if exist(opts.imdbPath, 'file')
  imdb = load(opts.imdbPath) ;
else
  imdb = getMnistImdb(opts) ;
  mkdir(opts.expDir) ;
  save(opts.imdbPath, '-struct', 'imdb') ;
end

net = ASP_cnn_mnist_init('useBnorm', opts.useBnorm) ;

% --------------------------------------------------------------------
%                                                                Train
% --------------------------------------------------------------------

[net, info] = cnn_train(net, imdb, @getBatch, ...
    opts.train, ...
    'val', find(imdb.images.set == 3)) ;

% --------------------------------------------------------------------
function [im, labels] = getBatch(imdb, batch)
% --------------------------------------------------------------------
im = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.labels(1,batch) ;

% --------------------------------------------------------------------
function imdb = getMnistImdb(opts)
% --------------------------------------------------------------------
% Preapre the imdb structure, returns image data with mean image subtracted
files = {'train-images-idx3-ubyte', ...
         'train-labels-idx1-ubyte', ...
         't10k-images-idx3-ubyte', ...
         't10k-labels-idx1-ubyte'} ;

if ~exist(opts.dataDir, 'dir')
  mkdir(opts.dataDir) ;
end

for i=1:4
  if ~exist(fullfile(opts.dataDir, files{i}), 'file')
    url = sprintf('http://yann.lecun.com/exdb/mnist/%s.gz',files{i}) ;
    fprintf('downloading %s\n', url) ;
    gunzip(url, opts.dataDir) ;
  end
end

f=fopen(fullfile(opts.dataDir, 'train-images-idx3-ubyte'),'r') ;
x1=fread(f,inf,'uint8');
fclose(f) ;
x1=permute(reshape(x1(17:end),28,28,60e3),[2 1 3]) ;

f=fopen(fullfile(opts.dataDir, 't10k-images-idx3-ubyte'),'r') ;
x2=fread(f,inf,'uint8');
fclose(f) ;
x2=permute(reshape(x2(17:end),28,28,10e3),[2 1 3]) ;

f=fopen(fullfile(opts.dataDir, 'train-labels-idx1-ubyte'),'r') ;
y1=fread(f,inf,'uint8');
fclose(f) ;
y1=double(y1(9:end)')+1 ;

f=fopen(fullfile(opts.dataDir, 't10k-labels-idx1-ubyte'),'r') ;
y2=fread(f,inf,'uint8');
fclose(f) ;
y2=double(y2(9:end)')+1 ;

set = [ones(1,numel(y1)) 3*ones(1,numel(y2))];
data = single(reshape(cat(3, x1, x2),28,28,1,[]));







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data [32 32 3 60000] --> ASP data [32 32 25 60000]
%load parameters
loadASPparam;

%Tile pattern

% no sin/cosine distinction
nfSort1=[47 35  5 43 13 19 ...
          1 29 39 27 ...
           9 21];
nfSort2= nfSort1 + 1;


% Construct ASP transform matrix, stored in ASP variable. GW_sj = makes
% gabor wavelet corresponding to each type of ASP 

M = 7;   % n samples in response (long axis)
A = 0.7; % peak amplitude of response
S = 2;   % Scale of gaussian window
I0 = 0.3; % base intensity scale


for zz=1:12
    lamb = nfSort1(zz);
    lamb2 = nfSort2(zz);
    temp = Gw(I0,M,A,par(lamb,:));   
    temp2 = Gw(I0,M,A,par(lamb2,:));
    ASP{zz} = temp-temp2;
    ASP{zz} = fliplr(ASP{zz});
    ASPtemp = ASP{zz};
    ASP{zz} = (ASPtemp-min(ASPtemp(:)))./(max(ASPtemp(:))-min(ASPtemp(:))); %add normalization
    ASP{zz} = ASP{zz} .* 2 - 1;
%     figure, imagesc(ASP{zz}); colormap gray;
end
% temp = Amp(M,S);
% ASP{25} = temp;
% ASP{25} = (temp-min(temp(:)))./(max(temp(:))-min(temp(:)));

wscale = 1;


for zz=1:12
    for c=1:1
        wASP(:,:,c,zz) = wscale.*ASP{zz};
    end
end

wASP = single(wASP);


% Perform ASP convolutions
data = vl_nnconv(data,wASP,[],'Pad',3);  
%The pad is added when the data size is less than 32x32
data = (data-min(data(:)))./(max(data(:))-min(data(:))).*255;
data = single(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % DIMENSIONALITY CHECK
% for c=4:25
% data(:,:,c,:) = squeeze(randn(32,32,size(data,4)));
% end

% remove mean in any case
dataMean = mean(data(:,:,:,set == 1), 4);
data = bsxfun(@minus, data, dataMean);


imdb.images.data = data ;
imdb.images.data_mean = dataMean;
imdb.images.labels = cat(2, y1, y2) ;
imdb.images.set = set ;
imdb.meta.sets = {'train', 'val', 'test'} ;
imdb.meta.classes = arrayfun(@(x)sprintf('%d',x),0:9,'uniformoutput',false) ;
