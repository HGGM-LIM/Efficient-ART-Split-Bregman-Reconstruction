function h = Plot2DMapsGridSolution(ValGrid,X,Y,Z,dimPlot)
% h = Plot2DMapsGridSolution(ValGrid,X,Y,Z,dimPlot)
%
% Plot 3D scalar function that corresponds to X,Y,Z from meshgrid
%
% ValGrid     = image, nx-ny-nz
% dimPlot     = 1,2 or 3
%
% Code downloaded from the repository
% https://github.com/HGGM-LIM/Efficient-ART-Split-Bregman-Reconstruction
%
% Juan FPJ Abascal, Judit Chamorro-Servent, Juan Aguirre
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es

switch dimPlot
    case 3
        xlin        = unique(X(:,:,1));
        ylin        = unique(Y(:,:,1));
    case 2
        xlin        = unique(X(:,:,1));
        ylin        = unique(Z(:,1,:));
    case 1
        xlin        = unique(Y(:,:,1));
        ylin        = unique(Z(1,:,:));
end % switch


valMax      = full(max(ValGrid(:)));
valMin      = full(min(ValGrid(:)));
numSubplot  = ceil(sqrt(size(ValGrid,dimPlot)));
%
switch dimPlot
    case 3
        numPlots = size(ValGrid,3);
    case 2
        numPlots = size(ValGrid,2);
    case 1
        numPlots = size(ValGrid,1);
end % switch
%
h = figure;
for iq = 1:numPlots
    subplot(numSubplot,numSubplot,iq);
    
    switch dimPlot
        case 3
            imagesc(xlin,ylin,squeeze(ValGrid(:,:,iq)));
            title(['z: ' num2str(unique(Z(:,:,iq)))]);
            axis image;
        case 2
            imagesc(ylin,xlin,squeeze(ValGrid(:,iq,:)));
            title(['y: ' num2str(unique(Y(iq,:,:)))]);
            axis image;
        case 1
            imagesc(ylin,xlin,squeeze(ValGrid(iq,:,:)));
            title(['x: ' num2str(unique(X(:,iq,:)))]);
            axis image;
    end % switch     
    if (valMin == valMax)
        valMin = 0;
    else
         caxis([valMin valMax]);    
    end
    %
end % iq
