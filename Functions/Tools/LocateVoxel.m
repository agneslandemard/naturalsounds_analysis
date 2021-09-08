function vx = LocateVoxel(type,infos,D)
%type can be 'choose' or 'locate'
global data_path

switch type
    case 'choose'
        X = load([data_path infos.hemi D.data_suffix '.mat'],'param','Anat');
       
        j = find(strcmp(D.hemis,infos.hemi));
        
        figure;
        imagesc(X.Anat(:,:,infos.slice));
        axis equal tight
        colormap hot
        
        % Select voxel in chosen hemi and slice
        px = round(ginput(1));
        PXID = X.Anat.*0;
        PXID(px(2),px(1),infos.slice) = 1;
        pix = find(Mat2Pixs(PXID,X.param.msk));
        vx = find(D.si==j,1)+pix-1;
        
    case 'locate'
        pixf = infos.voxel;
        
        vx.hemi = D.hemis{D.si(pixf)};
        X = load([data_path vx.hemi D.data_suffix '.mat'],'param','Anat');
       

        vx_in_hemi = pixf-find(D.si==D.si(pixf),1)+1;
        
        PXID = zeros(length(find(D.si==D.si(pixf))),1);
        PXID(vx_in_hemi) = 1;
        pix = find(Pixs2Mat(PXID,X.param.msk)==1);
        [vx.x,vx.y,vx.slice] = ind2sub(size(Pixs2Mat(PXID,X.param.msk)),pix);
      
        if infos.plot
            
            figure;
            imagesc(sqrt(X.Anat(:,:,vx.slice)));
            colormap gray
            hold on
            scatter(vx.y,vx.x,'*w')
            axis equal tight
            title([vx.hemi ', slice ' num2str(vx.slice)])

        end

end