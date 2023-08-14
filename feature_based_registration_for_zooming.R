library(DRIP)
library(matlab)
library(jpeg)


w1<-30


img<- brain     ##Reference image
img_zoom<- brain_z   ##Zoomed image




      img1<-matrix(NA,nrow=nrow(img),ncol=ncol(img))
      img_zoom1<-matrix(NA,nrow=nrow(img_zoom),ncol=ncol(img_zoom))
      
      img1<-img + matrix(rnorm(nrow(img1)*ncol(img1),0,0.01),nrow = nrow(img1),ncol = ncol(img1))   ##reference image (w-noise)
      img_zoom1<-img_zoom + matrix(rnorm(nrow(img_zoom1)*ncol(img_zoom1),0,0.01),nrow = nrow(img_zoom1),ncol = ncol(img_zoom1)) ##Zoomed image ( w-noise)
      
      edge_ref<-which(stepEdgeLC2K(image=img1,bandwidth=2,thresh=0.1,plot=FALSE)==1,arr.ind = TRUE)
      edge_zoomed<-which(stepEdgeLC2K(image=img_zoom1,bandwidth=2,thresh=0.1,plot=FALSE)==1,arr.ind = TRUE)
      
      
      mat1<- padarray(img_zoom1,c(w1,w1),"symmetric","both") ##Zoomed image (padded)
      mat2<- padarray(img1,c(w1,w1),"symmetric","both")   ##Reference image (padded)
      
      
      edg<-edge_zoomed+w1
      edge_reg_L1<-matrix(NA,nrow = nrow(edg),ncol = 2)
      edge_reg_L2<-matrix(NA,nrow = nrow(edg),ncol = 2)
      edge_reg_cor<-matrix(NA,nrow = nrow(edg),ncol = 2)
      
      edge_x<-edg[,1]
      edge_y<-edg[,2]
      
      for(i in 1:nrow(edg))
      {
        
        N1<-c()
        for(k1 in (edge_x[i]-r1):(edge_x[i]+r1))
        {
          for(l1 in (edge_y[i]-r1):(edge_y[i]+r1))
          {
            if( (k1-edge_x[i])^2 + (l1-edge_y[i])^2 <= (r1)^2)
            {
              N1<-c(N1,mat1[k1,l1])
            }
          }
        }
        msd<-c()
        mad<-c()
        
        xz<-c()
        yz<-c()
        
        for(m1 in (edge_x[i]-r2):(edge_x[i]+ r2))
        {
          for(n1 in (edge_y[i]-r2):(edge_y[i]+ r2))
          {
            if( ((m1-edge_x[i])^2 + (n1-edge_y[i])^2) <= (r2)^2)
            {
              N2<-c()
              N<-0
              for(s1 in (m1-r1):(m1+r1))                              
              {
                for(t1 in (n1-r1):(n1+r1))
                {
                  if( ((s1-m1)^2+(t1-n1)^2) <= (r1)^2)
                  {
                    N2<-c(N2,mat2[s1,t1])
                    N<-N+1
                  }
                }
              }
              
              msd<-c(msd,sum((N1-N2)^2)/N)
              mad<-c(mad,sum(abs(N1-N2))/N)
              xz<-c(xz,round(m1))
              yz<-c(yz,round(n1))
            }
          }
        }
        index_L2<-which.min(msd)
        index_L1<-which.min(mad)
        
        edge_reg_L1[i,1] <- xz[index_L1]
        edge_reg_L1[i,2] <- yz[index_L1]
        edge_reg_L2[i,1] <- xz[index_L2]
        edge_reg_L2[i,2] <- yz[index_L2]
        
        
        
      }
      edge_x_L1<-c()
      edge_x_L1<-edge_reg_L1[,1]
      edge_y_L1<-c()
      edge_y_L1<-edge_reg_L1[,2]
      
      edge_x_L2<-c()
      edge_x_L2<-edge_reg_L2[,1]
      edge_y_L2<-c()
      edge_y_L2<-edge_reg_L2[,2]
      
      # edge_x_cor<-edge_reg_cor[,1]
      # edge_y_cor<-edge_reg_cor[,2]
      
      s_L1<-0
      h_L1<-0
      m_L1<-0
      
      s_L2<-0
      h_L2<-0
      m_L2<-0
      
      mod_L1_x<-list()
      mod_L1_x<-lm(edge_x~edge_x_L1)
      mod_L1_y<-list()
      mod_L1_y<-lm(edge_y~edge_y_L1)
      s_L1= round(mod_L1_x$coefficients[2],2)
      h_L1= mod_L1_x$coefficients[1]/(1-s_L1)
      m_L1=mod_L1_y$coefficients[1]/(1-s_L1)
      reg_L1<- matrix(NA, nrow = nrow(mat1),ncol = ncol(mat1))
      for(i in 1:(nrow(mat1)))
      {
        for(j in 1:(nrow(mat1)))
        {
          i1<- ceiling((i-h_L1)/s_L1 +h_L1)
          j1<- ceiling((j-m_L1)/s_L1 +m_L1)
          if(i1<1)
            i1<-1
          else if(i1>nrow(mat1))
            i1<-nrow(mat1)
          
          if(j1<1)
            j1<-1
          else if(j1>nrow(mat1))
            j1<-nrow(mat1)
          reg_L1[i1,j1]<-mat1[i,j]
        }
      }
      
      mod_L2_x<-list()
      mod_L2_x<-lm(edge_x~edge_x_L2)
      mod_L2_y<-list()
      mod_L2_y<-lm(edge_y~edge_y_L2)
      s_L2= round(mod_L2_x$coefficients[2],2)
      h_L2= mod_L2_x$coefficients[1]/(1-s_L2)
      m_L2=mod_L2_y$coefficients[1]/(1-s_L2)
      reg_L2<- matrix(NA, nrow = nrow(mat1),ncol = ncol(mat1))
      for(i in 1:(nrow(mat1)))
      {
        for(j in 1:(nrow(mat1)))
        {
          i2<- ceiling((i-h_L2)/s_L2 +h_L2)
          j2<- ceiling((j-m_L2)/s_L2 +m_L2)
          if(i2<1)
            i2<-1
          else if(i2>nrow(mat1))
            i2<-nrow(mat1)
          
          if(j2<1)
            j2<-1
          else if(j2>nrow(mat1))
            j2<-nrow(mat1)
          reg_L2[i2,j2]<-mat1[i,j]
        }
      }
      
      
     
      
      t1<-0
      t2<-0
      # t3<-0
      n1<-0
      n2<-0
      # n3<-0
      for(i in (w1+r2):(nrow(mat2)-w1-r2+1))
      {
        for(j in (w1+r2):(nrow(mat2)-w1-r2+1))
        {
          if(is.na(reg_L1[i,j])==FALSE)
          {
            t1<-t1+(mat2[i,j]-reg_L1[i,j])^2
            n1<-n1+1
          }
          
          if(is.na(reg_L2[i,j])==FALSE)
          {
            t2<-t2+(mat2[i,j]-reg_L2[i,j])^2
            n2<-n2+1
            
          }
          
          # if(is.na(reg_cor[i,j])==FALSE)
          # {
          #   t3<-t3+(mat2[i,j]-reg_cor[i,j])^2
          #   n3<-n3+1
          # }
          
        }
      }
      
      msd_L1<- t1/n1
      msd_L2<-t2/n2
      
      L1_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
      L1_img[(r2+1):((nrow(L1_img)-r2)),(r2+1):((ncol(L1_img)-r2))]<-reg_L1[(w1+r2+1):((nrow(img)-w1-r2)),(w1+r2+1):((ncol(img)-w1-r2))]  ##Registered image under L1-norm
      
      L2_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
      L2_img[(r2+1):((nrow(L2_img)-r2)),(r2+1):((ncol(L2_img)-r2))]<-reg_L2[(w1+r2+1):((nrow(img)-w1-r2)),(w1+r2+1):((ncol(img)-w1-r2))]  ##Registered image under L2-norm
      
      ##Residual image
      
      res_L1<- img - L1_img  ##residual image under L1-norm
      
      res_L2<- img - L2_img  ##residual image under L2-norm
      
      ##Visualization
      image(rot90(L1_img,3),col = grey(seq(0,1,length=256)))
      