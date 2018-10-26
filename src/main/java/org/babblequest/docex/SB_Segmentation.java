package org.babblequest.docex;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/*				sbseg.c  by Tom Goldstein
 *   This code performs isotropic segmentation using the "Split Bregman" algorithm.
 * An image of dimensions mxn is denoised using the following command
 * 
 *   sbseg(f,edge,mu);
 * 
 * where:
 * 
 *   - "f" is the mxn noisy image to be segmented
 *   - "mu" is the weighting parameter for the fidelity term (mu should be about 1e-4 for images with pixels on the 0-255 scale)
 *   - "edge" is a 2D array of edge detector values at each pixel (e.g. the weights for TV regularizer).  For standard segmentation
 *         simply choose "edge = ones(size(f))"
 */

public class SB_Segmentation implements PlugInFilter
{

 private int flags = DOES_ALL;
   int GX[][] = { {-1, 0, 1}, {-2,0,2},{-1,0,1}};
   //GX[0][0] = -1; GX[0][1] = 0; GX[0][2] = 1;
   //GX[1][0] = -2; GX[1][1] = 0; GX[1][2] = 2;
   //GX[2][0] = -1; GX[2][1] = 0; GX[2][2] = 1;

   /* 3x3 GY Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
   int GY[][] = { {1,2,1},{0,0,0},{-1,-2,-1}};

  public int setup(String argv, ImagePlus imp)
  {
      if (IJ.versionLessThan("1.38x"))        // generates an error message for older versions
            return DONE;
        return flags;
  }

  public void run(ImageProcessor ip) {
    int pixels[][] = ip.getIntArray();
    int cols = ip.getWidth();
    int rows = ip.getHeight();

    double[][] f = newMatrixseg(rows,cols);

    for (int r=0;r<rows;r++)
       for (int c=0;c<cols;c++)
       {
          double val = (double)pixels[c][r];
          f[r][c] = val;
       }

    double[][] edges = newMatrixseg(rows,cols);
    edges = getEdges(pixels, cols, rows);

    for (int r=0;r<rows;r++)
       for (int c=0;c<cols;c++)
       {
          edges[r][c] = 1.0;
          //edges[r][c] = 1.0 - edges[r][c];
       }

    double[][] u = segment(f, edges, rows, cols);

    /* copy denoised image to output vector*/
        for(int r=0;r<rows;r++)
            for(int c=0;c<cols;c++){
                if (u[r][c] > 0.5)
                  ip.setf(c,r,(float)255.0);
                else
                  ip.setf(c,r,(float)0.0);
                //ip.setf(c,r,(float)u[r][c]);  // u is the denoised image
            }
    
  }

public double[][] getEdges(int[][] pixels, int width, int height)
{
   double edges[][] = new double[height][width];

   int xMin=1;
   int xMax=width-2;
   int yMin=1;
   int yMax=height-2;

   for (int y=yMin; y<=yMax; y++)
   {
        for (int x=xMax; x<=xMax; x++)
        {
            int sumX = 0;
            int sumY = 0;

            for(int I=-1; I<=1; I++)  {
                   for(int J=-1; J<=1; J++)  {
                      sumX = sumX + (pixels[x+I][y+J] * GX[I+1][J+1]);
                   }
               }

               /*-------Y GRADIENT APPROXIMATION-------*/
               for(int I=-1; I<=1; I++)  {
                   for(int J=-1; J<=1; J++)  {
                       sumY = sumY + (pixels[x+I][y+J] * GY[I+1][J+1]);
                   }
               }

               /*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
               //double sum = Math.abs(sumX) + Math.abs(sumY);
               double sum = Math.abs(sumX) + Math.abs(sumY);
/*
               double sum = 0.0;

               if (Math.abs(sumX) < Math.abs(sumY))
                  sum = Math.abs(sumY);
               else
                  sum = Math.abs(sumX);

*/
            sum = 1.0/(1.0+sum*sum);
/*
            if(sum>40) sum=0;
            else sum=1;
*/

System.out.println(sum);
            edges[y][x] = sum;
        }
     }

     return(edges);
}
	
public double[][] segment(double[][] f, double[][] edges,  int rows, int cols)
{
        /* get the fidelity and convergence parameters*/
	double mu =  0.001;
	double lambda = 0.3;
	/*double tol = (double)(mxGetScalar(prhs[3]));*/
	
	double[][] u = newMatrixseg(rows,cols);
	double[][] x = newMatrixseg(rows-1,cols);
	double[][] y = newMatrixseg(rows,cols-1);
	double[][] bx = newMatrixseg(rows-1,cols);
	double[][] by = newMatrixseg(rows,cols-1);
	
	double diff;
	int count=0;
	double c1,c2;
        double uOld[][];
        double change;
    
        int i,j;
    
	/* segment the image*/

	uOld = newMatrixseg(rows,cols);
	double rtn[] = guessC(f,u,rows,cols);
        c1 = rtn[0];
        c2 = rtn[1];
    
    
	do{
		rtn = sbseg(f,edges,u,x,y,bx,by,c1,c2,mu,lambda,5,rows,cols);
                c1 = rtn[0];
                c2 = rtn[1];
                change = rtn[2];

		diff = copyseg(u,uOld,rows,cols);
	  /*}while( count++<500 && diff>tol);*/
     }while( count++<500 && change>.025);

    return(u);
}

public void fillArray(double[][] a, double val, int rows, int cols){
	int r,c;
	for(r=0;r<rows;r++)
		for(c=0;c<cols;c++)
			a[r][c] = val;
	return;
}





/*                IMPLEMENTATION BELOW THIS LINE                         */

/******************Isotropic Segmentation**************/

public double[] sbseg(double[][] f, double[][] edge,  double[][] u, double[][] x,double[][] y, 
                    double[][] bx, double[][] by,
                    double c1, double c2, double mu, double lambda, int nIter, int width, int height)
{
	int j,n;
    double change=0;
	for(j=0;j<nIter;j++){	
        for(n=0;n<1;n++){
			change = gsUseg(f,u,x,y,bx,by,c1,c2,mu,lambda, width, height);
			gsSpaceseg(u,edge,x,y,bx,by,lambda, width, height);
            bregmanXseg(x, u, bx, width, height);
			bregmanYseg(y, u, by, width, height);
        }
			
			double rtn[] = updateC(f, u, .5, width, height);
                        c1 = rtn[0];
                        c2 = rtn[1];
		}

    double[] rtn = new double[3];

    rtn[0] = c1;
    rtn[1] = c2;
    rtn[2] = change;

    return (rtn);
}

public double[] updateC(double[][] f, double[][] u, double thresh, int rows, int cols){
	int r,c;
	double sum1=0, sum2=0;
	int n1=0,n2=0;
	for(r=0;r<rows;r++){
		for(c=0;c<cols;c++){
			if(u[r][c]>thresh){
				sum1+=f[r][c];
				n1++;
			}else{
				sum2+=f[r][c];
				n2++;
			}
		}
	}
	
	if(n1<1)n1=1;
	if(n2<1)n2=1;

        double[] rtn = new double[2];
	rtn[0] = sum1/(double)n1;
	rtn[1] = sum2/(double)n2;

	return(rtn);
}

public double[] guessC(double[][] f, double[][] u, int rows, int cols){
		int r,c;
		double av=0,sum1=0, sum2=0;
		int n1=0,n2=0;

		for(r=0;r<rows;r++){
			for(c=0;c<cols;c++){
				av+=f[r][c];
			}
		}
		av/=rows*cols;
		for(r=0;r<rows;r++){
			for(c=0;c<cols;c++){
				if(f[r][c]>av){
					sum1+=f[r][c];
					n1++;
					u[r][c]=1;
				}else{
					sum2+=f[r][c];
					n2++;
				}
			}
		}
		if(n1<1)n1=1;
		if(n2<1)n2=1;

                double[] rtn = new double[2];
		rtn[0] = sum1/(double)n1;
		rtn[1] = sum2/(double)n2;

		return (rtn);
	}

/****Relaxation operators****/

public double gsUseg(double[][] f, double[][] u, double[][] x, double[][] y, double[][] bx, double[][] by, double c1, double c2, double mu, double lambda, int rows, int cols){
		int r,c;
		double g,sum,temp1,temp2;
		  /* optimization varables*/
		int cm1;
        double old, diff, change=0;
		for(r=1;r<rows-1;r++){
			for(c=1;c<cols-1;c++){
                old = u[r][c];
				cm1 = c-1;
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-x[r][c]-bx[r-1][c]+bx[r][c]+y[r][cm1]-y[r][c]-by[r][cm1]+
                                        by[r][c]-mu/lambda*g;
				u[r][c] = 0.25*(u[r-1][c]+u[r+1][c]+u[r][c+1]+u[r][cm1]+sum);
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
                               diff = u[r][c]-old; diff*=diff; if(diff>change) change=diff;
			}
		}
		r=0;
			for(c=1;c<cols-1;c++){
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]+y[r][c-1]-y[r][c]-by[r][c-1]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c+1]+u[r][c-1]+sum)/3.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		r = rows-1;
			for(c=1;c<cols-1;c++){
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]+y[r][c-1]-y[r][c]-by[r][c-1]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c+1]+u[r][c-1]+sum)/3.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
			
		c = 0;
			for(r=1;r<rows-1;r++){
						temp1 = c1-f[r][c];
						temp2 = c2-f[r][c];
						g = temp1*temp1-temp2*temp2;
						sum = x[r-1][c]-x[r][c]-bx[r-1][c]+bx[r][c]-y[r][c]+by[r][c]-mu/lambda*g;
						u[r][c] = 0.25*(u[r-1][c]+u[r+1][c]+u[r][c+1]+sum);
						if(u[r][c]>1) u[r][c]=1;
						else if(u[r][c]<0) u[r][c]=0;
					}

		c = cols-1;
			for(r=1;r<rows-1;r++){
					temp1 = c1-f[r][c];
					temp2 = c2-f[r][c];
					g = temp1*temp1-temp2*temp2;
					sum = x[r-1][c]-x[r][c]-bx[r-1][c]+bx[r][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
					u[r][c] = (u[r-1][c]+u[r+1][c]+u[r][c-1]+sum)/3.0;
					if(u[r][c]>1) u[r][c]=1;
					else if(u[r][c]<0) u[r][c]=0;
				}
			
		r=c=0;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]-y[r][c]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c+1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=0;c=cols-1;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c-1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=rows-1;c=0;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]-y[r][c]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c+1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=rows-1;c=cols-1;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c-1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		return change;
}


public void gsSpaceseg(double[][] u,double[][] edge, double[][] x, double[][] y, double[][] bx, double[][] by, double lambda, int width, int height){
	int w,h;
	double a,b,s;
	double flux = 1.0/lambda;
    double base;
	for(w=0;w<width-1;w++){	
		for(h=0;h<height-1;h++){
			flux = edge[w][h]/lambda;
			a =  u[w+1][h]-u[w][h]+bx[w][h];
			b =  u[w][h+1]-u[w][h]+by[w][h];
			s = a*a+b*b;
			if(s<flux*flux){x[w][h]=0;y[w][h]=0;continue;}
			s = Math.sqrt(s);
			s=(s-flux)/s;
			x[w][h] = s*a;
			y[w][h] = s*b;
		}
	}		
	
	h = height-1;
	for(w=0;w<width-1;w++){	
			flux = edge[w][h]/lambda;
			base =  u[w+1][h]-u[w][h]+bx[w][h];
			if(base>flux) {x[w][h] = base-flux; continue;}
			if(base<-flux){x[w][h] = base+flux; continue;}
			x[w][h] = 0;
	}
	w = width-1;
	for(h=0;h<height-1;h++){	
		flux = edge[w][h]/lambda;
		base =  u[w][h+1]-u[w][h]+by[w][h];
		if(base>flux) {y[w][h] = base-flux; continue;}
		if(base<-flux){y[w][h] = base+flux; continue;}
		y[w][h] = 0;
	}
}



public void bregmanXseg(double[][] x,double[][] u, double[][] bx, int width, int height){
		int w,h;
		double d;
		for(w=0;w<width-1;w++){
			for(h=0;h<height;h++){
				d = u[w+1][h]-u[w][h];
				bx[w][h]+= d-x[w][h];		
			}
		}
	}


public void bregmanYseg(double[][] y,double[][] u, double[][] by, int width, int height){
		int w,h;
		double d;
		int hSent = height-1;
		for(w=0;w<width;w++){
			for(h=0;h<hSent;h++){
				d = u[w][h+1]-u[w][h];
				by[w][h]+= d-y[w][h];
			}
		}
	}
	
/************************memory****************/

public double copyseg(double[][] source, double[][] dest, int rows, int cols){
	int r,c;
	double temp,sumDiff=0;
	for(r=0;r<rows;r++)
		for(c=0;c<cols;c++){
			temp = dest[r][c]- source[r][c];
			sumDiff +=temp*temp;
			
			dest[r][c]=source[r][c];
		
		}
	return Math.sqrt(sumDiff/(rows*cols));
}



public double[][] newMatrixseg(int rows, int cols){
                return(new double[rows+1][cols+1]);
}

}
