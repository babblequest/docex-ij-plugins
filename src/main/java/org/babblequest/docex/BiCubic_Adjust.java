package org.babblequest.docex;

import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.process.*;
import ij.plugin.filter.GaussianBlur;
import ij.measure.Calibration;
import ij.gui.Roi;
import java.util.Vector;
import java.util.*;
import java.awt.*;

import java.awt.image.*;

import org.apache.mahout.math.*;
//import cern.colt.matrix.*;
//import cern.colt.matrix.impl.*;

//import edu.umbc.cs.maple.utils.*;

/***********************************************
   Lighting adjustment for text image taken from a camera.
   For more information on this see http://randomvectors.com/?p=121

   Enhancement just do bight colors and dark colors seperatly. As
the white balance will effect the fonts as well. Perhaps take 
auto threshold calculation as cutoff. See autothreshold plugin
************************************************/

public class BiCubic_Adjust implements PlugInFilter 
{
  private int flags = DOES_ALL;
  protected int[][] pixels;
  private int binaryCount, binaryBackground;
  int width;
  int height;
  private ImageProcessor ip;

  double sample = 300.0;

  private class Point
  {
     int x;
     int y;
     double value;

     public Point(int x, int y, double value)
     {
        this.x = x;
        this.y = y;
        this.value = value;
     }
   
     public int getX()
     {
        return(x);
     }

     public int getY()
     {
        return(y);
     }

     public double getValue()
     {
        return(value);
     } 

     public String toString()
     {
       return(x+","+y+"="+value);
     }
  }


  public int setup(String argv, ImagePlus imp)
  {
      if (IJ.versionLessThan("1.38x"))        // generates an error message for older versions
            return DONE;
        return flags;
  }

  public void run(ImageProcessor ip) {

    //ip = ip.convertToByte(false);
    //ip.autoThreshold();
    this.ip = ip;

    pixels = ip.getIntArray();
    width = ip.getWidth();
    height = ip.getHeight();
    filter(pixels, height, width, ip);
  }
 
   public double pow(double x, int pow)
   {
      double p;
      if (pow == 0)
         p = 1;
      else
      {
        p = x;
        for (int i=1;i<pow;i++)
           p = p * p;
      }

      return(p);
   }
 
   public Matrix calcFit(ArrayList<Point> measures)
   {
     int cols = 16;
     int rows = measures.size();;

     Matrix X = new DenseMatrix(rows,cols);
     Matrix Y = new DenseMatrix(rows,1);

     for (int k=0;k<measures.size();k++)
     {
        Point measure = measures.get(k);
        Y.set(k,0,measure.getValue());

        double x = measure.getX();
        double y = measure.getY();
      
        int col=0;
        for (int i=0;i<4;i++)
          for (int j=0;j<4;j++)
          {
             double cal = pow(x,i)*pow(y,j);
             X.set(k,col++,cal);
          }
     } 

     // solve SVD
      SingularValueDecomposition svd = new SingularValueDecomposition(X);
      Matrix U = svd.getU();
      Matrix V = svd.getV();
      Matrix S = svd.getS();

      // S is the diagonal matrix 
      for (int i=0;i<16;i++)
      {
        double val = 1.0/S.get(i,i);
        //if (val < 1.0E-5)
            //val = 0.0;
        if (S.get(i,i) == 0)
          val = 0.0;
          
        S.set(i,i, val); 
      }

      // X phsuedoinverse
      Matrix VS = V.times(S);
      System.out.println(VS);
      Matrix Xinv = VS.times(U.transpose());

      // matrix has coeffients
      Matrix C = Xinv.times(Y);

      Matrix nY = X.times(C);

      Matrix dY = Y.minus(nY);
      //double sum = ColtUtils.dotproduct(ColtUtils.getcol(dY,0),ColtUtils.getcol(dY,0));
      double sum = dY.viewColumn(0).dot(dY.viewColumn(0));

      return(C);
   }
 
   public double calc(Matrix C, double x, double y)
   {
      double value=0;
      int col = 0;

      for (int i=0;i<4;i++)
        for (int j=0;j<4;j++)
        {
           value = value + (C.get(col++,0) * (pow(x,i)*pow(y,j)));
        }

      return(value);
   }

   public void filter(int [][] pixels, int rows, int cols, ImageProcessor ip) 
   {
     int offset, sum1, sum2=0, sum=0;
     int rowOffset = width;
     int count;

     ArrayList<Point> maxMeasures = new ArrayList<Point>();
     ArrayList<Point> minMeasures = new ArrayList<Point>();
     ArrayList<Point> meanMeasures = new ArrayList<Point>();

     int ys = 0;
     int xs = 0;

     double max = 0;
     double min = 99999;
     double mean=0.0;

     for (int y=0; y<rows-(int)sample; y=y+(int)sample) 
     {
        xs = 0;
        for (int x=0; x<cols-(int)sample; x=x+(int)sample) 
        {
            max = 0;
            min = 99999;
            mean=0.0;

            count = 0;

            for(int yi=y; yi<=y+sample; yi++)  {
		   for(int xi=x; xi<=x+sample; xi++)  {

                      int pix = pixels[xi][yi];
 
                      if (max < pix)
                        max = pix;
                        
                      if (min > pix)
                        min = pix;

                      mean = mean+(double)pix;
                      count++;
		   }
	    }

            mean = mean/(double)count;

            // put max and min at center points in sample square
            maxMeasures.add(new Point(xs,ys,max));
            //maxMeasures.add(new Point(x+(sample),y+(sample),(double)max));
            minMeasures.add(new Point(xs,ys,min));
            meanMeasures.add(new Point(xs++,ys,mean));
            //minMeasures.add(new Point(x+(sample),y+(sample),(double)min));

/***
for(int yi=y; yi<=y+sample; yi++)  {
                   for(int xi=x; xi<=x+sample; xi++)  {
        ip.putPixel(xi,yi,max);
}}
***/
        }
        ys++;
     }

/**
     maxMeasures.add(new Point(xs,ys,(double)max));
     minMeasures.add(new Point(xs,ys,(double)min));
     meanMeasures.add(new Point(xs,ys,(double)mean));
**/


     Matrix Cmax = calcFit(maxMeasures);
     Matrix Cmin = calcFit(minMeasures);
     Matrix Cmean = calcFit(meanMeasures);

System.out.println(maxMeasures);
 
     for (int y=0;y<rows;y++)
        for (int x=0;x<cols;x++)
        {
           int pix = ip.getPixel(x,y);

           min = calc(Cmin,((double)x)/sample,((double)y)/sample);
           max = calc(Cmax,((double)x)/sample,((double)y)/sample);
           mean = calc(Cmean,((double)x)/sample,((double)y)/sample);

           // Linear strech
           //OUTVAL = (INVAL - INLO) * ((OUTUP-OUTLO)/(INUP-INLO)) + OUTLO

           double newPix = pix;

           if (pix >85)
             newPix = (pix - min) * ((255.0-0.0)/(max-min)) + 0.0; 

           if (newPix < 0)
              newPix = 0;
           else if (newPix > 255)
              newPix = 255;

/**
           if ((pix < mean) && (pix-mean >= 0))
             newPix = pix-min;

           if ((pix > mean) && (pix+max <= 255))
             newPix = pix+max;
****/
           //double scale = (max-min)/255.0;

               //ip.putPixelValue(x,y,calc(Cmax,x/150,y/150));
               ip.putPixelValue(x,y,newPix);
        }
   }
}
