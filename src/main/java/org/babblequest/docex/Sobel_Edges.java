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
import java.awt.*;

import java.awt.image.*;

/******************************************************
  from http://www.pages.drexel.edu/~weg22/edge.html
******************************************************/

public class Sobel_Edges implements PlugInFilter 
{
  private int flags = DOES_ALL|DOES_8G|SNAPSHOT;
  protected int[][] pixels;
  private int binaryCount, binaryBackground;
  int width;
  int height;
  private ImageProcessor ip;


   /* 3x3 GX Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
   int GX[][] = { {-1, 0, 1}, {-2,0,2},{-1,0,1}};
   //GX[0][0] = -1; GX[0][1] = 0; GX[0][2] = 1;
   //GX[1][0] = -2; GX[1][1] = 0; GX[1][2] = 2;
   //GX[2][0] = -1; GX[2][1] = 0; GX[2][2] = 1;

   /* 3x3 GY Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
   int GY[][] = { {1,2,1},{0,0,0},{-1,-2,-1}};
   //GY[0][0] =  1; GY[0][1] =  2; GY[0][2] =  1;
   //GY[1][0] =  0; GY[1][1] =  0; GY[1][2] =  0;
   //GY[2][0] = -1; GY[2][1] = -2; GY[2][2] = -1;


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
    filter();
  }

   public void filter() 
   {
     int p1, p2, p3, p4, p5, p6, p7, p8, p9;
     int offset, sum1, sum2=0, sum=0;
     int rowOffset = width;
     int count;

     int xMin=1;
     int xMax=width-2; 
     int yMin=1; 
     int yMax=height-2;

     for (int y=yMin; y<=yMax; y++) 
     {
        for (int x=xMin; x<=xMax; x++) 
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
               //sum = Math.abs(sumX) + Math.abs(sumY);

               if (Math.abs(sumX) < Math.abs(sumY))
                  sum = Math.abs(sumY);
               else
                  sum = Math.abs(sumX);

            if(sum>255) sum=255;
             if(sum<0) sum=0;
           
           ip.set(x,y,sum&0xff);
        }
     }
   }
}
