package org.babblequest.docex;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.lang.Math;

import java.util.HashMap;

import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import org.apache.mahout.math.DenseMatrix;
import org.apache.mahout.math.Matrix;
import org.apache.mahout.math.SingularValueDecomposition;

/**
 * Plugin class for Lighting adjustment of text image taken from a camera. 
 * This plugin can be used to correct documents with multiple folds or lighting sources.
 * 
 **/

// CHECKSTYLE:OFF
public class L2_Adjust implements PlugInFilter {
  // CHECKSTYLE:ON

  /** Supported imagej image types. */
  private int flags = DOES_16;

  /** Image pixels. */
  protected int[][] pixels;

  /** Image width. */
  int width;

  /** Image height. */
  int height;

  /** Image samples. */
  final int window = 300;

  /*
   * (non-Javadoc)
   * 
   * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
   */
  public int setup(String argv, ImagePlus imp) {
    if (IJ.versionLessThan("1.38x")) {
      return DONE;
    }
    return flags;
  }

  /*
   * (non-Javadoc)
   * 
   * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
   */
  public void run(ImageProcessor ip) {

    // ip = ip.convertToByte(false);
    // ip.autoThreshold();
    pixels = ip.getIntArray();
    width = ip.getWidth();
    height = ip.getHeight();
    filter(pixels, height, width, ip);
  }

//return the integer between 0 and 255 closest to c
  public static int truncate(double c) {
      if (c <= 0.0) return 0;
      if (c >= 255.0) return 254;
      return (int) (Math.round(c));
  }

  // From https://introcs.cs.princeton.edu/java/95linear/
  //    9.5 Numerical Linear Algebra
  public static Matrix KL(Matrix A, int r) {
      int m = A.numRows();
      int n = A.numCols();
      SingularValueDecomposition svd = new SingularValueDecomposition(A);
      
      Matrix Ur = svd.getU().viewPart(0, m, 0, r);  // first r columns of U
      Matrix Vr = svd.getV().viewPart(0, n, 0, r);  // first r columns of V
      Matrix Sr = svd.getS().viewPart(0, r, 0, r);  // first r rows and columns of S
      return Ur.times(Sr).times(Vr.transpose());
  }

  /**
   * Filter interface for imageJ.
   *
   * @param pixels image pixels
   * @param rows
   * @param cols
   * @param ImageProcessor
   */
  public void filter(int[][] pixels, int rows, int cols, ImageProcessor ip) {
    
    // Need to make scaled image for darkest and lightest windows
        int scale = 30;
        
        // create two scaled images one taking min values the other max values
        DenseMatrix AMin = new DenseMatrix((cols/scale)+1, (rows/scale)+1);
        DenseMatrix AMax = new DenseMatrix((cols/scale)+1, (rows/scale)+1);
        int edgeCount[][] = new int[cols/scale+1][rows/scale+1];
        
        System.out.println("cols = " + cols/scale + " rows = " + rows/scale);
        for (int x=0;x<cols/scale;x++)
          for (int y=0;y<rows/scale;y++)
            AMin.set(x,y,Integer.MAX_VALUE);
        
        Sobel_Edges sobel = new Sobel_Edges();
        
        int edges[][] = sobel.filter(pixels, cols, rows);
        
        for (int x=0;x<cols;x++)
           for (int y=0;y<rows;y++)
           {
              int value = pixels[x][y];
              int xscaled = x/scale;
              int yscaled = y/scale;
              
              if (value < AMin.get(xscaled,yscaled))
                AMin.set(xscaled,yscaled,value);
              if (value > AMax.get(xscaled,yscaled))
                AMax.set(xscaled,yscaled,value);
              
              if (edges[x][y] == 255)
                edgeCount[xscaled][yscaled]++;              
           }
        
        System.out.println("Calculating SVD");
        Matrix ArMax = KL(AMax, 3);
        Matrix ArMin = KL(AMin, 3);
        System.out.println("Setting pixels");
               
        // find sobel edges
        
        
        int threshold = ip.getAutoThreshold();
        System.out.println("autothreshold  = " + threshold);
        for (int x=0;x<cols;x++) // was - scale+1
          for (int y=0;y<rows;y++) {
            
              double pix = pixels[x][y];
              double min = ArMin.get(x/scale,y/scale);
              double max = ArMax.get(x/scale,y/scale);
              
              int newPix = (int)Math.round((pix/(max)) * 255.0);
              
              if (edgeCount[x/scale][y/scale] > 20)
              {
                //newPix = (int)Math.round(((pix - min)/(max - min)) * 255.0);  // minMax[0.0-1.0] * max pixel value
                newPix = (int)Math.round(((pix)/(max)) * 255.0);  // minMax[0.0-1.0] * max pixel value

                //int newPix = (int)Math.round(pix-min);
                //int newPix = (int)Math.round(pix-min);
                //if (pix > threshold)
                // newPix = (int)Math.round(pix+max);
              }
              
              ip.putPixelValue(x, y, truncate(newPix));
          }
  }
}
