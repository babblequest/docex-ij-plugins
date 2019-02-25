package org.babblequest.docex;

/*******************************************************************************
 *
 * Licensed under the GNU General Public License, Version 2 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   https://www.gnu.org/licenses/gpl-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */ 

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.lang.Math;
import org.apache.commons.math3.analysis.BivariateFunction;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;

import org.apache.commons.math3.analysis.interpolation.BivariateGridInterpolator;

import org.apache.mahout.math.DenseMatrix;
import org.apache.mahout.math.Matrix;
import org.apache.mahout.math.SingularValueDecomposition;

/**
 * Plugin class for Lighting adjustment of text image taken from a camera. 
 * Correction of flash or mounted lighting that is not perpendicular (off axis) to a flat document.
 * 
 * This will not correct photographs of pocket trash with multiple folds.
 * 
 * Discribed in:
 *      Image Processing Handbook by John Russ 2002 pg. 47
 *      http://www.geo.uzh.ch/microsite/rsl-documents/research/SARlab/GMTILiterature/PDF/1142_CH03.pdf
 **/    
 

// CHECKSTYLE:OFF
public class BiCubic_Adjust implements PlugInFilter {
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

  /**
   * Calculate least squares fit for the bicubic polynomial corresponding to lighting that is off axis.
   *
   * @param measures the measures
   * @return the matrix
   */
  private BivariateFunction calcFit(int xSamples[], int ySamples[], double measures[][]) {
    
    BivariateGridInterpolator interpolator = new BicubicInterpolator();
    
    double xValues[] = new double[xSamples.length];
    for (int i=0;i<xSamples.length;i++)
      xValues[i] = (double)xSamples[i];
    
    double yValues[] = new double[ySamples.length];
    for (int i=0;i<ySamples.length;i++)
      yValues[i] = (double)ySamples[i];
    
    BivariateFunction function = interpolator.interpolate(xValues, yValues,measures);
    
    return function;
  }

  /**
   * Create bicubic factors to fit via SVD
   *
   * @param coefficientMatrix the c
   * @param xPixel 
   * @param yPixel
   * @return bicubic value from matrix points
   */
  private double calc(BivariateFunction function, double xPixel, double yPixel) {
    return(function.value(xPixel,  yPixel));
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
    int count;

    double max = 0;
    double min = Double.MAX_VALUE;
    double mean = 0.0;

    // FIXME need to handle edge case where window is exact divisor of windows. Now
    // adding max values to avoid out of bounds issues with interpolation functions
    
    // Take samples over window size
    int[] xSamples = new int[(cols/window)+1];
    int[] ySamples = new int[(rows/window)+1];
    
    //define sorted list of grid points
    int yindex = 0;
    for (int y = 0; y < rows - window; y = y + window) {
      ySamples[yindex++] = y;
    }
    ySamples[yindex] = rows-1;
    
    for (int i=0;i<ySamples.length;i++) {
      System.out.println(ySamples[i]);
    }
    
    int xindex = 0;
    for (int x = 0; x < cols - window; x = x + window) {
      xSamples[xindex++] = x;
    }
    xSamples[xindex] = cols-1;
    
    double maxMeasures[][] = new double[xSamples.length][ySamples.length];
    double minMeasures[][] = new double[xSamples.length][ySamples.length];
    double meanMeasures[][] = new double[xSamples.length][ySamples.length];
    
    for (yindex = 0; yindex < ySamples.length; yindex++) {
      int y = ySamples[yindex];
      for (xindex = 0; xindex < xSamples.length; xindex++) {
        int x = xSamples[xindex];
        max = 0.0;
        min = Double.MAX_VALUE;
        mean = 0.0;

        count = 0;

        int yi = y;
        int ystep = window;
        if (yindex == ySamples.length-1) { // handle edge case at max of image
          yi = y-window;
          ystep = 0;
        }  
        for (; yi <= y + ystep; yi++) {
          
          int xi = x;
          int xstep = window;
          if (xindex == xSamples.length-1) {
            xi = x-window;
            xstep = 0;
          }
          
          for (; xi <= x + xstep; xi++) {

            int pix = pixels[xi][yi];
            if (max < pix) {
              max = pix;
            }

            if (min > pix) {
              min = pix;
            }

            mean = mean + (double) pix;
            count++;
          }
        }

        mean = mean / (double) count;

        // put max and min at center points in sample square
        maxMeasures[xindex][yindex] =  max;
        
        // maxMeasures.add(new Point(x+(sample),y+(sample),(double)max));
        minMeasures[xindex][yindex] = min;
        
        meanMeasures[xindex][yindex] = mean;
        // minMeasures.add(new Point(x+(sample),y+(sample),(double)min));
        
      }
    }

    BivariateFunction cMax = calcFit(xSamples, ySamples, maxMeasures);
    BivariateFunction cMin = calcFit(xSamples, ySamples, minMeasures);
    BivariateFunction cMean = calcFit(xSamples, ySamples, meanMeasures);

    double newPixels[][] = new double[pixels.length][pixels[0].length];
    for (int y = 0; y < rows; y++) {
      for (int x = 0; x < cols; x++) {

        // Calculate point values keeping in mind samples correspond to windows not actual point values
        min = calc(cMin, (double) x, (double) y);
        max = calc(cMax, (double) x, (double) y);;
        mean = calc(cMean, (double) x, (double) y);

        // Linear strech
        // OUTVAL = (INVAL - INLO) * ((OUTUP-OUTLO)/(INUP-INLO)) + OUTLO

        int pix = pixels[x][y];
        double newPix = (double)pix;

       /* if (pix > 85) {
          newPix = (pix - min) * ((255.0 - 0.0) / (max - min)) + 0.0;
        }*/

        /*if (newPix < 0) {
          newPix = 0;
        } else if (newPix > 255) {
          newPix = 255;
        }*/

        // min max scaling
        newPix = ((pix - min)/(max - min)) * 255.0;  // minMax[0.0-1.0] * max pixel value
        
        //System.out.println("min " + min + " max " + max + " pix "+ pix + " min-max " + newPix);

        // simple lighting correction
        //newPix = pix+min;
        
        /***
         * if ((pix < mean) && (pix-mean >= 0)) newPix = pix-min;
         * 
         * if ((pix > mean) && (pix+max <= 255)) newPix = pix+max;
         ****/
        // double scale = (max-min)/255.0;

        ip.putPixelValue(x, y, newPix);
      }
    }
  }
}
