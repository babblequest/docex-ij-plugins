package org.babblequest.docex;

/*******************************************************************************
 * Steven Parker 2018. All rights reserved.
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
import org.apache.mahout.math.DenseMatrix;
import org.apache.mahout.math.Matrix;
import org.apache.mahout.math.SingularValueDecomposition;

/***
 * Plugin class for Lighting adjustment of text image taken from a camera. 
 * Correction of flash or mounted lighting that is not perpendicular (off axis) to a flat document.
 * 
 * Enhancement just do bight colors and dark colors separately. As the white balance will effect the
 * fonts as well. Perhaps take auto threshold calculation as cutoff. See autothreshold plugin
 ************************************************/

// CHECKSTYLE:OFF
public class BiCubic_Adjust implements PlugInFilter {
  // CHECKSTYLE:ON


  /** Supported imagej image types. */
  private int flags = DOES_ALL;

  /** Image pixels. */
  protected int[][] pixels;

  /** Image width. */
  int width;

  /** Image height. */
  int height;

  /** Image samples. */
  double sample = 300.0;

  /**
   * The Class Point.
   */
  private class Point {

    int xp;
    int yp;
    double value;

    private Point(int x, int y, double value) {
      this.xp = x;
      this.yp = y;
      this.value = value;
    }

    /**
     * Gets the x.
     *
     * @return the x
     */
    public int getX() {
      return (xp);
    }

    /**
     * Gets the y.
     *
     * @return the y
     */
    public int getY() {
      return (yp);
    }

    /**
     * Gets the value.
     *
     * @return the value
     */
    public double getValue() {
      return (value);
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#toString()
     */
    public String toString() {
      return (xp + "," + yp + "=" + value);
    }
  }


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
   * Pow.
   *
   * @param x the x
   * @param pow the pow
   * @return the double
   */
  private double pow(double x, int pow) {
    double p;
    if (pow == 0) {
      p = 1;
    } else {
      p = x;
      for (int i = 1; i < pow; i++) {
        p = p * p;
      }
    }

    return (p);
  }

  /**
   * Calc fit.
   *
   * @param measures the measures
   * @return the matrix
   */
  private Matrix calcFit(ArrayList<Point> measures) {
    int cols = 16;
    int rows = measures.size();;

    Matrix xMatrix = new DenseMatrix(rows, cols);
    Matrix yMatrix = new DenseMatrix(rows, 1);

    for (int k = 0; k < measures.size(); k++) {
      Point measure = measures.get(k);
      yMatrix.set(k, 0, measure.getValue());

      double x = measure.getX();
      double y = measure.getY();

      int col = 0;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          double cal = pow(x, i) * pow(y, j);
          xMatrix.set(k, col++, cal);
        }
      }
    }

    // solve SVD
    SingularValueDecomposition svd = new SingularValueDecomposition(xMatrix);
    Matrix uMatrix = svd.getU();
    Matrix vMatrix = svd.getV();
    Matrix sMatrix = svd.getS();

    // S is the diagonal matrix
    for (int i = 0; i < 16; i++) {
      double val = 1.0 / sMatrix.get(i, i);
      // if (val < 1.0E-5)
      // val = 0.0;
      if (sMatrix.get(i, i) == 0) {
        val = 0.0;
      }

      sMatrix.set(i, i, val);
    }

    // X phsuedoinverse
    Matrix vsMatrix = vMatrix.times(sMatrix);
    System.out.println(vsMatrix);
    Matrix xInverseMatrix = vsMatrix.times(uMatrix.transpose());

    // matrix has coeffients
    Matrix coefficientMatrix = xInverseMatrix.times(yMatrix);

    Matrix nYMatrix = xMatrix.times(coefficientMatrix);

    Matrix dYMatrix = yMatrix.minus(nYMatrix);
    // double sum = ColtUtils.dotproduct(ColtUtils.getcol(dY,0),ColtUtils.getcol(dY,0));
    double sum = dYMatrix.viewColumn(0).dot(dYMatrix.viewColumn(0));

    return (coefficientMatrix);
  }

  /**
   * Calc.
   *
   * @param cMatrix the c
   * @param xPixel 
   * @param yPixel
   * @return bicubic value from matrix points
   */
  private double calc(Matrix cMatrix, double xPixel, double yPixel) {
    double value = 0;
    int col = 0;

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        value = value + (cMatrix.get(col++, 0) * (pow(xPixel, i) * pow(yPixel, j)));
      }
    }

    return (value);
  }

  /**
   * Filter.
   *
   * @param pixels image pixels
   * @param rows
   * @param cols
   * @param ImageProcessor
   */
  public void filter(int[][] pixels, int rows, int cols, ImageProcessor ip) {
    int offset;
    int sum1;
    int sum2 = 0;
    int sum = 0;
    int rowOffset = width;
    int count;

    ArrayList<Point> maxMeasures = new ArrayList<Point>();
    ArrayList<Point> minMeasures = new ArrayList<Point>();
    ArrayList<Point> meanMeasures = new ArrayList<Point>();

    int ys = 0;
    int xs = 0;

    double max = 0;
    double min = 99999;
    double mean = 0.0;

    for (int y = 0; y < rows - (int) sample; y = y + (int) sample) {
      xs = 0;
      for (int x = 0; x < cols - (int) sample; x = x + (int) sample) {
        max = 0;
        min = 99999;
        mean = 0.0;

        count = 0;

        for (int yi = y; yi <= y + sample; yi++) {
          for (int xi = x; xi <= x + sample; xi++) {

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
        maxMeasures.add(new Point(xs, ys, max));
        // maxMeasures.add(new Point(x+(sample),y+(sample),(double)max));
        minMeasures.add(new Point(xs, ys, min));
        meanMeasures.add(new Point(xs++, ys, mean));
        // minMeasures.add(new Point(x+(sample),y+(sample),(double)min));

        /***
         * for(int yi=y; yi<=y+sample; yi++) { for(int xi=x; xi<=x+sample; xi++) {
         * ip.putPixel(xi,yi,max); }}
         ***/
      }
      ys++;
    }

    /**
     * maxMeasures.add(new Point(xs,ys,(double)max)); minMeasures.add(new Point(xs,ys,(double)min));
     * meanMeasures.add(new Point(xs,ys,(double)mean));
     **/


    Matrix cMax = calcFit(maxMeasures);
    Matrix cMin = calcFit(minMeasures);
    Matrix cMean = calcFit(meanMeasures);

    System.out.println(maxMeasures);

    for (int y = 0; y < rows; y++) {
      for (int x = 0; x < cols; x++) {

        min = calc(cMin, ((double) x) / sample, ((double) y) / sample);
        max = calc(cMax, ((double) x) / sample, ((double) y) / sample);
        mean = calc(cMean, ((double) x) / sample, ((double) y) / sample);

        // Linear strech
        // OUTVAL = (INVAL - INLO) * ((OUTUP-OUTLO)/(INUP-INLO)) + OUTLO

        int pix = ip.getPixel(x, y);
        double newPix = pix;

        if (pix > 85) {
          newPix = (pix - min) * ((255.0 - 0.0) / (max - min)) + 0.0;
        }

        if (newPix < 0) {
          newPix = 0;
        } else if (newPix > 255) {
          newPix = 255;
        }

        /***
         * if ((pix < mean) && (pix-mean >= 0)) newPix = pix-min;
         * 
         * if ((pix > mean) && (pix+max <= 255)) newPix = pix+max;
         ****/
        // double scale = (max-min)/255.0;

        // ip.putPixelValue(x,y,calc(Cmax,x/150,y/150));
        ip.putPixelValue(x, y, newPix);
      }
    }
  }
}
