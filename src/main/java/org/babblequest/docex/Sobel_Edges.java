package org.babblequest.docex;

/*******************************************************************************
 * Implementation of Sobel Edges there are many examples.
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

public class Sobel_Edges implements PlugInFilter {
  private int flags = DOES_ALL | DOES_8G | SNAPSHOT;
  int width;
  int height;

  /* 3x3 GX Sobel mask. Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
  int GX[][] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};


  /* 3x3 GY Sobel mask. Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
  int GY[][] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};


  public int setup(String argv, ImagePlus imp) {
    if (IJ.versionLessThan("1.38x")) // generates an error message for older versions
      return DONE;
    return flags;
  }

  public void run(ImageProcessor ip) {
    //ip = ip.convertToByte(false);
    int[][] pixels = ip.getIntArray();
    int[][] newPixels = filter(pixels, ip.getWidth(), ip.getHeight());
    ip.setIntArray(newPixels);
  }

  public int[][] filter(int[][] pixels, int width, int height) {    
    int newPixels[][] = new int[width][height];
    int xMin = 1;
    int xMax = width - 2;
    int yMin = 1;
    int yMax = height - 2;

    for (int y = yMin; y <= yMax; y++) {
      for (int x = xMin; x <= xMax; x++) {
        int sumX = 0;
        int sumY = 0;

        for (int I = -1; I <= 1; I++) {
          for (int J = -1; J <= 1; J++) {
            sumX = sumX + (pixels[x + I][y + J] * GX[I + 1][J + 1]);
          }
        }

        /*-------Y GRADIENT APPROXIMATION-------*/
        for (int I = -1; I <= 1; I++) {
          for (int J = -1; J <= 1; J++) {
            sumY = sumY + (pixels[x + I][y + J] * GY[I + 1][J + 1]);
          }
        }

        /*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
        // sum = Math.abs(sumX) + Math.abs(sumY);
        int sum = Math.abs(sumX);
        if (Math.abs(sumX) < Math.abs(sumY))
          sum = Math.abs(sumY);

        if (sum > 255)
          sum = 255;
        if (sum < 0)
          sum = 0;

        newPixels[x][y] = sum & 0xff;
      }
    }
    return(newPixels);
  }
}
