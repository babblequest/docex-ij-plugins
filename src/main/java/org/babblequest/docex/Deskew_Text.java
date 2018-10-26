package org.babblequest.docex;

/*******************************************************************************
 * Derived from Sergey Zolotaryov's deskew code 2010
 *    Original code is at https://anydoby.com/jblog/en/java/1990.
 * 
 * Modified by Steven Lee on 2013
 *    https://gist.github.com/witwall/5565179
 * 
 * Integration into ImageJ Steven Parker 2015
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
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;

/**********************************************
 * Text image deskew using code derived from an abandoned Gimp autodeskew plugin - i it uses GPL
 * version of the Radon transform. You may find the original auther here for the initial java port.
 * https://gist.github.com/witwall/5565179
 **********************************************/

public class Deskew_Text implements PlugInFilter {
  private int flags = DOES_ALL | CONVERT_TO_FLOAT;

  public int setup(String argv, ImagePlus imp) {
    if (IJ.versionLessThan("1.38x"))
    {
      return DONE;
    }
    return flags;
  }

  public void run(ImageProcessor ip) {
    Image image = ip.createImage();
    final double skewRadians;
    BufferedImage black = new BufferedImage(image.getWidth(null), image.getHeight(null),
        BufferedImage.TYPE_BYTE_BINARY);
    final Graphics2D g = black.createGraphics();
    g.drawImage(image, 0, 0, null);
    g.dispose();

    skewRadians = findSkew(black);
    ip.rotate(-1 * -57.295779513082320876798154814105 * skewRadians);
    // System.out.println(-57.295779513082320876798154814105 * skewRadians);
  }

  static int getByteWidth(final int width) {
    return (width + 7) / 8;
  }

  static int next_pow2(final int n) {
    int retval = 1;
    while (retval < n) {
      retval <<= 1;
    }
    return retval;
  }

  static class BitUtils {
    static int[] bitcount = new int[256];
    static int[] invbits_ = new int[256];

    static {
      for (int i = 0; i < 256; i++) {
        int j = i;
        int cnt = 0;
        do {
          cnt += j & 1;
        } while ((j >>= 1) != 0);
        int x = (i << 4) | (i >> 4);
        x = ((x & 0xCC) >> 2) | ((x & 0x33) << 2);
        x = ((x & 0xAA) >> 1) | ((x & 0x55) << 1);
        bitcount[i] = cnt;
        invbits_[i] = x;
      }
    }
  }

  static double findSkew(final BufferedImage img) {
    final DataBuffer buffer = img.getRaster().getDataBuffer();
    final int byteWidth = getByteWidth(img.getWidth());
    final int padmask = 0xFF << ((img.getWidth() + 7) % 8);
    int elementIndex = 0;
    for (int row = 0; row < img.getHeight(); row++) {
      for (int col = 0; col < byteWidth; col++) {
        int elem = buffer.getElem(elementIndex);
        elem ^= 0xff;// invert colors
        elem = BitUtils.invbits_[elem]; // Change the bit order
        buffer.setElem(elementIndex, elem);
        elementIndex++;
      }
      final int lastElement = buffer.getElem(elementIndex - 1) & padmask;
      buffer.setElem(elementIndex - 1, lastElement); // Zero trailing bits
    }
    final int w2 = next_pow2(byteWidth);
    final int ssize = 2 * w2 - 1; // Size of sharpness table
    final int[] sharpness = new int[ssize];
    radon(img.getWidth(), img.getHeight(), buffer, 1, sharpness);
    radon(img.getWidth(), img.getHeight(), buffer, -1, sharpness);
    int i;
    int imax = 0;
    int vmax = 0;
    double sum = 0.;
    for (i = 0; i < ssize; i++) {
      final int s = sharpness[i];
      if (s > vmax) {
        imax = i;
        vmax = s;
      }
      sum += s;
    }
    final int h = img.getHeight();
    if (vmax <= 3 * sum / h) { // Heuristics !!!
      return 0;
    }
    final double iskew = imax - w2 + 1;
    return Math.atan(iskew / (8 * w2));
  }

  static void radon(final int width, final int height, final DataBuffer buffer, final int sign,
      final int[] sharpness) {

    int[] p1;
    int[] p2; // Stored columnwise

    final int w2 = next_pow2(getByteWidth(width));
    final int w = getByteWidth(width);
    final int h = height;

    final int s = h * w2;
    p1 = new int[s];
    p2 = new int[s];
    // Fill in the first table
    int row;
    int column;
    int scanlinePosition = 0;
    for (row = 0; row < h; row++) {
      scanlinePosition = row * w;
      for (column = 0; column < w; column++) {
        if (sign > 0) {
          int b = buffer.getElem(0, scanlinePosition + w - 1 - column);
          p1[h * column + row] = BitUtils.bitcount[b];
        } else {
          int b = buffer.getElem(0, scanlinePosition + column);
          p1[h * column + row] = BitUtils.bitcount[b];
        }
      }
    }

    int[] x1 = p1;
    int[] x2 = p2;
    // Iterate
    int step = 1;
    for (;;) {
      int i;
      for (i = 0; i < w2; i += 2 * step) {
        int j;
        for (j = 0; j < step; j++) {
          // Columns-sources:
          final int s1 = h * (i + j);// x1 pointer
          final int s2 = h * (i + j + step); // x1 pointer

          // Columns-targets:
          final int t1 = h * (i + 2 * j); // x2 pointer
          final int t2 = h * (i + 2 * j + 1); // x2 pointer
          int m;
          for (m = 0; m < h; m++) {
            x2[t1 + m] = x1[s1 + m];
            x2[t2 + m] = x1[s1 + m];
            if (m + j < h) {
              x2[t1 + m] += x1[s2 + m + j];
            }
            if (m + j + 1 < h) {
              x2[t2 + m] += x1[s2 + m + j + 1];
            }
          }
        }
      }

      // Swap the tables:
      final int[] aux = x1;
      x1 = x2;
      x2 = aux;
      // Increase the step:
      step *= 2;
      if (step >= w2) {
        break;
      }
    }
    // Now, compute the sum of squared finite differences:
    for (column = 0; column < w2; column++) {
      int acc = 0;
      final int col = h * column;
      for (row = 0; row + 1 < h; row++) {
        final int diff = x1[col + row] - x1[col + row + 1];
        acc += diff * diff;
      }
      sharpness[w2 - 1 + sign * column] = acc;
    }
  }
}
