# docex-ij2-plugins

## Document Exploitation Plugins for ImageJ

# Image processing tasks for document exploitation

* Page Alignment (Radno Transform)
* Lighting Correction
   * Off axis lighting (Bicubic)
   * Folds and multiple light sources (L2)
* Image Segmentation (Split Bregman)

Recommend using the Fiji or imageJ2 of imagej as maven integration is as simple as:
   * mvn -Dimagej.app.directory=/path/to/ImageJ.app/

Plugins can be found under plugins->DOCEX

# Dependencies
  * imagej
  * mahout math (SVD)
  * apache commons math3 (Bicubic)

#Sample Images

![Sample Image](data/sampleImage.jpg "Sample Page")
![Off axis adjustment](data/offaxisAdjust.jpg "Sample Page")
![Multi axis adjustment](data/multiAxisAdjust.jpg "Sample Page")
![Split Bregman segmentation](data/SplitBregman.jpg "Sample Page")
![Deskew](data/deskew.jpg "Sample Page")

All code is freely availabile under the GNU General Public License.

For more information email Steven Parker (parker@babblequest.org)

