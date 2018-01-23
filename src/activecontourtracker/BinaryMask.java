/*
 * A copy of ImageJ's Huang Dark auto-thresholder
 */
package activecontourtracker;

import java.util.Arrays;

/**
 *
 * @author edwardharry
 */
public
        class BinaryMask
{

    static private
            int[] MakeHistogram(byte[] image)
    {
        int[] histogram = new int[256];
        for (int index = 0;
             index < image.length;
             ++index)
        {
            ++histogram[image[index] & 0xff];
        }

        return histogram;
    }

    static
            byte[] MakeMask(short[] image)
    {
        // Make a scaled byte image first
        short min = Short.MAX_VALUE;
        short max = Short.MIN_VALUE;
        for (int index = 0;
             index < image.length;
             ++index)
        {
            short value = image[index];
            if (value > max)
            {
                max = value;
            }
            if (value < min)
            {
                min = value;
            }
        }
        double scale = 256.0 / (max - min + 1);
        byte[] scaledByteImage = new byte[image.length];
        for (int index = 0;
             index < image.length;
             ++index)
        {
            int value = (image[index] & 0xffff) - min;
            if (value < 0)
            {
                value = 0;
            }
            value = (int) ((value * scale) + 0.5);
            if (value > 255)
            {
                value = 255;
            }
            scaledByteImage[index] = (byte) value;
        }

        // Apply Huang's algorithm
        int threshold = Huang(MakeHistogram(scaledByteImage));

        // Set and scale the thresholds
        double lower = threshold + 1;
        double upper = max;
        if (lower > 255)
        {
            lower = 255;
        }
        if (max > min)
        {
            lower = min + (lower / 255.0) * (max - min);
        }
        else
        {
            lower = upper = min;
        }
        if (lower < 0.0)
        {
            lower = 0.0;
        }
        if (upper > 65535.0)
        {
            upper = 65535.0;
        }

        // Make the mask
        byte[] mask = new byte[image.length];
        for (int index = 0;
             index < image.length;
             ++index)
        {
            double value = image[index] & 0xffff;
            if (value >= lower && value <= upper)
            {
                mask[index] = (byte) 255;
            }
        }

        return mask;
    }

    // Copied from ImageJ
    static
            int Huang(int[] data)
    {
        // Implements Huang's fuzzy thresholding method
        // Uses Shannon's entropy function (one can also use Yager's entropy function)
        // Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing
        // the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
        // M. Emre Celebi  06.15.2007
        // Ported to ImageJ plugin by G. Landini from E Celebi's fourier_0.8 routines
        int threshold;
        int ih, it;
        int first_bin;
        int last_bin;
        double sum_pix;
        double num_pix;
        double term;
        double ent;  // entropy
        double min_ent; // min entropy
        double mu_x;

        /* Determine the first non-zero bin */
        first_bin = 0;
        for (ih = 0; ih < 256; ih++)
        {
            if (data[ih] != 0)
            {
                first_bin = ih;
                break;
            }
        }

        /* Determine the last non-zero bin */
        last_bin = 255;
        for (ih = 255; ih >= first_bin; ih--)
        {
            if (data[ih] != 0)
            {
                last_bin = ih;
                break;
            }
        }
        term = 1.0 / (double) (last_bin - first_bin);
        double[] mu_0 = new double[256];
        sum_pix = num_pix = 0;
        for (ih = first_bin; ih < 256; ih++)
        {
            sum_pix += (double) ih * data[ih];
            num_pix += data[ih];
            /* NUM_PIX cannot be zero ! */
            mu_0[ih] = sum_pix / num_pix;
        }

        double[] mu_1 = new double[256];
        sum_pix = num_pix = 0;
        for (ih = last_bin; ih > 0; ih--)
        {
            sum_pix += (double) ih * data[ih];
            num_pix += data[ih];
            /* NUM_PIX cannot be zero ! */
            mu_1[ih - 1] = sum_pix / (double) num_pix;
        }

        /* Determine the threshold that minimizes the fuzzy entropy */
        threshold = -1;
        min_ent = Double.MAX_VALUE;
        for (it = 0; it < 256; it++)
        {
            ent = 0.0;
            for (ih = 0; ih <= it; ih++)
            {
                /* Equation (4) in Ref. 1 */
                mu_x = 1.0 / (1.0 + term * Math.abs(ih - mu_0[it]));
                if (!((mu_x < 1e-06) || (mu_x > 0.999999)))
                {
                    /* Equation (6) & (8) in Ref. 1 */
                    ent += data[ih] * (-mu_x * Math.log(mu_x) - (1.0 - mu_x) * Math.log(1.0 - mu_x));
                }
            }

            for (ih = it + 1; ih < 256; ih++)
            {
                /* Equation (4) in Ref. 1 */
                mu_x = 1.0 / (1.0 + term * Math.abs(ih - mu_1[it]));
                if (!((mu_x < 1e-06) || (mu_x > 0.999999)))
                {
                    /* Equation (6) & (8) in Ref. 1 */
                    ent += data[ih] * (-mu_x * Math.log(mu_x) - (1.0 - mu_x) * Math.log(1.0 - mu_x));
                }
            }
            /* No need to divide by NUM_ROWS * NUM_COLS * LOG(2) ! */
            if (ent < min_ent)
            {
                min_ent = ent;
                threshold = it;
            }
        }
        return threshold;
    }

    static
            void SobelFilter(byte[] mask, Dims dims)
    {
        KernelFilter(mask, dims, null, 0, FilterType.Sobel);
    }

    static
            void MedianFilter(byte[] mask, Dims dims, FilterKernel kernel)
    {
        int halfArea = ((kernel.kNPoints + 1) & ~1) / 2;
        KernelFilter(mask, dims, kernel, halfArea, FilterType.Sum);
    }

    static
            void ErodeFilter(byte[] mask, Dims dims, FilterKernel kernel)
    {
        KernelFilter(mask, dims, kernel, kernel.kNPoints, FilterType.Sum);
    }

    static
            class FilterKernel
    {

        int[] lineRadii;
        int kHeight;
        int kRadius;
        int kNPoints;

        //Copied from ImageJ
        /**
         * Create a circular kernel (structuring element) of a given radius.
         *
         * @param radius Radius = 0.5 includes the 4 neighbors of the pixel in
         * the center, radius = 1 corresponds to a 3x3 kernel size.
         * @return the circular kernel The output is an array that gives the
         * length of each line of the structuring element (kernel) to the left
         * (negative) and to the right (positive): [0] left in line 0, [1] right
         * in line 0, [2] left in line 2, ... The maximum (absolute) value
         * should be kernelRadius. Array elements at the end: length-2: nPoints,
         * number of pixels in the kernel area length-1: kernelRadius in x
         * direction (kernel width is 2*kernelRadius+1) Kernel height can be
         * calculated as (array length - 1)/2 (odd number); Kernel radius in y
         * direction is kernel height/2 (truncating integer division). Note that
         * kernel width and height are the same for the circular kernels used
         * here, but treated separately for the case of future extensions with
         * non-circular kernels.
         */
        FilterKernel(double radius)
        {
            if (radius >= 1.5 && radius < 1.75) //this code creates the same sizes as the previous RankFilters
            {
                radius = 1.75;
            }
            else if (radius >= 2.5 && radius < 2.85)
            {
                radius = 2.85;
            }
            int r2 = (int) (radius * radius) + 1;

            kRadius = (int) (Math.sqrt(r2 + 1e-10));
            kHeight = 2 * kRadius + 1;

            lineRadii = new int[2 * kHeight];
            lineRadii[2 * kRadius] = -kRadius;
            lineRadii[2 * kRadius + 1] = kRadius;

            kNPoints = 2 * kRadius + 1;

            for (int y = 1; y <= kRadius; y++)
            {                //lines above and below center together
                int dx = (int) (Math.sqrt(r2 - y * y + 1e-10));
                lineRadii[2 * (kRadius - y)] = -dx;
                lineRadii[2 * (kRadius - y) + 1] = dx;
                lineRadii[2 * (kRadius + y)] = -dx;
                lineRadii[2 * (kRadius + y) + 1] = dx;
                kNPoints += 4 * dx + 2;      //2*dx+1 for each line, above&below
            }
        }
    }

    enum FilterType
    {
        Sum, Sobel;
    }

    // A lot of this is copied from ImageJ
    private static
            void KernelFilter(byte[] mask, Dims dims, FilterKernel kernel,
                              int ForegroundMinCount, FilterType filterType)
    {
        if (filterType == FilterType.Sobel)
        {
            kernel = new FilterKernel(1.0);
        }

        int[] sobelKernel =
        {
            -1, 1, 0, 1, 1, 1, -1, 0, 0, 0, 1,
            0, -1, -1, 0, -1, 1, -1
        };

        //int halfArea = ((kernel.kNPoints + 1) & ~1) / 2;
        int cacheWidth = dims.width + (2 * kernel.kRadius);
        byte[] cache = new byte[cacheWidth * kernel.kHeight];

        int[] cachePointers = new int[2 * kernel.kHeight];
        for (int i = 0;
             i < kernel.kHeight;
             ++i)
        {
            cachePointers[2 * i] = i * cacheWidth + kernel.kRadius
                                   + kernel.lineRadii[2 * i];
            cachePointers[2 * i + 1] = i * cacheWidth + kernel.kRadius
                                       + kernel.lineRadii[2 * i + 1];
        }

        int previousY = kernel.kHeight / 2 - kernel.kHeight;

        for (int y = 0;
             y < dims.height;
             ++y)
        {
            for (int i = 0;
                 i < cachePointers.length;
                 ++i)	//shift kernel pointers to new line
            {
                cachePointers[i] = (cachePointers[i] + cacheWidth
                                                       * (y - previousY)) % cache.length;
            }
            previousY = y;

            int yStartReading = y == 0 ? 0 : y + kernel.kHeight / 2;
            for (int yNew = yStartReading;
                 yNew <= y + kernel.kHeight / 2;
                 ++yNew)
            { //only 1 line except at start
                ReadLineToCacheOrPad(mask, dims, cache, cacheWidth,
                                     kernel.kRadius, kernel.kHeight, yNew);
            }

            int maskP = y * dims.width;

            for (int x = 0;
                 x < dims.width;
                 ++x, ++maskP)
            {
                int total = 0;
                int total_2 = 0;
                int sobelCounter = -1;

                for (int kk = 0;
                     kk < cachePointers.length;
                     ++kk)
                {
                    for (int p = cachePointers[kk++] + x;
                         p <= cachePointers[kk] + x;
                         ++p)
                    {
                        if (filterType == FilterType.Sum)
                        {
                            if ((cache[p] & 0xff) > 0)
                            {
                                total += 1;
                            }
                        }
                        else if (filterType == FilterType.Sobel)
                        {
                            total += ((cache[p] & 0xff) * sobelKernel[++sobelCounter]);
                            total_2 += ((cache[p] & 0xff) * sobelKernel[++sobelCounter]);
                        }
                    }
                }

                if (filterType == FilterType.Sum)
                {
                    if (total >= ForegroundMinCount)
                    {
                        mask[maskP] = (byte) 255;
                    }
                    else
                    {
                        mask[maskP] = (byte) 0;
                    }
                }
                else if (filterType == FilterType.Sobel)
                {
                    total *= total;
                    total_2 *= total_2;
                    
                    if ((total | total_2) > 0)
                    {
                        mask[maskP] = (byte) 255;
                    }
                    else
                    {
                        mask[maskP] = (byte) 0;
                    }
                }
            }
        }
    }

    /**
     * Read a line into the cache (including padding in x). If y>=height,
     * instead of reading new data, it duplicates the line y=height-1. If y==0,
     * it also creates the data for y<0, as far as necessary, thus filling the
     * cache with more than one line (padding by duplicating the y=0 row).
     */
    private static
            void ReadLineToCacheOrPad(byte[] mask, Dims dims, byte[] cache,
                                      int cacheWidth, int pad, int kHeight, int y)
    {
        int lineInCache = y % kHeight;
        if (y < dims.height)
        {
            ReadLineToCache(mask, y * dims.width, dims.width, cache,
                            lineInCache * cacheWidth, pad);
            if (y == 0)
            {
                for (int prevY = -(kHeight / 2);
                     prevY < 0;
                     ++prevY)
                {	//for y<0, pad with y=0 border pixels
                    int prevLineInCache = kHeight + prevY;
                    System.arraycopy(cache, 0, cache, prevLineInCache
                                                      * cacheWidth, cacheWidth);
                }
            }
        }
        else
        {
            System.arraycopy(cache, cacheWidth * ((dims.height - 1) % kHeight),
                             cache, lineInCache * cacheWidth, cacheWidth);
        }
    }

    /**
     * Read a line into the cache. Pad with edge pixels in x if necessary
     */
    private static
            void ReadLineToCache(byte[] mask, int pixelLineP, int width,
                                 byte[] cache, int cacheLineP, int pad)
    {
        System.arraycopy(mask, pixelLineP, cache, cacheLineP + pad, width);

        for (int cp = cacheLineP;
             cp < cacheLineP + pad;
             ++cp)
        {
            cache[cp] = cache[cacheLineP + pad];
        }
        for (int cp = cacheLineP + pad + width;
             cp < cacheLineP + (2 * pad) + width;
             ++cp)
        {
            cache[cp] = cache[cacheLineP + pad + width - 1];
        }
    }

    static
            void InvertMask(byte[] mask)
    {
        for (int index = 0;
             index < mask.length;
             ++index)
        {
            mask[index] = (byte) (((byte) 255) - mask[index]);
        }
    }

    static
            void Watershed(byte[] mask, Dims dims, WatershedParams params)
    {
        float[] floatEdm = new float[mask.length];

        for (int index = 0;
             index < mask.length;
             ++index)
        {
            if (mask[index] != 0)
            {
                floatEdm[index] = Float.MAX_VALUE;
            }
        }

        int[][] pointBufs = new int[2][dims.width];  //two buffers for two passes; low short contains x, high short y

        // pass 1 & 2: increasing y
        for (int x = 0;
             x < dims.width;
             ++x)
        {
            pointBufs[0][x] = -1;
            pointBufs[1][x] = -1;
        }
        for (int y = 0;
             y < dims.height;
             ++y)
        {
            EdmLine(mask, floatEdm, pointBufs, dims.width, y * dims.width, y);
        }
        //pass 3 & 4: decreasing y
        for (int x = 0;
             x < dims.width;
             ++x)
        {
            pointBufs[0][x] = -1;
            pointBufs[1][x] = -1;
        }
        for (int y = dims.height - 1;
             y >= 0;
             --y)
        {
            EdmLine(mask, floatEdm, pointBufs, dims.width, y * dims.width, y);
        }

        pointBufs = null;

        for (int index = 0;
             index < floatEdm.length;
             ++index)
        {
            floatEdm[index] = (float) Math.sqrt(floatEdm[index]);
        }

        FindMaxima(mask, floatEdm, dims, params);
    }

    // Handle a line; two passes: left-to-right and right-to-left
    static private
            void EdmLine(byte[] bPixels, float[] fPixels, int[][] pointBufs, int width,
                         int offset, int y)
    {
        int[] points = pointBufs[0];        // the buffer for the left-to-right pass
        int pPrev = -1;
        int pDiag = -1;               // point at (-/+1, -/+1) to current one (-1,-1 in the first pass)
        int pNextDiag;
        int distSqr = Integer.MAX_VALUE;    // this value is used only if edges are not background
        for (int x = 0;
             x < width;
             ++x, ++offset)
        {
            pNextDiag = points[x];
            if (bPixels[offset] == 0)
            {
                points[x] = x | y << 16;      // remember coordinates as a candidate for nearest background point
            }
            else
            {                        // foreground pixel:
                float dist2 = MinDist2(points, pPrev, pDiag, x, y, distSqr);
                if (fPixels[offset] > dist2)
                {
                    fPixels[offset] = dist2;
                }
            }
            pPrev = points[x];
            pDiag = pNextDiag;
        }
        --offset; //now points to the last pixel in the line
        points = pointBufs[1];              // the buffer for the right-to-left pass. Low short contains x, high short y
        pPrev = -1;
        pDiag = -1;
        for (int x = width - 1;
             x >= 0;
             --x, --offset)
        {
            pNextDiag = points[x];
            if (bPixels[offset] == 0)
            {
                points[x] = x | y << 16;      // remember coordinates as a candidate for nearest background point
            }
            else
            {
                float dist2 = MinDist2(points, pPrev, pDiag, x, y, distSqr);
                if (fPixels[offset] > dist2)
                {
                    fPixels[offset] = dist2;
                }
            }
            pPrev = points[x];
            pDiag = pNextDiag;
        }
    }

    // Calculates minimum distance^2 of x,y from the following three points:
    //  - points[x] (nearest point found for previous line, same x)
    //  - pPrev (nearest point found for same line, previous x), and
    //  - pDiag (nearest point found for diagonal, i.e., previous line, previous x)
    // Sets array element points[x] to the coordinates of the point having the minimum distance to x,y
    // If the distSqr parameter is lower than the distance^2, then distSqr is used
    // Returns to the minimum distance^2 obtained
    static private
            float MinDist2(int[] points, int pPrev, int pDiag,
                           int x, int y, int distSqr)
    {
        int p0 = points[x];              // the nearest background point for the same x in the previous line
        int nearestPoint = p0;
        if (p0 != -1)
        {
            int x0 = p0 & 0xffff;
            int y0 = (p0 >> 16) & 0xffff;
            int dist1Sqr = (x - x0) * (x - x0) + (y - y0) * (y - y0);
            if (dist1Sqr < distSqr)
            {
                distSqr = dist1Sqr;
            }
        }
        if (pDiag != p0 && pDiag != -1)
        {
            int x1 = pDiag & 0xffff;
            int y1 = (pDiag >> 16) & 0xffff;
            int dist1Sqr = (x - x1) * (x - x1) + (y - y1) * (y - y1);
            if (dist1Sqr < distSqr)
            {
                nearestPoint = pDiag;
                distSqr = dist1Sqr;
            }
        }
        if (pPrev != pDiag && pPrev != -1)
        {
            int x1 = pPrev & 0xffff;
            int y1 = (pPrev >> 16) & 0xffff;
            int dist1Sqr = (x - x1) * (x - x1) + (y - y1) * (y - y1);
            if (dist1Sqr < distSqr)
            {
                nearestPoint = pPrev;
                distSqr = dist1Sqr;
            }
        }
        points[x] = nearestPoint;
        return (float) distSqr;
    }

    public static
            class WatershedParams
    {

        int intEncodeXMask;
        int intEncodeYMask;
        int intEncodeShift;
        int[] dirOffset;
        int[] fateTable;

        public
                WatershedParams(int width)
        {
            int shift = 0, mult = 1;
            do
            {
                ++shift;
                mult *= 2;
            }
            while (mult < width);
            intEncodeXMask = mult - 1;
            intEncodeYMask = ~intEncodeXMask;
            intEncodeShift = shift;
            dirOffset = new int[]
            {
                -width, -width + 1, +1, +width + 1, +width, +width - 1, -1, -width - 1
            };

            fateTable = MakeFateTable();
        }

    }

    static private
            void FindMaxima(byte[] mask, float[] edm, Dims dims,
                            WatershedParams params)
    {
        byte[] types = new byte[mask.length];
        float globalMin = Float.MAX_VALUE;
        float globalMax = -Float.MAX_VALUE;
        for (int index = 0;
             index < edm.length;
             ++index)
        {
            float v = edm[index];
            if (globalMin > v)
            {
                globalMin = v;
            }
            if (globalMax < v)
            {
                globalMax = v;
            }

        }

        long[] maxPoints = GetSortedMaxPoints(edm, dims, types, globalMin, globalMax, params.dirOffset);
        AnalyzeAndMarkMaxima(edm, dims, types, maxPoints, globalMin, params.dirOffset);

        // Segmentation required, convert to 8bit (also for 8-bit images, since the calibration
        // may have a negative slope). outIp has background 0, maximum areas 255
        Make8bit(mask, edm, dims, types, globalMax);
        CleanupMaxima(mask, dims, types, maxPoints, params.dirOffset);     //eliminate all the small maxima (i.e. those outside MAX_AREA)
        WatershedSegment(mask, dims, params);
        for (int i = 0;
             i < mask.length;
             ++i)
        {
            if ((mask[i] & 255) < 255)
            {
                mask[i] = (byte) 0;
            }
        }                //levels to binary image
    } //

    static private
            long[] GetSortedMaxPoints(float[] edm, Dims dims, byte[] types, float globalMin,
                                      float globalMax, int[] dirOffset)
    {
        int nMax = 0;
        for (int y = 0;
             y < dims.height;
             ++y)
        {
            for (int x = 0, index = x + y * dims.width;
                 x < dims.width;
                 ++x, ++index)
            {
                float v = edm[index];
                float vTrue = TrueEdmHeight(x, y, dims, edm, dirOffset);
                if (v == globalMin)
                {
                    continue;
                }
                boolean isMax = true;
                /* check wheter we have a local maximum.
                 Note: For an EDM, we need all maxima: those of the EDM-corrected values
                 (needed by findMaxima) and those of the raw values (needed by cleanupMaxima) */
                boolean isInner = (y != 0 && y != dims.height - 1)
                                  && (x != 0 && x != dims.width - 1); //not necessary, but faster than isWithin
                for (int d = 0;
                     d < 8;
                     ++d)
                {                         // compare with the 8 neighbor pixels
                    if (isInner || IsWithin(x, y, d, dims))
                    {
                        float vNeighbor = edm[x + DIR_X_OFFSET[d]
                                              + ((y + DIR_Y_OFFSET[d]) * dims.width)];
                        float vNeighborTrue = TrueEdmHeight(x + DIR_X_OFFSET[d], y + DIR_Y_OFFSET[d],
                                                            dims, edm, dirOffset);
                        if (vNeighbor > v && vNeighborTrue > vTrue)
                        {
                            isMax = false;
                            break;
                        }
                    }
                }
                if (isMax)
                {
                    types[index] = MAXIMUM;
                    nMax++;
                }
            } // for x
        } // for y

        float vFactor = (float) (2e9 / (globalMax - globalMin)); //for converting float values into a 32-bit int
        long[] maxPoints = new long[nMax];                  //value (int) is in the upper 32 bit, pixel offset in the lower
        int iMax = 0;
        for (int y = 0;
             y < dims.height;
             ++y)           //enter all maxima into an array
        {
            for (int x = 0, index = x + y * dims.width;
                 x < dims.width;
                 ++x, ++index)
            {
                if (types[index] == MAXIMUM)
                {
                    float fValue = TrueEdmHeight(x, y, dims, edm, dirOffset);
                    int iValue = (int) ((fValue - globalMin) * vFactor); //32-bit int, linear function of float value
                    maxPoints[iMax++] = (long) iValue << 32 | index;
                }
            }
        }
        Arrays.sort(maxPoints); //sort the maxima by value
        return maxPoints;
    } //getSortedMaxPoints

    /**
     * Get estimated "true" height of a maximum or saddle point of a Euclidian
     * Distance Map. This is needed since the point sampled is not necessarily
     * at the highest position. For simplicity, we don't care about the Sqrt(5)
     * distance here although this would be more accurate
     *
     * @param x x-position of the point
     * @param y y-position of the point
     * @param ip the EDM (FloatProcessor)
     * @return estimated height
     */
    static private
            float TrueEdmHeight(int x, int y, Dims dims, float[] pixels, int[] dirOffset)
    {
        int xmax = dims.width - 1;
        int ymax = dims.height - 1;
        int offset = x + y * dims.width;
        float v = pixels[offset];
        if (x == 0 || y == 0 || x == xmax || y == ymax || v == 0)
        {
            return v;  //don't recalculate for edge pixels or background
        }
        else
        {
            float trueH = v + 0.5f * SQRT2;   //true height can never by higher than this
            boolean ridgeOrMax = false;
            for (int d = 0;
                 d < 4;
                 ++d)
            {               //for all directions halfway around:
                int d2 = (d + 4) % 8;  //get the opposite direction and neighbors
                float v1 = pixels[offset + dirOffset[d]];
                float v2 = pixels[offset + dirOffset[d2]];
                float h;
                if (v >= v1 && v >= v2)
                {
                    ridgeOrMax = true;
                    h = (v1 + v2) / 2;
                }
                else
                {
                    h = Math.min(v1, v2);
                }
                h += (d % 2 == 0) ? 1 : SQRT2;  //in diagonal directions, distance is sqrt2
                if (trueH > h)
                {
                    trueH = h;
                }
            }
            if (!ridgeOrMax)
            {
                trueH = v;
            }
            return trueH;
        }
    }

    static private
            boolean IsWithin(int x, int y, int direction, Dims dims)
    {
        int xmax = dims.width - 1;
        int ymax = dims.height - 1;
        switch (direction)
        {
            case 0:
                return (y > 0);
            case 1:
                return (x < xmax && y > 0);
            case 2:
                return (x < xmax);
            case 3:
                return (x < xmax && y < ymax);
            case 4:
                return (y < ymax);
            case 5:
                return (x > 0 && y < ymax);
            case 6:
                return (x > 0);
            case 7:
                return (x > 0 && y > 0);
        }
        return false;   //to make the compiler happy :-)
    } // isWithin

    static private
            void AnalyzeAndMarkMaxima(float[] edm, Dims dims, byte[] types, long[] maxPoints,
                                      float globalMin, int[] dirOffset)
    {
        float maxSortingError = 0.5f * 1.1f * SQRT2;
        double tolerance = 0.5;

        int nMax = maxPoints.length;
        int[] pList = new int[edm.length];       //here we enter points starting from a maximum

        for (int iMax = nMax - 1;
             iMax >= 0;
             --iMax)
        {    //process all maxima now, starting from the highest
            int offset0 = (int) maxPoints[iMax];     //type cast gets 32 lower bits, where pixel index is encoded
            //int offset0 = maxPoints[iMax].offset;
            if ((types[offset0] & PROCESSED) != 0)      //this maximum has been reached from another one, skip it
            {
                continue;
            }
            //we create a list of connected points and start the list at the current maximum
            int x0 = offset0 % dims.width;
            int y0 = offset0 / dims.width;
            float v0 = TrueEdmHeight(x0, y0, dims, edm, dirOffset);
            boolean sortingError;
            do
            {                                    //repeat if we have encountered a sortingError
                pList[0] = offset0;
                types[offset0] |= (EQUAL | LISTED);   //mark first point as equal height (to itself) and listed
                int listLen = 1;                    //number of elements in the list
                int listI = 0;                      //index of current element in the list
                sortingError = false;       //if sorting was inaccurate: a higher maximum was not handled so far
                boolean maxPossible = true;         //it may be a true maximum
                double xEqual = x0;                 //for creating a single point: determine average over the
                double yEqual = y0;                 //  coordinates of contiguous equal-height points
                int nEqual = 1;                     //counts xEqual/yEqual points that we use for averaging
                do
                {                                //while neigbor list is not fully processed (to listLen)
                    int offset = pList[listI];
                    int x = offset % dims.width;
                    int y = offset / dims.width;
                    boolean isInner = (y != 0 && y != dims.height - 1)
                                      && (x != 0 && x != dims.width - 1); //not necessary, but faster than isWithin
                    for (int d = 0;
                         d < 8;
                         ++d)
                    {       //analyze all neighbors (in 8 directions) at the same level
                        int offset2 = offset + dirOffset[d];
                        if ((isInner || IsWithin(x, y, d, dims))
                            && (types[offset2] & LISTED) == 0)
                        {
                            if (edm[offset2] <= 0)
                            {
                                continue;   //ignore the background (non-particles)
                            }
                            if ((types[offset2] & PROCESSED) != 0)
                            {
                                maxPossible = false; //we have reached a point processed previously, thus it is no maximum now
                                break;
                            }
                            int x2 = x + DIR_X_OFFSET[d];
                            int y2 = y + DIR_Y_OFFSET[d];
                            float v2 = TrueEdmHeight(x2, y2, dims, edm, dirOffset);
                            if (v2 > v0 + maxSortingError)
                            {
                                maxPossible = false;    //we have reached a higher point, thus it is no maximum
                                break;
                            }
                            else if (v2 >= v0 - (float) tolerance)
                            {
                                if (v2 > v0)
                                {          //maybe this point should have been treated earlier
                                    sortingError = true;
                                    offset0 = offset2;
                                    v0 = v2;
                                    x0 = x2;
                                    y0 = y2;

                                }
                                pList[listLen] = offset2;
                                listLen++;              //we have found a new point within the tolerance
                                types[offset2] |= LISTED;
                                if (v2 == v0)
                                {           //prepare finding center of equal points (in case single point needed)
                                    types[offset2] |= EQUAL;
                                    xEqual += x2;
                                    yEqual += y2;
                                    ++nEqual;
                                }
                            }
                        } // if isWithin & not LISTED
                    } // for directions d
                    ++listI;
                }
                while (listI < listLen);

                if (sortingError)
                {				  //if x0,y0 was not the true maximum but we have reached a higher one
                    for (listI = 0;
                         listI < listLen;
                         ++listI)
                    {
                        types[pList[listI]] = 0;	//reset all points encountered, then retry
                    }
                }
                else
                {
                    int resetMask = ~(maxPossible ? LISTED : (LISTED | EQUAL));
                    xEqual /= nEqual;
                    yEqual /= nEqual;
                    double minDist2 = 1e20;
                    int nearestI = 0;
                    for (listI = 0;
                         listI < listLen;
                         ++listI)
                    {
                        int offset = pList[listI];
                        int x = offset % dims.width;
                        int y = offset / dims.width;
                        types[offset] &= resetMask;		//reset attributes no longer needed
                        types[offset] |= PROCESSED;		//mark as processed
                        if (maxPossible)
                        {
                            types[offset] |= MAX_AREA;
                            if ((types[offset] & EQUAL) != 0)
                            {
                                double dist2 = (xEqual - x) * (double) (xEqual - x)
                                               + (yEqual - y) * (double) (yEqual - y);
                                if (dist2 < minDist2)
                                {
                                    minDist2 = dist2;	//this could be the best "single maximum" point
                                    nearestI = listI;
                                }
                            }
                        }
                    } // for listI
                    if (maxPossible)
                    {
                        int offset = pList[nearestI];
                        types[offset] |= MAX_POINT;
                    }
                } //if !sortingError
            }
            while (sortingError);	//redo if we have encountered a higher maximum: handle it now.
        } // for all maxima iMax
    } //void analyzeAndMarkMaxima

    static private
            void Make8bit(byte[] pixels, float[] edm,
                          Dims dims, byte[] types,
                          float globalMax)
    {
        double threshold = 0.5;
        double minValue = 1.;

        double offset = minValue - (globalMax - minValue)
                                   * (1. / 253 / 2 - 1e-6); //everything above minValue should become >(byte)0
        double factor = 253 / (globalMax - minValue);

        if (factor > 1)
        {
            factor = 1;   // with EDM, no better resolution
        }

        long v;
        for (int y = 0, index = 0;
             y < dims.height;
             ++y)
        {
            for (int x = 0;
                 x < dims.width;
                 ++x, ++index)
            {
                float rawValue = edm[index];
                if (rawValue < threshold)
                {
                    pixels[index] = (byte) 0;
                }
                else if ((types[index] & MAX_AREA) != 0)
                {
                    pixels[index] = (byte) 255;  //prepare watershed by setting "true" maxima+surroundings to 255
                }
                else
                {
                    v = 1 + Math.round((rawValue - offset) * factor);
                    if (v < 1)
                    {
                        pixels[index] = (byte) 1;
                    }
                    else if (v <= 254)
                    {
                        pixels[index] = (byte) (v & 255);
                    }
                    else
                    {
                        pixels[index] = (byte) 254;
                    }
                }
            }
        }
    } // byteProcessor make8bit

    private static
            void CleanupMaxima(byte[] pixels, Dims dims, byte[] types,
                               long[] maxPoints, int[] dirOffset)
    {
        int nMax = maxPoints.length;
        int[] pList = new int[pixels.length];
        for (int iMax = nMax - 1;
             iMax >= 0;
             --iMax)
        {
            int offset0 = (int) maxPoints[iMax];     //type cast gets lower 32 bits where pixel offset is encoded
            if ((types[offset0] & (MAX_AREA | ELIMINATED)) != 0)
            {
                continue;
            }
            int level = pixels[offset0] & 255;
            int loLevel = level + 1;
            pList[0] = offset0;                     //we start the list at the current maximum

            types[offset0] |= LISTED;               //mark first point as listed
            int listLen = 1;                        //number of elements in the list
            int lastLen = 1;
            int listI;                          //index of current element in the list
            boolean saddleFound = false;
            while (!saddleFound && loLevel > 0)
            {
                loLevel--;
                lastLen = listLen;                  //remember end of list for previous level
                listI = 0;                          //in each level, start analyzing the neighbors of all pixels
                do
                {                                //for all pixels listed so far
                    int offset = pList[listI];
                    int x = offset % dims.width;
                    int y = offset / dims.width;
                    boolean isInner = (y != 0 && y != dims.height - 1)
                                      && (x != 0 && x != dims.width - 1); //not necessary, but faster than isWithin
                    for (int d = 0;
                         d < 8;
                         ++d)
                    {       //analyze all neighbors (in 8 directions) at the same level
                        int offset2 = offset + dirOffset[d];
                        if ((isInner || IsWithin(x, y, d, dims))
                            && (types[offset2] & LISTED) == 0)
                        {
                            if ((types[offset2] & MAX_AREA) != 0
                                || (((types[offset2] & ELIMINATED) != 0)
                                    && (pixels[offset2] & 255) >= loLevel))
                            {
                                saddleFound = true; //we have reached a point touching a "true" maximum...
                                break;              //...or a level not lower, but touching a "true" maximum
                            }
                            else if ((pixels[offset2] & 255) >= loLevel
                                     && (types[offset2] & ELIMINATED) == 0)
                            {
                                pList[listLen] = offset2;
                                ++listLen;          //we have found a new point to be processed
                                types[offset2] |= LISTED;
                            }
                        } // if isWithin & not LISTED
                    } // for directions d
                    if (saddleFound)
                    {
                        break;         //no reason to search any further
                    }
                    listI++;
                }
                while (listI < listLen);
            } // while !levelFound && loLevel>=0
            for (listI = 0;
                 listI < listLen;
                 ++listI)   //reset attribute since we may come to this place again
            {
                types[pList[listI]] &= ~LISTED;
            }
            for (listI = 0;
                 listI < lastLen;
                 ++listI)
            { //for all points higher than the level of the saddle point
                int offset = pList[listI];
                pixels[offset] = (byte) loLevel;     //set pixel value to the level of the saddle point
                types[offset] |= ELIMINATED;        //mark as processed: there can't be a local maximum in this area
            }
        } // for all maxima iMax
    } // void cleanupMaxima

    static private
            void WatershedSegment(byte[] pixels, Dims dims,
                                  WatershedParams params)
    {
        // Create an array with the coordinates of all points between value 1 and 254
        // This method, suggested by Stein Roervik (stein_at_kjemi-dot-unit-dot-no),
        // greatly speeds up the watershed segmentation routine.
        int[] histogram = MakeHistogram(pixels);
        int arraySize = pixels.length - histogram[0] - histogram[255];
        int[] coordinates = new int[arraySize];    //from pixel coordinates, low bits x, high bits y
        int highestValue = 0;
        int maxBinSize = 0;
        int offset = 0;
        int[] levelStart = new int[256];
        for (int v = 1;
             v < 255;
             ++v)
        {
            levelStart[v] = offset;
            offset += histogram[v];
            if (histogram[v] > 0)
            {
                highestValue = v;
            }
            if (histogram[v] > maxBinSize)
            {
                maxBinSize = histogram[v];
            }
        }
        int[] levelOffset = new int[highestValue + 1];
        for (int y = 0, index = 0;
             y < dims.height;
             ++y)
        {
            for (int x = 0;
                 x < dims.width;
                 ++x, ++index)
            {
                int v = pixels[index] & 255;
                if (v > 0 && v < 255)
                {
                    offset = levelStart[v] + levelOffset[v];
                    coordinates[offset] = x | y << params.intEncodeShift;
                    levelOffset[v]++;
                }
            } //for x
        } //for y
        // Create an array of the points (pixel offsets) that we set to 255 in one pass.
        // If we remember this list we need not create a snapshot of the ImageProcessor.
        int[] setPointList = new int[Math.min(maxBinSize, (pixels.length + 2) / 3)];
        // now do the segmentation, starting at the highest level and working down.
        // At each level, dilate the particle (set pixels to 255), constrained to pixels
        // whose values are at that level and also constrained (by the fateTable)
        // to prevent features from merging.
        final
                int[] directionSequence = new int[]
                {
                    7, 3, 1, 5, 0, 4, 2, 6
                }; // diagonal directions first
        for (int level = highestValue;
             level >= 1;
             --level)
        {
            int remaining = histogram[level];  //number of points in the level that have not been processed
            int idle = 0;
            while (remaining > 0 && idle < 8)
            {
                int dIndex = 0;
                do
                {                        // expand each level in 8 directions
                    int n = ProcessLevel(directionSequence[dIndex % 8], pixels, dims,
                                         levelStart[level], remaining, coordinates, setPointList,
                                         params);

                    remaining -= n;         // number of points processed
                    if (n > 0)
                    {
                        idle = 0;    // nothing processed in this direction?
                    }
                    ++dIndex;
                }
                while (remaining > 0 && idle++ < 8);
            }
            if (remaining > 0 && level > 1)
            {   // any pixels that we have not reached?
                int nextLevel = level;      // find the next level to process
                do
                {
                    --nextLevel;
                }
                while (nextLevel > 1 && histogram[nextLevel] == 0);
                // in principle we should add all unprocessed pixels of this level to the
                // tasklist of the next level. This would make it very slow for some images,
                // however. Thus we only add the pixels if they are at the border (of the
                // image or a thresholded area) and correct unprocessed pixels at the very
                // end by CleanupExtraLines
                if (nextLevel > 0)
                {
                    int newNextLevelEnd = levelStart[nextLevel] + histogram[nextLevel];
                    for (int index = 0, p = levelStart[level];
                         index < remaining;
                         ++index, ++p)
                    {
                        int xy = coordinates[p];
                        int x = xy & params.intEncodeXMask;
                        int y = (xy & params.intEncodeYMask) >> params.intEncodeShift;
                        int pOffset = x + y * dims.width;

                        boolean addToNext = false;
                        if (x == 0 || y == 0 || x == dims.width - 1 || y == dims.height - 1)
                        {
                            addToNext = true;           //image border
                        }
                        else
                        {
                            for (int d = 0; d < 8; d++)
                            {
                                if (IsWithin(x, y, d, dims) && pixels[pOffset
                                                                      + params.dirOffset[d]] == 0)
                                {
                                    addToNext = true;       //border of area below threshold
                                    break;
                                }
                            }
                        }
                        if (addToNext)
                        {
                            coordinates[newNextLevelEnd++] = xy;
                        }
                    }

                    histogram[nextLevel] = newNextLevelEnd - levelStart[nextLevel];
                }
            }
        }
    } // boolean watershedSegment

    static private
            int ProcessLevel(int pass, byte[] pixels, Dims dims,
                             int levelStart, int levelNPoints,
                             int[] coordinates, int[] setPointList,
                             WatershedParams params)
    {
        int xmax = dims.width - 1;
        int ymax = dims.height - 1;
        int nChanged = 0;
        int nUnchanged = 0;
        for (int i = 0, p = levelStart;
             i < levelNPoints;
             ++i, ++p)
        {
            int xy = coordinates[p];
            int x = xy & params.intEncodeXMask;
            int y = (xy & params.intEncodeYMask) >> params.intEncodeShift;
            int offset = x + y * dims.width;
            int index = 0;      //neighborhood pixel ocupation: index in fateTable
            if (y > 0 && (pixels[offset - dims.width] & 255) == 255)
            {
                index ^= 1;
            }
            if (x < xmax && y > 0 && (pixels[offset
                                             - dims.width + 1] & 255) == 255)
            {
                index ^= 2;
            }
            if (x < xmax && (pixels[offset + 1] & 255) == 255)
            {
                index ^= 4;
            }
            if (x < xmax && y < ymax && (pixels[offset
                                                + dims.width + 1] & 255) == 255)
            {
                index ^= 8;
            }
            if (y < ymax && (pixels[offset + dims.width] & 255) == 255)
            {
                index ^= 16;
            }
            if (x > 0 && y < ymax && (pixels[offset
                                             + dims.width - 1] & 255) == 255)
            {
                index ^= 32;
            }
            if (x > 0 && (pixels[offset - 1] & 255) == 255)
            {
                index ^= 64;
            }
            if (x > 0 && y > 0 && (pixels[offset
                                          - dims.width - 1] & 255) == 255)
            {
                index ^= 128;
            }
            int mask = 1 << pass;
            if ((params.fateTable[index] & mask) == mask)
            {
                setPointList[nChanged++] = offset;  //remember to set pixel to 255
            }
            else
            {
                coordinates[levelStart + (nUnchanged++)] = xy; //keep this pixel for future passes
            }
        } // for pixel i
        for (int i = 0;
             i < nChanged;
             ++i)
        {
            pixels[setPointList[i]] = (byte) 255;
        }
        return nChanged;
    } //processLevel

    /**
     * Creates the lookup table used by the watershed function for dilating the
     * particles. The algorithm allows dilation in both straight and diagonal
     * directions. There is an entry in the table for each possible 3x3
     * neighborhood: x-1 x x+1 y-1 128 1 2 y 64 pxl_unset_yet 4 y+1 32 16 8 (to
     * find throws entry, sum up the numbers of the neighboring pixels set; e.g.
     * entry 6=2+4 if only the pixels (x,y-1) and (x+1, y-1) are set. A pixel is
     * added on the 1st pass if bit 0 (2^0 = 1) is set, on the 2nd pass if bit 1
     * (2^1 = 2) is set, etc. pass gives the direction of rotation, with 0 = to
     * top left (x--,y--), 1 to top, and clockwise up to 7 = to the left (x--).
     * E.g. 4 = add on 3rd pass, 3 = add on either 1st or 2nd pass.
     */
    static private
            int[] MakeFateTable()
    {
        int[] table = new int[256];
        boolean[] isSet = new boolean[8];
        for (int item = 0;
             item < 256;
             ++item)
        {        //dissect into pixels
            for (int i = 0, mask = 1;
                 i < 8;
                 ++i)
            {
                isSet[i] = (item & mask) == mask;
                mask *= 2;
            }
            for (int i = 0, mask = 1;
                 i < 8;
                 ++i)
            {       //we dilate in the direction opposite to the direction of the existing neighbors
                if (isSet[(i + 4) % 8])
                {
                    table[item] |= mask;
                }
                mask *= 2;
            }
            for (int i = 0;
                 i < 8;
                 i += 2)                //if side pixels are set, for counting transitions it is as good as if the adjacent edges were also set
            {
                if (isSet[i])
                {
                    isSet[(i + 1) % 8] = true;
                    isSet[(i + 7) % 8] = true;
                }
            }
            int transitions = 0;
            for (int i = 0;
                 i < 8;
                 ++i)
            {
                if (isSet[i] != isSet[(i + 1) % 8])
                {
                    transitions++;
                }
            }
            if (transitions >= 4)
            {                   //if neighbors contain more than one region, dilation ito this pixel is forbidden
                table[item] = 0;
            }
        }
        return table;
    } // int[] makeFateTable

    static
            byte[] AddMasks(byte[] in1, byte[] in2)
    {
        byte[] out = new byte[in1.length];
        for (int index = 0;
             index < in1.length;
             ++index)
        {
            out[index] = (byte) (in1[index] | in2[index]);
        }
        return out;
    }

    static
            byte[] SubtractMasks(byte[] in1, byte[] in2)
    {
        byte[] out = new byte[in1.length];
        for (int index = 0;
             index < in1.length;
             ++index)
        {
            out[index] = (byte) ((in1[index] ^ in2[index])
                                 & in1[index]);
        }
        return out;
    }

    static
            byte[] CopyMask(byte[] in)
    {
        byte[] out = new byte[in.length];
        System.arraycopy(in, 0, out, 0, in.length);
        return out;
    }

    final static private
            float SQRT2 = 1.4142135624f;

    final static private
            byte MAXIMUM = (byte) 1;            // marks local maxima (irrespective of noise tolerance)
    final static private
            byte LISTED = (byte) 2;             // marks points currently in the list
    final static private
            byte PROCESSED = (byte) 4;          // marks points processed previously
    final static private
            byte MAX_AREA = (byte) 8;           // marks areas near a maximum, within the tolerance
    final static private
            byte EQUAL = (byte) 16;             // marks contigous maximum points of equal level
    final static private
            byte MAX_POINT = (byte) 32;         // marks a single point standing for a maximum
    final static private
            byte ELIMINATED = (byte) 64;        // marks maxima that have been eliminated before watershed

    final static private
            int[] DIR_X_OFFSET = new int[]
            {
                0, 1, 1, 1, 0, -1, -1, -1
            };
    final static private
            int[] DIR_Y_OFFSET = new int[]
            {
                -1, -1, 0, 1, 1, 1, 0, -1
            };
}
