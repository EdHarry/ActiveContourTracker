/*
 * Implamentation of ImageJ's BackGroundSubtractor
* 
* ORIGIANL HEADER:
//* Implements ImageJ's Subtract Background command. Based on the concept of the
//rolling ball algorithm described in Stanley Sternberg's article, "Biomedical Image
//Processing", IEEE Computer, January 1983.
//
//Imagine that the 2D grayscale image has a third (height) dimension by the image
//value at every point in the image, creating a surface. A ball of given radius is
//rolled over the bottom side of this surface; the hull of the volume reachable by
//the ball is the background.
//
//With "Sliding Parabvoloid", the rolling ball is replaced by a sliding paraboloid
//of rotation with the same curvature at its apex as a ball of a given radius.
//A paraboloid has the advantage that suitable paraboloids can be found for any image
//values, even if the pixel values are much larger than a typical object size (in pixels).
//The paraboloid of rotation is approximated as parabolae in 4 directions: x, y and
//the two 45-degree directions. Lines of the image in these directions are processed
//by sliding a parabola against them. Obtaining the hull needs the parabola for a
//given direction to be applied multiple times (after doing the other directions);
//in this respect the current code is a compromise between accuracy and speed.
//
//For noise rejection, with the sliding paraboloid algorithm, a 3x3 maximum of the
//background is applied. With both, rolling ball and sliding paraboloid,
//the image used for calculating the background is slightly smoothened (3x3 average).
//This can result in negative values after background subtraction. This preprocessing
//can be disabled.
//
//In the sliding paraboloid algorithm, additional code has been added to avoid
//subtracting corner objects as a background (note that a paraboloid or ball would
//always touch the 4 corner pixels and thus make them background pixels).
//This code assumes that corner particles reach less than 1/4 of the image size
//into the image.
//
//Rolling ball code based on the NIH Image Pascal version by Michael Castle and Janice 
//Keller of the University of Michigan Mental Health Research Institute.
//Sliding Paraboloid by Michael Schmid, 2007.
//
//Version 10-Jan-2008
 */
package activecontourtracker;

/**
 *
 * @author edwardharry
 */
public
        class RollingBallBackgroundSubtractor
{

    static
            void SubtractBackground(short[] pixels, Dims dims,
                                    RollingBall ball)
    {
        float[] floats = new float[pixels.length];
        for (int index = 0;
             index < floats.length;
             ++index)
        {
            floats[index] = (pixels[index] & 0xffff);
        }

        SubtractBackground(floats, dims, ball);

        for (int index = 0;
             index < floats.length;
             ++index)
        {
            float value = floats[index];
            if (value < 0f)
            {
                value = 0f;
            }
            if (value > 65535f)
            {
                value = 65535f;
            }
            pixels[index] = (short) value;
        }
    }

    static
            void SubtractBackground(float[] pixels, Dims dims,
                                    RollingBall ball)
    {
        Img img = new Img(dims, pixels);
        Img imgBg = img.MakeCopy();

        RollingBallFloatBackground(imgBg, ball);

        for (int index = 0;
             index < img.pixels.length;
             ++index)
        {
            img.pixels[index] -= (imgBg.pixels[index] - 0.5f);
        }
    }

    static private
            class Img
    {

        Dims dims;
        float[] pixels;

        Img(Dims d)
        {
            dims = d;
            pixels = new float[d.width * d.height];
        }

        Img(Dims d, float[] p)
        {
            dims = d;
            pixels = p;
        }

        Img MakeCopy()
        {
            Img newImg = new Img(dims);
            System.arraycopy(pixels, 0, newImg.pixels, 0, pixels.length);
            return newImg;
        }
    }

    /**
     * Create background for a float image by rolling a ball over the image.
     */
    private static
            void RollingBallFloatBackground(Img img, RollingBall ball)
    {
        boolean shrink = ball.shrinkFactor > 1;

        Filter3x3(img);
        Img smallImage = shrink ? ShrinkImage(img, ball.shrinkFactor) : img;
        RollBall(ball, smallImage);
        if (shrink)
        {
            EnlargeImage(smallImage, img, ball.shrinkFactor);
        }
    }

    /**
     * Creates a lower resolution image for ball-rolling.
     */
    private static
            Img ShrinkImage(Img img, int shrinkFactor)
    {
        int width = img.dims.width;
        int height = img.dims.height;
        float[] pixels = img.pixels;

        int sWidth = (width + shrinkFactor - 1) / shrinkFactor;
        int sHeight = (height + shrinkFactor - 1) / shrinkFactor;

        Img smallImg = new Img(new Dims(sWidth, sHeight));
        float[] sPixels = smallImg.pixels;

        float min, thispixel;
        for (int ySmall = 0; ySmall < sHeight; ySmall++)
        {
            for (int xSmall = 0; xSmall < sWidth; xSmall++)
            {
                min = Float.MAX_VALUE;
                for (int j = 0, y = shrinkFactor * ySmall; j < shrinkFactor && y < height; j++, y++)
                {
                    for (int k = 0, x = shrinkFactor * xSmall; k < shrinkFactor && x < width; k++, x++)
                    {
                        thispixel = pixels[x + y * width];
                        if (thispixel < min)
                        {
                            min = thispixel;
                        }
                    }
                }
                sPixels[xSmall + ySmall * sWidth] = min; // each point in small image is minimum of its neighborhood
            }
        }
        return smallImg;
    }

    /**
     * 'Rolls' a filtering object over a (shrunken) image in order to find the
     * image's smooth continuous background. For the purpose of explaining this
     * algorithm, imagine that the 2D grayscale image has a third (height)
     * dimension defined by the intensity value at every point in the image. The
     * center of the filtering object, a patch from the top of a sphere having
     * radius BallRadius, is moved along each scan line of the image so that the
     * patch is tangent to the image at one or more points with every other
     * point on the patch below the corresponding (x,y) point of the image. Any
     * point either on or below the patch during this process is considered part
     * of the background. Shrinking the image before running this procedure is
     * advised for large ball radii because the processing time increases with
     * ball radius^2.
     */
    private static
            void RollBall(RollingBall ball, Img img)
    {
        int width = img.dims.width;
        int height = img.dims.height;
        float[] pixels = img.pixels;

        float[] zBall = ball.data;
        int ballWidth = ball.width;
        int radius = ballWidth / 2;
        float[] cache = new float[width * ballWidth]; //temporarily stores the pixels we work on

        for (int y = -radius; y < height + radius; y++)
        { //for all positions of the ball center:

            int nextLineToWriteInCache = (y + radius) % ballWidth;
            int nextLineToRead = y + radius;        //line of the input not touched yet
            if (nextLineToRead < height)
            {
                System.arraycopy(pixels, nextLineToRead * width, cache, nextLineToWriteInCache * width, width);
                for (int x = 0, p = nextLineToRead * width; x < width; x++, p++)
                {
                    pixels[p] = -Float.MAX_VALUE;   //unprocessed pixels start at minus infinity
                }
            }
            int y0 = y - radius;                      //the first line to see whether the ball touches
            if (y0 < 0)
            {
                y0 = 0;
            }
            int yBall0 = y0 - y + radius;               //y coordinate in the ball corresponding to y0
            int yend = y + radius;                    //the last line to see whether the ball touches
            if (yend >= height)
            {
                yend = height - 1;
            }
            for (int x = -radius; x < width + radius; x++)
            {
                float z = Float.MAX_VALUE;          //the height of the ball (ball is in position x,y)
                int x0 = x - radius;
                if (x0 < 0)
                {
                    x0 = 0;
                }
                int xBall0 = x0 - x + radius;
                int xend = x + radius;
                if (xend >= width)
                {
                    xend = width - 1;
                }
                for (int yp = y0, yBall = yBall0; yp <= yend; yp++, yBall++)
                { //for all points inside the ball
                    int cachePointer = (yp % ballWidth) * width + x0;
                    for (int xp = x0, bp = xBall0 + yBall * ballWidth; xp <= xend; xp++, cachePointer++, bp++)
                    {
                        float zReduced = cache[cachePointer] - zBall[bp];
                        if (z > zReduced)           //does this point imply a greater height?
                        {
                            z = zReduced;
                        }
                    }
                }
                for (int yp = y0, yBall = yBall0; yp <= yend; yp++, yBall++) //raise pixels to ball surface
                {
                    for (int xp = x0, p = xp + yp * width, bp = xBall0 + yBall * ballWidth; xp <= xend; xp++, p++, bp++)
                    {
                        float zMin = z + zBall[bp];
                        if (pixels[p] < zMin)
                        {
                            pixels[p] = zMin;
                        }
                    }
                }
            }
        }
    }

    /**
     * Uses bilinear interpolation to find the points in the full-scale
     * background given the points from the shrunken image background. (At the
     * edges, it is actually extrapolation.)
     */
    private static
            void EnlargeImage(Img smallImg, Img img, int shrinkFactor)
    {
        int width = img.dims.width;
        int height = img.dims.height;
        int smallWidth = smallImg.dims.width;
        int smallHeight = smallImg.dims.height;
        float[] pixels = img.pixels;
        float[] sPixels = smallImg.pixels;
        int[] xSmallIndices = new int[width];         //index of first point in smallImage
        float[] xWeights = new float[width];        //weight of this point
        MakeInterpolationArrays(xSmallIndices, xWeights, width, smallWidth, shrinkFactor);
        int[] ySmallIndices = new int[height];
        float[] yWeights = new float[height];
        MakeInterpolationArrays(ySmallIndices, yWeights, height, smallHeight, shrinkFactor);
        float[] line0 = new float[width];
        float[] line1 = new float[width];
        for (int x = 0; x < width; x++)                 //x-interpolation of the first smallImage line
        {
            line1[x] = sPixels[xSmallIndices[x]] * xWeights[x]
                       + sPixels[xSmallIndices[x] + 1] * (1f - xWeights[x]);
        }
        int ySmallLine0 = -1;                       //line0 corresponds to this y of smallImage
        for (int y = 0; y < height; y++)
        {
            if (ySmallLine0 < ySmallIndices[y])
            {
                float[] swap = line0;               //previous line1 -> line0
                line0 = line1;
                line1 = swap;                       //keep the other array for filling with new data
                ySmallLine0++;
                int sYPointer = (ySmallIndices[y] + 1) * smallWidth; //points to line0 + 1 in smallImage
                for (int x = 0; x < width; x++)         //x-interpolation of the new smallImage line -> line1
                {
                    line1[x] = sPixels[sYPointer + xSmallIndices[x]] * xWeights[x]
                               + sPixels[sYPointer + xSmallIndices[x] + 1] * (1f - xWeights[x]);
                }
            }
            float weight = yWeights[y];
            for (int x = 0, p = y * width; x < width; x++, p++)
            {
                pixels[p] = line0[x] * weight + line1[x] * (1f - weight);
            }
        }
    }

    /**
     * Create arrays of indices and weigths for interpolation.
     * <pre>
     * Example for shrinkFactor = 4:
     * small image pixel number         |       0       |       1       |       2       | ...
     * full image pixel number          | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |11 | ...
     * smallIndex for interpolation(0)  | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 2 | 2 | ...
     * (0) Note: This is smallIndex for the left pixel; for the right pixel used for interpolation
     * it is higher by one
     * </pre>
     */
    private static
            void MakeInterpolationArrays(int[] smallIndices, float[] weights, int length, int smallLength, int shrinkFactor)
    {
        for (int i = 0; i < length; i++)
        {
            int smallIndex = (i - shrinkFactor / 2) / shrinkFactor;
            if (smallIndex >= smallLength - 1)
            {
                smallIndex = smallLength - 2;
            }
            smallIndices[i] = smallIndex;
            float distance = (i + 0.5f) / shrinkFactor - (smallIndex + 0.5f); //distance of pixel centers (in smallImage pixels)
            weights[i] = 1f - distance;
        }
    }

    //   C O M M O N   S E C T I O N   F O R   B O T H   A L G O R I T H M S
    /**
     * Replace the pixels by the mean or maximum in a 3x3 neighborhood. No
     * snapshot is required (less memory needed than e.g., fp.smooth()). When
     * used as maximum filter, it returns the average change of the pixel value
     * by this operation
     */
    private static
            void Filter3x3(Img img)
    {
        int width = img.dims.width;
        int height = img.dims.height;
        float[] pixels = img.pixels;

        for (int y = 0; y < height; y++)
        {
            Filter3(pixels, width, y * width, 1);
        }
        for (int x = 0; x < width; x++)
        {
            Filter3(pixels, height, x, width);
        }
    }

    /**
     * Filter a line: maximum or average of 3-pixel neighborhood
     */
    private static
            void Filter3(float[] pixels, int length, int pixel0, int inc)
    {
        float v3 = pixels[pixel0];  //will be pixel[i+1]
        float v2 = v3;              //will be pixel[i]
        float v1;                   //will be pixel[i-1]
        for (int i = 0, p = pixel0; i < length; i++, p += inc)
        {
            v1 = v2;
            v2 = v3;
            if (i < length - 1)
            {
                v3 = pixels[p + inc];
            }
            pixels[p] = (v1 + v2 + v3) * 0.33333333f;
        }
    }

    //  C L A S S   R O L L I N G B A L L
// Copied from ImageJ
    /**
     * A rolling ball (or actually a square part thereof) Here it is also
     * determined whether to shrink the image
     */
    static
            class RollingBall
    {

        float[] data;
        int width;
        int shrinkFactor;

        RollingBall(double radius)
        {
            int arcTrimPer;
            if (radius <= 10)
            {
                shrinkFactor = 1;
                arcTrimPer = 24; // trim 24% in x and y
            }
            else if (radius <= 30)
            {
                shrinkFactor = 2;
                arcTrimPer = 24; // trim 24% in x and y
            }
            else if (radius <= 100)
            {
                shrinkFactor = 4;
                arcTrimPer = 32; // trim 32% in x and y
            }
            else
            {
                shrinkFactor = 8;
                arcTrimPer = 40; // trim 40% in x and y
            }

            double rsquare;     // rolling ball radius squared
            int xtrim;          // # of pixels trimmed off each end of ball to make patch
            int xval, yval;     // x,y-values on patch relative to center of rolling ball
            double smallballradius; // radius of rolling ball (downscaled in x,y and z when image is shrunk)
            int halfWidth;      // distance in x or y from center of patch to any edge (patch "radius")

            smallballradius = radius / shrinkFactor;
            if (smallballradius < 1)
            {
                smallballradius = 1;
            }
            rsquare = smallballradius * smallballradius;
            xtrim = (int) (arcTrimPer * smallballradius) / 100; // only use a patch of the rolling ball
            halfWidth = (int) Math.round(smallballradius - xtrim);
            width = 2 * halfWidth + 1;
            data = new float[width * width];

            for (int y = 0, p = 0; y < width; y++)
            {
                for (int x = 0; x < width; x++, p++)
                {
                    xval = x - halfWidth;
                    yval = y - halfWidth;
                    double temp = rsquare - xval * xval - yval * yval;
                    data[p] = temp > 0. ? (float) (Math.sqrt(temp)) : 0f;
                }
            }
        }
    }
}
