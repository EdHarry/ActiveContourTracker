/*
Copies of the Icy plugins "LabelExtractor" and "ActiveContours" by Alexandre Dufour.
 */
package activecontourtracker;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

/**
 *
 * @author edwardharry
 */
public
        class ROIExtractor
{

    public static
            class ROI
    {

        Vector2i origin;
        Dims dims;
        int[] points;
        int t;
        int trackID;
        int previousTrackID;
        double perimeter;
        double area;
        Vector2 centre;

        ROI(Vector2i origin, Dims dims, int[] points)
        {
            this.origin = origin;
            this.dims = dims;
            this.points = points;
            t = -1;
            trackID = -1;
            previousTrackID = -1;
            perimeter = -1;
            area = -1;
            centre = Vector2.One().Multiply(-1);
        }

        boolean[] MakeMask()
        {
            boolean[] mask = new boolean[dims.width * dims.height];

            for (int index : points)
            {
                mask[index] = true;
            }

            return mask;
        }
    }

    private static
            class ROIBuilder
    {

        Vector2i origin;
        Dims dims;
        LinkedList<Integer> pointXs;
        LinkedList<Integer> pointYs;
        int nPoints;

        public
                ROIBuilder()
        {
            origin = Vector2i.Zero();
            dims = Dims.Zero();

            nPoints = 0;
            pointXs = new LinkedList<>();
            pointYs = new LinkedList<>();
        }

        void AddPoint(int newX, int newY)
        {
            if ((dims.width | dims.height) == 0)
            {
                dims.width = dims.height = 1;
                origin.x = newX;
                origin.y = newY;
            }
            else
            {
                int tmp = origin.x + dims.width - 1;
                if (newX > tmp)
                {
                    dims.width += (newX - tmp);
                }
                else if (newX < origin.x)
                {
                    dims.width += (origin.x - newX);
                    origin.x = newX;
                }

                tmp = origin.y + dims.height - 1;
                if (newY > tmp)
                {
                    dims.height += (newY - tmp);
                }
                else if (newY < origin.y)
                {
                    dims.height += (origin.y - newY);
                    origin.y = newY;
                }
            }

            pointXs.add(newX);
            pointYs.add(newY);
            ++nPoints;
        }

        ROI ToROI()
        {
            Iterator<Integer> itX = pointXs.descendingIterator();
            Iterator<Integer> itY = pointYs.descendingIterator();
            int[] points = new int[nPoints];
            int index = -1;

            while (itX.hasNext())
            {
                points[++index] = ((itY.next() - origin.y) * dims.width)
                                  + itX.next() - origin.x;
            }

            return new ROI(origin, dims, points);
        }
    }

    private static
            class ConnectedComponent
    {

        int targetLabel;
        ConnectedComponent targetComponent;

        ConnectedComponent(int label)
        {
            targetLabel = label;
        }

        int getTargetLabel()
        {
            return targetComponent == null ? targetLabel
                   : targetComponent.getTargetLabel();
        }
    }

    static
            ROI[] SingleImageROIs(byte[] pixels, Dims dims)
    {
        final
                HashMap<Integer, ConnectedComponent> ccs = new HashMap<>();

        final
                HashMap<Integer, ROIBuilder> roiMap = new HashMap<>();

        int[] neighbourLabels = new int[4];
        int nbNeighbours = 0;

        // temporary label buffer for the current slice
        int[] _labelsHere = new int[pixels.length];

        // first image pass: naive labeling with simple backward neighborhood
        int highestKnownLabel = 0;

        for (int y = 0, index = 0;
             y < dims.height;
             ++y)
        {
            for (int x = 0;
                 x < dims.width;
                 ++x, ++index)
            {
                if (pixels[index] == (byte) 0)
                {
                    continue;
                }

                // from here on, the current pixel should be labeled
                // -> look for existing labels in its neighborhood
                // 1) define the neighborhood of interest here
                // NB: this is a single pass method, so backward neighborhood is sufficient
                // legend:
                // e = edge
                // x = current pixel
                // n = valid neighbor
                // . = other neighbor
                if (y == 0)
                {
                    if (x == 0)
                    {
                        // e e e
                        // e x .
                        // e . .

                        // do nothing
                    }
                    else
                    {
                        // e e e
                        // n x .
                        // . . .

                        neighbourLabels[0] = _labelsHere[index - 1];
                        nbNeighbours = 1;
                    }
                }
                else
                {
                    int north = index - dims.width;

                    if (x == 0)
                    {
                        // e n n
                        // e x .
                        // e . .

                        neighbourLabels[0] = _labelsHere[north];
                        neighbourLabels[1] = _labelsHere[north + 1];
                        nbNeighbours = 2;
                    }
                    else if (x == dims.width - 1)
                    {
                        // n n e
                        // n x e
                        // . . e

                        neighbourLabels[0] = _labelsHere[north - 1];
                        neighbourLabels[1] = _labelsHere[north];
                        neighbourLabels[2] = _labelsHere[index - 1];
                        nbNeighbours = 3;
                    }
                    else
                    {
                        // n n n
                        // n x .
                        // . . .

                        neighbourLabels[0] = _labelsHere[north - 1];
                        neighbourLabels[1] = _labelsHere[north];
                        neighbourLabels[2] = _labelsHere[north + 1];
                        neighbourLabels[3] = _labelsHere[index - 1];
                        nbNeighbours = 4;
                    }
                }

                // 2) the neighborhood is ready, move to the labeling step
                // to avoid creating too many labels and fuse them later on,
                // find the minimum non-zero label in the neighborhood
                // and assign that minimum label right now
                int currentLabel = Integer.MAX_VALUE;

                for (int iNeighbour = 0;
                     iNeighbour < nbNeighbours;
                     ++iNeighbour)
                {
                    int neighborLabel = neighbourLabels[iNeighbour];

                    // "zero" neighbors belong to the background...
                    if (neighborLabel == 0)
                    {
                        continue;
                    }

                    // here, the neighbor label is valid
                    // => check if it is lower
                    if (neighborLabel < currentLabel)
                    {
                        currentLabel = neighborLabel;
                    }
                }

                if (currentLabel == Integer.MAX_VALUE)
                {
                    // currentLabel didn't change
                    // => there is no lower neighbor
                    // => register a new connected component
                    ++highestKnownLabel;
                    currentLabel = highestKnownLabel;
                    ccs.put(currentLabel, new ConnectedComponent(currentLabel));
                }
                else
                {
                    // currentLabel has been modified
                    // -> browse the neighborhood again
                    // --> fuse high labels with low labels

                    ConnectedComponent currentCC = ccs.get(currentLabel);
                    int currentTargetLabel = currentCC.getTargetLabel();

                    for (int iNeighbour = 0;
                         iNeighbour < nbNeighbours;
                         ++iNeighbour)
                    {
                        int neighborLabel = neighbourLabels[iNeighbour];

                        if (neighborLabel == 0)
                        {
                            continue; // no object in this pixel
                        }
                        ConnectedComponent neighborCC = ccs.get(neighborLabel);
                        int neighborTargetLabel = neighborCC.getTargetLabel();

                        if (neighborTargetLabel == currentTargetLabel)
                        {
                            continue;
                        }

                        // fuse the highest with the lowest
                        if (neighborTargetLabel > currentTargetLabel)
                        {
                            ccs.get(neighborTargetLabel).targetComponent
                            = ccs.get(currentTargetLabel);
                        }
                        else
                        {
                            ccs.get(currentTargetLabel).targetComponent
                            = ccs.get(neighborTargetLabel);
                        }
                    }
                }

                // -> store this label in the labeled image
                if (currentLabel != 0)
                {
                    _labelsHere[index] = currentLabel;
                }
            }
        }

        // fusion strategy: fuse higher labels with lower ones
        // "highestKnownLabel" holds the highest known label
        // -> loop downwards from there to accumulate object size recursively
        int finalLabel = 0;

        for (int currentLabel = highestKnownLabel;
             currentLabel > 0;
             --currentLabel)
        {
            ConnectedComponent currentCC = ccs.get(currentLabel);

            // if the target label is higher than or equal to the current label
            if (currentCC.targetLabel >= currentLabel)
            {
                // label has same labelValue and targetLabelValue
                // -> it cannot be fused to anything
                // -> this is a terminal label

                // -> assign its final labelValue (for the final image labeling pass)
                ++finalLabel;
                currentCC.targetLabel = finalLabel;
            }
            else
            {
                // the current label should be fused to targetLabel
                currentCC.targetComponent = ccs.get(currentCC.targetLabel);
            }
        }

        // 3) second image pass: replace all labels by their final values
        for (int y = 0, index = 0;
             y < dims.height;
             ++y)
        {
            for (int x = 0;
                 x < dims.width;
                 ++x, ++index)
            {
                int targetLabel = _labelsHere[index];

                if (targetLabel == 0)
                {
                    continue;
                }

                // if a fusion was indicated, retrieve the final label value
                targetLabel = ccs.get(targetLabel).getTargetLabel();

                // store the current pixel in the component
                if (!roiMap.containsKey(targetLabel))
                {
                    roiMap.put(targetLabel, new ROIBuilder());
                }

                roiMap.get(targetLabel).AddPoint(x, y);

            }
        }

        ROI[] rois = new ROI[roiMap.size()];
        int index = -1;

        for (ROIBuilder ROIB : roiMap.values())
        {
            rois[++index] = ROIB.ToROI();
        }

        return rois;
    }
}
