package activecontourtracker;

import activecontourtracker.ContourTracker.Contour.ContourException;
import activecontourtracker.Tracking.Detection;
import activecontourtracker.Tracking.TrackGroup;
import activecontourtracker.Tracking.TrackSegment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.RejectedExecutionException;
import java.awt.geom.Line2D;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;

class ContourTracker
{

    interface ContourUpdate
    {

        void Refresh();

        void AddContour(Contour c);

        void UpdateContour(Contour c);

        void RemoveContour(Contour c);
    }

    static
            class ContourMultiUpdate implements ContourUpdate
    {

        private final
                ContourUpdate[] updaters;

        public
                ContourMultiUpdate(ContourUpdate[] updaters)
        {
            this.updaters = updaters;
        }

        @Override
        public
                void Refresh()
        {
            for (ContourUpdate contourUpdate : updaters)
            {
                contourUpdate.Refresh();
            }
        }

        @Override
        public
                void AddContour(Contour c)
        {
            for (ContourUpdate contourUpdate : updaters)
            {
                contourUpdate.AddContour(c);
            }
        }

        @Override
        public
                void UpdateContour(Contour c)
        {
            for (ContourUpdate contourUpdate : updaters)
            {
                contourUpdate.UpdateContour(c);
            }
        }

        @Override
        public
                void RemoveContour(Contour c)
        {
            for (ContourUpdate contourUpdate : updaters)
            {
                contourUpdate.RemoveContour(c);
            }
        }
    }

    static
            class Contour implements Cloneable, Iterable<Vector2>
    {

        private static final
                class Segment implements Iterable<Vector2>
        {

            final
                    ArrayList<Vector2> points;

            Segment(Vector2 head, Vector2 tail)
            {
                points = new ArrayList<>(2);
                points.add(head);
                points.add(tail);
            }

            final
                    Vector2 GetHead()
            {
                return points.get(0);
            }

            final
                    Vector2 GetTail()
            {
                return points.get(points.size() - 1);
            }

            final
                    void AddHead(Vector2 p)
            {
                points.add(0, p);
            }

            final
                    void AddHead(Segment s)
            {
                for (int index = 0;
                     index < s.points.size();
                     ++index)
                {
                    points.add(index, s.points.get(index));
                }
            }

            final
                    void AddTail(Vector2 p)
            {
                points.add(p);
            }

            @Override
            public
                    Iterator<Vector2> iterator()
            {
                return points.iterator();
            }
        }

        // Sliding Window, based on coefficient of variation
        static
                class SlidingWindow
        {

            private final
                    int windowLength;
            private final
                    double covCritSq;
            private final
                    double[] winx;
            private
                    double sumx;
            private
                    double sumx2;
            private
                    int addPtr;
            private
                    int subPtr;
            private
                    boolean init;
            private
                    boolean hasConv;
            private
                    boolean switchFlag;

            public
                    SlidingWindow(int length, double covCritSq)
            {
                windowLength = length;
                this.covCritSq = covCritSq;
                winx = new double[windowLength + 1];
                sumx = 0;
                sumx2 = 0;
                addPtr = windowLength;
                subPtr = 0;
                init = false;
                hasConv = false;
                switchFlag = true;
            }

            public
                    SlidingWindow(int length, double covCritSq,
                                  double[] initVals)
            {
                this(length, covCritSq);
                System.arraycopy(initVals, 0, winx, 0, initVals.length);
                subPtr = (initVals.length + 1) % windowLength;
                addPtr = initVals.length % windowLength;
                for (double d : initVals)
                {
                    sumx += d;
                    sumx2 += (d * d);
                }
            }

            public
                    SlidingWindow(SlidingWindow window)
            {
                windowLength = window.windowLength;
                covCritSq = window.covCritSq;

                winx = new double[windowLength + 1];
                System.arraycopy(window.winx, 0, winx, 0, winx.length);

                sumx = window.sumx;
                sumx2 = window.sumx2;
                addPtr = window.addPtr;
                subPtr = window.subPtr;
                init = window.init;
                hasConv = window.hasConv;
                switchFlag = window.switchFlag;
            }

            void Add(double val)
            {
                switchFlag = !switchFlag;
                if (switchFlag)
                {
                    return;
                }

                double subVal = winx[subPtr++];
                sumx -= subVal;
                sumx += val;
                winx[addPtr++] = val;

                sumx2 += (val * val);
                sumx2 -= (subVal * subVal);

                if (addPtr > windowLength)
                {
                    addPtr = 0;
                }
                else if (subPtr > windowLength)
                {
                    init = true;
                    subPtr = 0;
                }

                hasConv = init && ((sumx2 * windowLength / (sumx
                                                            * sumx)) - 1) < covCritSq;
            }

            boolean HasConverged()
            {
                return hasConv;
            }

            int GetLength()
            {
                return windowLength;
            }

            SlidingWindow Copy()
            {
                return new SlidingWindow(this);
            }

            SlidingWindow Expand(double factor)
            {
                return new SlidingWindow((int) (windowLength * factor),
                                         covCritSq, winx);
            }
        }

        /**
         * Copied from Icy plugins.adufour.morphology
         *
         * Fill holes in the specified boolean mask. The algorithm extracts the
         * background by flooding from the mask edges, then fills all remaining
         * background pixels (i.e. the holes in the mask)
         *
         * @param mask
         * @return <code>true</code> if the mask has changed, <code>false</code>
         * otherwise (i.e. no hole was detected)
         */
        private

                void FillHoles(ROIExtractor.ROI roi)
        {
            int width = roi.dims.width;
            int height = roi.dims.height;

            // avoid single line / column masks
            if (width == 1 || height == 1)
            {
                return;
            }

            boolean[] isObject = roi.MakeMask();
            int slice = isObject.length;

            byte[] labels = new byte[slice];
            byte UNKNOWN = 0;
            byte TO_VISIT = 1;
            byte BACKGROUND = 2;

            // hold a list of pixels to visit
            int[] pixelOffsetsToVisit = new int[slice];

            int n = 0;

            // 1) extract background pixels on the image edges
            // 1.a) top + bottom edges
            for (int top = 0, bottom = slice - width;
                 top < width;
                 ++top, ++bottom)
            {
                if (!isObject[top])
                {
                    labels[top] = BACKGROUND;

                    int neighbor = top + width;
                    if (labels[neighbor] == UNKNOWN)
                    {
                        pixelOffsetsToVisit[n++] = neighbor;
                        labels[neighbor] = TO_VISIT;
                    }
                }

                if (!isObject[bottom])
                {
                    labels[bottom] = BACKGROUND;

                    int neighbor = bottom - width;
                    if (labels[neighbor] == UNKNOWN)
                    {
                        pixelOffsetsToVisit[n++] = neighbor;
                        labels[neighbor] = TO_VISIT;
                    }
                }
            }

            // 1.b) left + right edges (minus the corners, done above)
            for (int left = width, right = 2 * width - 1;
                 left < slice - width;
                 left += width, right += width)
            {
                if (!isObject[left])
                {
                    labels[left] = BACKGROUND;

                    int neighbor = left + 1;
                    if (labels[neighbor] == UNKNOWN)
                    {
                        pixelOffsetsToVisit[n++] = neighbor;
                        labels[neighbor] = TO_VISIT;
                    }
                }

                if (!isObject[right])
                {
                    labels[right] = BACKGROUND;

                    int neighbor = right - 1;
                    if (labels[neighbor] == UNKNOWN)
                    {
                        pixelOffsetsToVisit[n++] = neighbor;
                        labels[neighbor] = TO_VISIT;
                    }
                }
            }

            // process only if there is something to process...
            if (n > 0)
            {
                // 2) flood the image from the list of (background) border pixels

                int[] neighbors = new int[4];

                for (int index = 0;
                     index < n;
                     ++index)
                {
                    int offset = pixelOffsetsToVisit[index];

                    // if the current pixel is already tagged 'background'
                    // it's neighbors are already in the list of pixels to visit
                    // => don't go further
                    if (labels[offset] == BACKGROUND)
                    {
                        continue;
                    }

                    // if background => flood
                    if (!isObject[offset])
                    {
                        // tag the current pixel
                        labels[offset] = BACKGROUND;

                        // calculate neighbor offsets
                        neighbors[0] = offset - 1;
                        neighbors[1] = offset + 1;
                        neighbors[2] = offset - width;
                        neighbors[3] = offset + width;

                        // and add its immediate neighbor to the list of pixels to visit
                        for (int neighbor : neighbors)
                        {
                            if (labels[neighbor] == UNKNOWN && !isObject[neighbor])
                            {
                                pixelOffsetsToVisit[n++] = neighbor;
                                labels[neighbor] = TO_VISIT;
                            }
                        }
                    }
                }
            }

            // 3) fill holes by filling all "UNKNOWN" pixels out of the mask
            LinkedList<Integer> roiPoints = new LinkedList<>();
            int nNew = 0;

            for (int index = 0;
                 index < slice;
                 ++index)
            {
                if (labels[index] == UNKNOWN && !isObject[index])
                {
                    roiPoints.add(index);
                    ++nNew;
                }
            }

            if (nNew > 0)
            {
                int index = roi.points.length;
                int[] newPoints = new int[index + nNew];
                System.arraycopy(roi.points, 0, newPoints, 0, index);
                for (int point : roiPoints)
                {
                    newPoints[index++] = point;
                }
                roi.points = newPoints;
            }
        }

        final
                ArrayList<Vector2> points = new ArrayList<>();

        private
                Vector2[] modelForces;

        Vector2[] contourNormals;

        private
                Vector2[] feedbackForces;

        private
                Vector2[] volumeConstraintForces;

        private
                double sampling;

        private
                SlidingWindow convergenceWindow;

        private
                BoundingCircle boundingCircle;

        private
                double divisionSensitivity;

        private
                boolean clockWise;

        /**
         * Creates a clone of the specified contour
         *
         * @param contour
         */
        public
                Contour(Contour contour)
        {
            sampling = contour.sampling;
            divisionSensitivity = contour.divisionSensitivity;
            convergenceWindow = contour.convergenceWindow
                    .Copy();

            int n = contour.points.size();

            points.ensureCapacity(n);
            contourNormals = new Vector2[n];
            modelForces = new Vector2[n];
            feedbackForces = new Vector2[n];
            volumeConstraintForces = new Vector2[n];

            for (int index = 0;
                 index < n;
                 ++index)
            {
                contourNormals[index] = Vector2.Zero();
                modelForces[index] = Vector2.Zero();
                feedbackForces[index] = Vector2.Zero();
                volumeConstraintForces[index] = Vector2.Zero();
                AddPoint(contour.points.get(index).Copy());
            }

            UpdateNormals();
            MakeBoundingCircle();
        }

        /**
         * Creates a clone of the specified contour
         *
         * @param contour
         */
        public
                Contour(Contour contour, double windowExpandFactor)
        {
            sampling = contour.sampling;
            divisionSensitivity = contour.divisionSensitivity;
            convergenceWindow = contour.convergenceWindow
                    .Expand(windowExpandFactor);

            int n = contour.points.size();

            points.ensureCapacity(n);
            contourNormals = new Vector2[n];
            modelForces = new Vector2[n];
            feedbackForces = new Vector2[n];
            volumeConstraintForces = new Vector2[n];

            for (int index = 0;
                 index < n;
                 ++index)
            {
                contourNormals[index] = Vector2.Zero();
                modelForces[index] = Vector2.Zero();
                feedbackForces[index] = Vector2.Zero();
                volumeConstraintForces[index] = Vector2.Zero();
                AddPoint(contour.points.get(index).Copy());
            }

            UpdateNormals();
            MakeBoundingCircle();
        }

        public
                Contour(double sampling, double divisionSensitivity,
                        int slidingWindowSize,
                        double convCritSq, ROIExtractor.ROI roi)
        {
            this.sampling = sampling;
            this.divisionSensitivity = divisionSensitivity;
            convergenceWindow = new SlidingWindow(slidingWindowSize, convCritSq);

            FillHoles(roi);

            try
            {
                Triangulate(roi);
                MakeBoundingCircle();
                ReSample(0.8, 1.4);
            }
            catch (ContourException e)
            {
                points.clear();

                double halfWidth = roi.dims.width / 2;
                double halfHeight = roi.dims.height / 2;
                double currentAngle = 0;

                Vector2 offset = new Vector2(halfWidth, halfHeight).
                        Add(roi.origin.ToVector2());

//                while (currentAngle < 6.283185307179586) // 2PI
                while (currentAngle < 6.282) // 2PI - 0.001...
                {
                    double factor1 = halfHeight * Math.cos(currentAngle);
                    double factor2 = halfWidth * Math.sin(currentAngle);
                    double factor3 = Math.sqrt((factor1 * factor1)
                                               + (factor2 * factor2));
                    double newX = factor1 * halfWidth / factor3;
                    double newY = factor2 * halfHeight / factor3;

                    AddPoint(new Vector2(newX, newY).Add(offset));

                    //currentAngle += 0.157079632679490; // PI/20
                    currentAngle += 1.57079632679490; // PI/2
                }

                contourNormals = null;

                try
                {
                    ReSample(0.8, 1.4);
                }
                catch (ContourException ee)
                {
                    throw new RuntimeException("Couldn't create a contour from ROI: "
                                               + roi.toString() + " (located at " + offset.x + " "
                                               + "; " + offset.y + ")");
                }
            }

            MakeBoundingCircle();
        }

        public
                Contour(double sampling, double divisionSensitivity,
                        SlidingWindow window, LinkedList<Vector2> points)
        {
            this.sampling = sampling;
            this.divisionSensitivity = divisionSensitivity;
            convergenceWindow = window;

            points.stream().forEach((point)
                    ->
            {
                AddPoint(point);
            });

            int n = points.size();

            contourNormals = new Vector2[n];
            modelForces = new Vector2[n];
            feedbackForces = new Vector2[n];
            volumeConstraintForces = new Vector2[n];

            for (int index = 0;
                 index < n;
                 ++index)
            {
                contourNormals[index] = Vector2.Zero();
                modelForces[index] = Vector2.Zero();
                feedbackForces[index] = Vector2.Zero();
                volumeConstraintForces[index] = Vector2.Zero();
            }

            UpdateNormals();
            MakeBoundingCircle();
        }

        private

                void AddPoint(Vector2 p)
        {
            points.add(p);
        }

        private
                class BoundingCircle
        {

            private
                    double radiusSq;
            private
                    Vector2 centre;

            public
                    BoundingCircle()
            {
                centre = Vector2.Zero();
                points.stream().forEach((point)
                        ->
                {
                    centre.AddTo(point);
                });
                centre.DivideBy(points.size());

                radiusSq = 0;
                points.stream().forEach((point)
                        ->
                {
                    double disSq = point.DistanceSq(centre);
                    if (disSq > radiusSq)
                    {
                        radiusSq = disSq;
                    }
                });
            }

            double GetRadius()
            {
                return Math.sqrt(radiusSq);
            }

            double GetRadiusSq()
            {
                return radiusSq;
            }

            Vector2 GetCentre()
            {
                return centre;
            }

            Vector2 GetCentre_OnlyPositive(Dims dims)
            {
                Vector2 c = centre.Copy();

                c.x = c.x < 0 ? 0 : c.x;
                c.x = c.x >= dims.width ? dims.width : c.x;
                c.y = c.y < 0 ? 0 : c.y;
                c.y = c.y >= dims.height ? dims.height : c.y;

                return c;
            }

            boolean DoesIntersect(BoundingCircle other)
            {
                return centre.DistanceSq(other.centre) < (radiusSq + other.radiusSq
                                                          + (2 * GetRadius() * other.GetRadius()));
            }
        }

        private
                void MakeBoundingCircle()
        {
            boundingCircle = new BoundingCircle();
        }

        BoundingCircle GetBoundingCircle()
        {
            return boundingCircle;
        }

        /**
         * Checks whether the contour is self-intersecting. Depending on the
         * given parameters, a self-intersection can be considered as a loop or
         * as a contour division.
         *
         * @param minDistance the distance threshold between non-neighboring
         * points to detect self-intersection
         * @return <code>null</code> if either no self-intersection is detected
         * or if one of the new contours is too small, or an array of contours
         * with 0 elements if both contours are too small, and 2 elements if
         * both contours are viable
         */
        protected
                Contour[] CheckSelfIntersection(double minDistance)
        {
            int index_0, index_1 = 0, n = points.size();
            Vector2 p_i = Vector2.Zero(), p_j = Vector2.Zero();

            double divisionDistQ = boundingCircle.GetRadiusSq()
                                   * 4 * divisionSensitivity
                                   * divisionSensitivity;

            double minDistanceQ = minDistance;
            minDistanceQ *= minDistanceQ;

            boolean selfIntersection = false;

            loop:
            for (index_0 = 0;
                 index_0 < n;
                 ++index_0)
            {
                p_i = points.get(index_0);
                Vector2 n_i = contourNormals[index_0];

                for (index_1 = index_0 + 2;
                     index_1 < n - 1;
                     ++index_1)
                {
                    p_j = points.get(index_1);

                    double distQ = p_i.DistanceSq(p_j);

                    if (distQ < minDistanceQ)
                    {
                        // local self-intersection
                        // deal with the special case that i and j are 2 points away

                        if (index_0 == 0 && index_1 == n - 2)
                        {
                            points.remove(--n);
                            continue;
                        }
                        else if (index_0 == 1 && index_1 == n - 1)
                        {
                            points.remove(0);
                            --n;
                            continue;
                        }
                        else if (index_1 == index_0 + 2)
                        {
                            points.remove(index_0 + 1);
                            --n;
                            continue;
                        }

                        // a real self-intersection is happening
                        selfIntersection = true;
                        break loop;
                    }

                    // self-intersection always involves opposite normals
                    if (n_i.Dot(contourNormals[index_1]) > -0.5)
                    {
                        continue;
                    }

                    // look for division
                    // what about the intensity profile between v1 and v2?
                    // ROI2DLine line = new ROI2DLine(points.get(i).x, points.get(i).y,
                    // points.get(j).x, points.get(j).y);
                    // int[] linePoints = line.getBooleanMask(true).getPointsAsIntArray();
                    // are points sufficiently close?
                    // => use the bounding radius / 2
                    if (distQ < divisionDistQ)
                    {
                        // are points located "in front" of each other?
                        // j - i ~= n/2 ?
                        // take a loose guess (+/- n/7)
                        if ((index_1 - index_0) > (2 * n / 5) && (index_1
                                                                  - index_0) < (3 * n / 5))
                        {
                            // check the local curvature on each side (4 points away)

                            Vector2 vi1 = points.get((index_0 + n - 4) % n).Copy();
                            Vector2 vi2 = points.get((index_0 + 4) % n).Copy();
                            vi1.SubtractFrom(p_i);
                            vi2.SubtractFrom(p_i);
                            vi1.Normalise();
                            vi2.Normalise();
                            double z = vi1.Cross(vi2);
                            // discard small of positive curvatures (i.e. z < 0)
                            if (z < 0.05)
                            {
                                continue;
                            }

                            double vi_lQ = vi1.LengthSq();
                            // System.out.println(vi1.lengthSquared());
                            Vector2 vj1 = points.get((index_1 + n - 4) % n).Copy();
                            Vector2 vj2 = points.get((index_1 + 4) % n).Copy();
                            vj1.SubtractFrom(p_j);
                            vj2.SubtractFrom(p_j);
                            vj1.Normalise();
                            vj2.Normalise();
                            z = vj1.Cross(vj2);

                            // discard small of positive curvatures (i.e. z < 0)
                            if (z < 0.05)
                            {
                                continue;
                            }

                            double vj_lQ = vj1.LengthSq();
                            // System.out.println(vj1.lengthSquared());

                            // curvature has to be at a relative minimum
                            if (vi_lQ < 0.5 && vj_lQ < 0.5)
                            {
                                continue;
                            }

                            // System.out.println(vi_lQ + " | " + vj_lQ);
                            // a real self-intersection is happening
                            selfIntersection = true;
                            break loop;
                        }
                    }

                }
            }

            if (!selfIntersection)
            {
                return null;
            }

            int nPoints = index_1 - index_0;
            LinkedList<Vector2> childPoints = new LinkedList<>();
            for (int index = 0;
                 index < nPoints;
                 ++index)
            {
                childPoints.add(points.get(index + index_0));
            }

            Contour child1 = new Contour(sampling, divisionSensitivity,
                                         convergenceWindow.Copy(), childPoints);

            childPoints.clear();
            nPoints = index_0 + n - index_1;
            for (int index = 0, pj = index + index_1;
                 index < nPoints;
                 ++index, ++pj)
            {
                childPoints.add(points.get(pj % n));
            }
            Contour child2 = new Contour(sampling, divisionSensitivity,
                                         convergenceWindow.Copy(), childPoints);

            // determine whether the intersection is a loop or a division
            // rationale: check the normal of the two colliding points (i & j)
            // if they point away from the junction => division
            // if they point towards the junction => loop
            Vector2 n_i = contourNormals[index_0];
            Vector2 i_j = p_j.Subtract(p_i);

            if (n_i.Dot(i_j) < 0)
            {
                // division => keep c1 and c2 if their size is ok

                double c1area = Math.abs(child1.GetAlgebraicInterior());
                double c2area = Math.abs(child2.GetAlgebraicInterior());

                // if only one of the two children has a size lower than minArea, then the division
                // should be considered as an artifact loop, the other child thus is the new contour
//                if (child1.points.size() < 10 || c1area < c2area / 5)
                if (c1area < c2area / 5)
                {
                    // remove c1 (too small)
                    points.clear();
                    points.addAll(child2.points);
                    return null;
                }

//                if (child2.points.size() < 10 || c2area < c1area / 5)
                if (c2area < c1area / 5)
                {
                    // remove c2 (too small)
                    points.clear();
                    points.addAll(child1.points);
                    return null;
                }

                // keep both then...
                return new Contour[]
                {
                    child1, child2
                };
            }

//            // loop => keep only the contour with correct orientation
//            // => the contour with a positive algebraic area
//            if (child1.counterClockWise == counterClockWise)
            if (child1.points.size() > child2.points.size())
            {
                // c1 is the outer loop => keep it
                points.clear();
                points.addAll(child1.points);
                return null;
            }

            // c1 is the inner loop => keep c2
            points.clear();
            points.addAll(child2.points);
            return null;
        }

        @Override
        public
                Contour clone() throws CloneNotSupportedException
        {
            super.clone();
            return new Contour(this);
        }

        /**
         * Update the axis constraint force, which adjusts the takes the final
         * forces and normalize them to keep the contour shape along its
         * principal axis <br>
         * WARNING: this method directly update the array of final forces used
         * to displace the contour points. It should be used among the last to
         * keep it most effective
         *
         * @param weight
         */
        void ComputeAxisForces(double weight)
        {
            Vector2 axis = Vector2.Unit();
            int s = points.size();

            // Compute the object axis as the vector between the two most distant
            // contour points
            // TODO this is not optimal, geometric moments should be used
            {
                double maxDistSq = 0;
                Vector2 vec;

                for (int index_0 = 0;
                     index_0 < s;
                     ++index_0)
                {
                    Vector2 vi = points.get(index_0);

                    for (int index_1 = index_0 + 1;
                         index_1 < s;
                         ++index_1)
                    {
                        Vector2 vj = points.get(index_1);

                        vec = vi.Subtract(vj);
                        double dSq = vec.LengthSq();

                        if (dSq > maxDistSq)
                        {
                            maxDistSq = dSq;
                            axis = vec.Copy();
                        }
                    }
                }
                axis.Normalise();
            }

            // To drive the contour along the main object axis, each displacement
            // vector is scaled by the scalar product between its normal and the main axis.
            {
                for (int index = 0;
                     index < s;
                     ++index)
                {
                    Vector2 normal = contourNormals[index];

                    // dot product between normalized vectors ranges from -1 to 1
                    double colinearity = Math.abs(normal.Dot(axis)); // now from 0 to 1

                    // goal: adjust the minimum using the weight, but keep max to 1
                    double threshold = Math.max(colinearity, 1 - weight);

                    modelForces[index].MultiplyBy(threshold);
                }
            }
        }

        void ComputeBalloonForces(double weight)
        {
            int n = points.size();

            for (int index = 0;
                 index < n;
                 ++index)
            {
                modelForces[index].AddTo(contourNormals[index].Multiply(weight));
            }
        }

        /**
         * Update edge term of the contour evolution according to the image
         * gradient
         *
         * @param weight
         * @param edge_data
         */
        void ComputeEdgeForces(byte[] edges, Dims dims, double weight)
        {
            int n = points.size();

            for (int index = 0;
                 index < n;
                 ++index)
            {
                Vector2 p = points.get(index);
                Vector2 force = modelForces[index];

                // compute the gradient (2nd order)
                double nextX = GetPixelValue(edges, dims, p.Add(new Vector2(0.5, 0)));
//                if (nextX == 0)
//                {
//                    continue;
//                }
                double prevX = GetPixelValue(edges, dims, p.Subtract(new Vector2(0.5, 0)));
//                if (prevX == 0)
//                {
//                    continue;
//                }
                double nextY = GetPixelValue(edges, dims, p.Add(new Vector2(0, 0.5)));
//                if (nextY == 0)
//                {
//                    continue;
//                }
                double prevY = GetPixelValue(edges, dims, p.Subtract(new Vector2(0, 0.5)));
//                if (prevY == 0)
//                {
//                    continue;
//                }

                Vector2 grad = new Vector2(nextX - prevX, nextY - prevY);
                grad.MultiplyBy(weight);

                force.AddTo(grad);
            }
        }

        void ComputeRegionForces(byte[] image, Dims dims, double weight,
                                 double sensitivity)
        {
            // sensitivity should be high for dim objects, low for bright objects...
            // ... but none of the following options work properly
            // sensitivity *= 1/(1+cin);
            // sensitivity = sensitivity * cin / (0.01 + cout);
            // sensitivity = sensitivity / (2 * Math.max(cout, cin));
            // sensitivity = sensitivity / (Math.log10(cin / cout));

            Vector2 p, force, norm, cvms;
            double val, inDiff, outDiff, forceFactor;
            int n = points.size();

            weight *= sampling;

            for (int index = 0;
                 index < n;
                 ++index)
            {
                p = points.get(index);
                force = modelForces[index];
                norm = contourNormals[index];

                // bounds check
                // if (p.x <= 1 || p.y <= 1 || p.x >= width - 2 || p.y >= height - 2) continue;
                val = GetPixelValue(image, dims, p);

                inDiff = val - 255;
                inDiff *= inDiff;

                outDiff = val;
                outDiff *= outDiff;

                forceFactor = weight * ((sensitivity * outDiff) - (inDiff / sensitivity));

                cvms = norm.Multiply(forceFactor);

                force.AddTo(cvms);
            }
        }

        void ComputeInternalForces(double weight)
        {
            if (feedbackForces == null)
            {
                return;
            }

            int n = points.size();

            if (n < 3)
            {
                return;
            }

            Vector2 force;
            Vector2 prev, curr, next;

            weight /= sampling;

            // first point
            prev = points.get(n - 1);
            curr = points.get(0);
            next = points.get(1);

            force = feedbackForces[0];

            force.AddTo(prev.Add(next).Subtract(curr.Multiply(2)).
                    Multiply(weight));

            // middle points
            for (int index = 1;
                 index < n - 1;
                 ++index)
            {
                force = feedbackForces[index];
                prev = points.get(index - 1);
                curr = points.get(index);
                next = points.get(index + 1);

                force.AddTo(prev.Add(next).Subtract(curr.Multiply(2)).
                        Multiply(weight));
            }

            // last point
            force = feedbackForces[n - 1];
            prev = points.get(n - 2);
            curr = points.get(n - 1);
            next = points.get(0);

            force.AddTo(prev.Add(next).Subtract(curr.Multiply(2)).
                    Multiply(weight));
        }

        void ComputeVolumeConstraint(double targetVolume)
        {
            // 1) compute the difference between target and current volume
            double volumeDiff = targetVolume - Math.abs(GetAlgebraicInterior());
            // if (volumeDiff > 0): contour too small, should no longer shrink
            // if (volumeDiff < 0): contour too big, should no longer grow

            int n = points.size();

            Vector2 avgFeedback = Vector2.Zero();
            int nbFeedbackForces = 0;

            for (int index = 0;
                 index < n;
                 ++index)
            {
                Vector2 mf = modelForces[index];

                // 2) check whether the final force has same direction as the outer normal
                double forceNorm = mf.Dot(contourNormals[index]);

                // if forces have same direction (forceNorm > 0): contour is growing
                // if forces have opposite direction (forceNorm < 0): contour is shrinking
                // if (forceNorm * volumeDiff < 0)
                // {
                // // forceNorm and volumeDiff have opposite signs because:
                // // - contour too small (volumeDiff > 0) and shrinking (forceNorm < 0)
                // // or
                // // - contour too large (volumeDiff < 0) and growing (forceNorm > 0)
                // // => in both cases, constrain the final force accordingly
                // mf.scale(1.0 / (1.0 + Math.abs(volumeDiff) / 10));
                // }
                // estimate an average feedback
                if (forceNorm > 0 && volumeDiff < 0)
                {
                    avgFeedback.AddTo(feedbackForces[index]);
                    ++nbFeedbackForces;
                }
            }

            if (avgFeedback.LengthSq() > 0)
            {
                avgFeedback.DivideBy(nbFeedbackForces);
                avgFeedback.MultiplyBy(Math.abs(volumeDiff / targetVolume) / 0.5);

                // move the entire mesh (ugly, but amazingly efficient!!)
                for (int index = 0;
                     index < n;
                     ++index)
                {
                    volumeConstraintForces[index].AddTo(avgFeedback);
                }
            }
        }

        /**
         * Computes the feedback forces yielded by the penetration of the
         * current contour into the target contour
         *
         * @param target the contour that is being penetrated
         * @return the number of actual point-mesh intersection tests
         */
        int ComputeFeedbackForces(Contour target)
        {
            Vector2 targetCenter = target.GetBoundingCircle().GetCentre();

            double targetRadiusSq = target.GetBoundingCircle().GetRadiusSq();

            double penetration;

            int tests = 0;
            int index = 0;

            for (Vector2 p : points)
            {
                double distanceSq = p.DistanceSq(targetCenter);

                if (distanceSq < targetRadiusSq)
                {
                    ++tests;

                    if ((penetration = target.GetDistanceToEdge(p)) > 0)
                    {
                        Vector2 feedbackForce = feedbackForces[index];

                        // feedbackForce.scale(-penetration, contourNormals[index]);
                        feedbackForce.AddTo(contourNormals[index].Multiply(
                                -penetration * 0.5));

                        modelForces[index].MultiplyBy(0.05);
                    }
                }

                ++index;
            }

            return tests;
        }

        private static
                void CreateEdge(ArrayList<Segment> segments, double xStart,
                                double yStart, double xEnd, double yEnd)
        {
            double EPSILON = 0.0000000001;

            Vector2 head = new Vector2(xStart, yStart);
            Vector2 tail = new Vector2(xEnd, yEnd);

            if (segments.isEmpty())
            {
                segments.add(new Segment(head, tail));
                return;
            }

            int insertAtTailOf = -1, insertAtHeadOf = -1;

            for (int index = 0;
                 index < segments.size();
                 ++index)
            {
                if (tail.DistanceSq(segments.get(index).
                        GetHead()) <= EPSILON)
                {
                    insertAtHeadOf = index;
                }
                else if (head.DistanceSq(segments.get(index).
                        GetTail()) <= EPSILON)
                {
                    insertAtTailOf = index;
                }
            }

            if (insertAtTailOf >= 0)
            {
                if (insertAtHeadOf >= 0)
                {
                    segments.get(insertAtHeadOf).
                            AddHead(segments.get(insertAtTailOf));
                    segments.remove(insertAtTailOf);
                }
                else
                {
                    segments.get(insertAtTailOf).AddTail(tail);
                }
            }
            else if (insertAtHeadOf >= 0)
            {
                segments.get(insertAtHeadOf).AddHead(head);
            }
            else
            {
                segments.add(new Segment(head, tail));
            }
        }

        /**
         * Calculates the 2D image value at the given real coordinates by
         * bilinear interpolation
         *
         * @param imageFloat the image to sample (must be of type
         * {@link DataType#DOUBLE})
         * @param x the X-coordinate of the point
         * @param y the Y-coordinate of the point
         * @return the interpolated image value at the given coordinates
         */
        private static
                float GetPixelValue(byte[] data, Dims dims, Vector2 point_in)
        {
            // "center" the coordinates to the center of the pixel
            //Vector2 point = point_in.Subtract(0.5);

            Vector2 point = point_in.Copy();

            int i = (int) Math.floor(point.x);
            int j = (int) Math.floor(point.y);

            // if (i <= 0 || i >= width - 1 || j <= 0 || j >= height - 1) return 0f;
            if (i < 0)
            {
                i = 0;
            }
            if (j < 0)
            {
                j = 0;
            }
            if (i > dims.width - 2)
            {
                i = dims.width - 2;
            }
            if (j > dims.height - 2)
            {
                j = dims.height - 2;
            }

            float value = 0;

            final
                    int offset = i + j * dims.width;
            final
                    int offset_plus_1 = offset + 1; // saves 1 addition

            point.SubtractFrom(new Vector2(i, j));

            final
                    Vector2 m = Vector2.One().Subtract(point);

            value += (m.x * m.y * (data[offset] & 0xff));
            value += (point.x * m.y * (data[offset_plus_1] & 0xff));
            value += (m.x * point.y * (data[offset + dims.width] & 0xff));
            value += (point.x * point.y * (data[offset_plus_1 + dims.width] & 0xff));

            return value;
        }

        /**
         * Computes the algebraic area of the current contour. The returned
         * value is negative if the contour points are order clockwise and
         * positive if ordered counter-clockwise. The contour's surface is just
         * the absolute value of this algebraic surface
         *
         * @return
         */
        private
                double GetAlgebraicInterior()
        {
            int nm1 = points.size() - 1;
            double area = 0;

            // all points but the last
            for (int index = 0;
                 index < nm1;
                 ++index)
            {
                area += points.get(index + 1).Cross(points.get(index));
            }

            // last point
            area += points.get(0).Cross(points.get(nm1));

            return area * 0.5;
        }

        private
                double GetAlgebraicInterior_OnlyPositive(Dims dims)
        {
            int nm1 = points.size() - 1;
            double area = 0;

            // all points but the last
            for (int index = 0;
                 index < nm1;
                 ++index)
            {
                Vector2 p0 = points.get(index).Copy();
                Vector2 p1 = points.get(index + 1).Copy();

                p0.x = p0.x < 0 ? 0 : p0.x;
                p0.x = p0.x >= dims.width ? dims.width : p0.x;
                p0.y = p0.y < 0 ? 0 : p0.y;
                p0.y = p0.y >= dims.height ? dims.height : p0.y;

                p1.x = p1.x < 0 ? 0 : p1.x;
                p1.x = p1.x >= dims.width ? dims.width : p1.x;
                p1.y = p1.y < 0 ? 0 : p1.y;
                p1.y = p1.y >= dims.height ? dims.height : p1.y;

                area += p1.Cross(p0);
            }

            // last point
            Vector2 p0 = points.get(nm1).Copy();
            Vector2 p1 = points.get(0).Copy();

            p0.x = p0.x < 0 ? 0 : p0.x;
            p0.x = p0.x >= dims.width ? dims.width : p0.x;
            p0.y = p0.y < 0 ? 0 : p0.y;
            p0.y = p0.y >= dims.height ? dims.height : p0.y;

            p1.x = p1.x < 0 ? 0 : p1.x;
            p1.x = p1.x >= dims.width ? dims.width : p1.x;
            p1.y = p1.y < 0 ? 0 : p1.y;
            p1.y = p1.y >= dims.height ? dims.height : p1.y;

            area += p1.Cross(p0);

            return area * 0.5;
        }

        double GetPerimeter()
        {
            int size = points.size();

            if (size == 0)
            {
                return 0;
            }

            Vector2 p1 = points.get(size - 1);
            Vector2 p2 = points.get(0);

            double perimeter = p1.Distance(p2);

            for (int index = 0;
                 index < size - 1;
                 ++index)
            {
                // shift pair of points by one index
                p1 = p2;
                p2 = points.get(index + 1);
                perimeter += p1.Distance(p2);
            }

            return perimeter;
        }

        double GetPerimeter_OnlyPositive(Dims dims)
        {
            int size = points.size();

            if (size == 0)
            {
                return 0;
            }

            Vector2 p1 = points.get(size - 1).Copy();
            Vector2 p2 = points.get(0).Copy();

            p1.x = p1.x < 0 ? 0 : p1.x;
            p1.x = p1.x >= dims.width ? dims.width : p1.x;
            p1.y = p1.y < 0 ? 0 : p1.y;
            p1.y = p1.y >= dims.height ? dims.height : p1.y;

            p2.x = p2.x < 0 ? 0 : p2.x;
            p2.x = p2.x >= dims.width ? dims.width : p2.x;
            p2.y = p2.y < 0 ? 0 : p2.y;
            p2.y = p2.y >= dims.height ? dims.height : p2.y;

            double perimeter = p1.Distance(p2);

            for (int index = 0;
                 index < size - 1;
                 ++index)
            {
                // shift pair of points by one index
                p1 = p2;
                p2 = points.get(index + 1).Copy();

                p2.x = p2.x < 0 ? 0 : p2.x;
                p2.x = p2.x >= dims.width ? dims.width : p2.x;
                p2.y = p2.y < 0 ? 0 : p2.y;
                p2.y = p2.y >= dims.height ? dims.height : p2.y;

                perimeter += p1.Distance(p2);
            }

            return perimeter;
        }

        /**
         * Tests whether the given point is inside the contour, and if so
         * returns the penetration depth of this point. <br>
         * This methods computes the number of intersections between the contour
         * and a semi-infinite line starting from the contour center and passing
         * through the given point. The point is thus considered inside if the
         * number of intersections is odd (Jordan curve theorem).<br>
         * Implementation note: the AWT Line2D class only provides a "segment to
         * segment" intersection test instead of a "semi-infinite line to
         * segment" test, meaning that one must "fake" a semi-infinite line
         * using a big segment. This is done by building a segment originating
         * from the given point and leaving in the opposite direction of the
         * contour center. The full segment can be written in algebraic
         * coordinates as
         *
         * <pre>
         * [PQ] where Q = P + n * CP
         * </pre>
         *
         * , where n is chosen arbitrarily large.
         *
         * @param c a contour
         * @param p a point to test
         * @return true if the point is inside the contour
         */
        public
                double GetDistanceToEdge(Vector2 p)
        {
            Vector2 centre = boundingCircle.GetCentre();

            Vector2 q = p.Subtract(centre).Multiply(10000).Add(p);

            int nb = 0;
            int nbPtsM1 = points.size() - 1;
            double distSq, minDistSq = Double.MAX_VALUE;

            // all points but the last
            for (int index = 0;
                 index < nbPtsM1;
                 ++index)
            {
                Vector2 p1 = points.get(index);
                Vector2 p2 = points.get(index + 1);

                if (Line2D.linesIntersect(p1.x, p1.y, p2.x, p2.y,
                                          p.x, p.y, q.x, q.y))
                {
                    ++nb;
                    distSq = Line2D.ptLineDistSq(p1.x, p1.y, p2.x,
                                                 p2.y, p.x, p.y);
                    if (distSq < minDistSq)
                    {
                        minDistSq = distSq;
                    }
                }
            }

            // last point
            Vector2 p1 = points.get(nbPtsM1);
            Vector2 p2 = points.get(0);
            if (Line2D.linesIntersect(p1.x, p1.y, p2.x, p2.y, p.x,
                                      p.y, q.x, q.y))
            {
                ++nb;
                distSq = Line2D.ptLineDistSq(p1.x, p1.y, p2.x, p2.y, p.x, p.y);
                if (distSq < minDistSq)
                {
                    minDistSq = distSq;
                }
            }

            return ((nb & 1) == 1) ? Math.sqrt(minDistSq) : 0.0;
        }

        @Override
        public
                Iterator<Vector2> iterator()
        {
            return points.iterator();
        }

        void Move(Dims dims, double timeStep)
        {
            Vector2 force = Vector2.Zero();
            double maxDispSq = sampling * timeStep;
            maxDispSq *= maxDispSq;

            int n = points.size();

            if (modelForces == null || modelForces.length != n)
            {
                return;
            }

            for (int index = 0;
                 index < n;
                 ++index)
            {
                Vector2 p = points.get(index);

                if (volumeConstraintForces[index].LengthSq() > 0)
                {
                    p.AddTo(volumeConstraintForces[index]);
                }

                if (modelForces[index] != null)
                {
                    if (p.x < 1 || p.x > dims.width - 2)
                    {
                        modelForces[index].MultiplyBy(0.1);
                    }
                    if (p.y < 1 || p.y > dims.height - 2)
                    {
                        modelForces[index].MultiplyBy(0.1);
                    }
                    force = modelForces[index].Copy();
                }
                else
                {
                    feedbackForces[index].MultiplyBy(0.25);
                }

                // apply feedback forces all the time
                force.AddTo(feedbackForces[index]);

                force.MultiplyBy(timeStep);

                double dispSq = force.LengthSq();

                if (dispSq > maxDispSq)
                {
                    force.MultiplyBy(Math.sqrt(maxDispSq / dispSq));
                }

                p.AddTo(force);

                force.MultiplyBy(0);

                modelForces[index].MultiplyBy(0);
                feedbackForces[index].MultiplyBy(0);
                volumeConstraintForces[index].MultiplyBy(0);
            }

            MakeBoundingCircle();

            // compute some convergence criterion
            if (convergenceWindow == null)
            {
                return;
            }

            convergenceWindow.Add(Math.abs(GetAlgebraicInterior()));
        }

        public
                boolean HasConverged()
        {
            return convergenceWindow.HasConverged();
        }

        private
                int cpt = 0;

        /**
         * Re-samples the Contour according to an 'average distance between
         * points' criterion. This method ensures that the distance between two
         * consecutive points is strictly comprised between a minimum value and
         * a maximum value. In order to avoid oscillatory behavior, 'max' and
         * 'min' should verify the following relations: min < 1, max > 1, 2*min
         * <= max.
         *
         * @param minFactor the minimum distance between two points.
         * @param maxFactor the maximum distance between two points.
         */
        private

                void ReSample(double minFactor, double maxFactor)
                throws ContourException
        {
//            if (GetPerimeter() < 10)
            if (GetPerimeter() < 1)
            {
                throw new ContourException(this, new Contour[0]);
            }

            double minDistanceSq = sampling * minFactor;
            minDistanceSq *= minDistanceSq;
            double maxDistanceSq = sampling * maxFactor;
            maxDistanceSq *= maxDistanceSq;

            // optimization to avoid multiple points.size() calls (WARNING: n must
            // be updated manually whenever points is changed)
            int n = points.size();

            if (contourNormals == null)
            {
                // first pass: update normals once
                contourNormals = new Vector2[n];
                for (int index = 0;
                     index < n;
                     ++index)
                {
                    contourNormals[index] = Vector2.Zero();
                }

                UpdateNormals();
            }

            Contour[] children = (cpt++ & 1) == 0 ? null
                                 : CheckSelfIntersection(sampling);

            if (children != null)
            {
                throw new ContourException(this, children);
            }

            // update the number of total points
            n = points.size();
            boolean noChange = false;

            int iterCount = 0;
            int maxIter = n * 10;

            while (noChange == false)
            {
                // Safeguard
                if (++iterCount > maxIter)
                {
                    System.err.println("[Contours] Warning: hitting "
                                       + "safeguard (preventing infinite resampling)");
                    break;
                }

                noChange = true;

                // all points but the last
                for (int index = 0;
                     index < n - 1;
                     ++index)
                {
                    if (n < 4)
                    {
                        throw new ContourException(this, new Contour[0]);
                    }

                    Vector2 pt1 = points.get(index);
                    Vector2 pt2 = points.get(index + 1);

                    double distanceSq = pt1.DistanceSq(pt2);

                    if (distanceSq < minDistanceSq)
                    {
                        noChange = false;
                        pt2.AddTo(pt1);
                        pt2.MultiplyBy(0.5);

                        points.remove(index--);
                        --n;
                    }
                    else if (distanceSq > maxDistanceSq)
                    {
                        noChange = false;
                        Vector2 pt3 = pt1.Add(pt2);
                        pt3.MultiplyBy(0.5);

                        points.add(++index, pt3);
                        ++n;
                    }
                }

                // last point
                Vector2 pt1 = points.get(n - 1);
                Vector2 pt2 = points.get(0);

                if (pt1.DistanceSq(pt2) < minDistanceSq)
                {
                    noChange = false;
                    pt2.AddTo(pt1);
                    pt2.MultiplyBy(0.5);

                    points.remove(--n);
                }
                else if (pt1.DistanceSq(pt2) > maxDistanceSq)
                {
                    noChange = false;
                    Vector2 pt3 = pt1.Add(pt2);
                    pt3.MultiplyBy(0.5);

                    points.add(pt3);
                    ++n;
                }
            }

            // re-sampling is done => update internal structures
            final
                    int nbPoints = n;
            if (modelForces == null || modelForces.length != nbPoints)
            {
                modelForces = new Vector2[nbPoints];
                contourNormals = new Vector2[nbPoints];
                feedbackForces = new Vector2[nbPoints];
                volumeConstraintForces = new Vector2[nbPoints];

                for (int index = 0;
                     index < nbPoints;
                     ++index)
                {
                    modelForces[index] = Vector2.Zero();
                    contourNormals[index] = Vector2.Zero();
                    feedbackForces[index] = Vector2.Zero();
                    volumeConstraintForces[index] = Vector2.Zero();
                }
            }

            UpdateNormals();
            MakeBoundingCircle();
        }

        private
                void Triangulate(ROIExtractor.ROI roi)// throws TopologyException
        {
            ArrayList<Segment> segments = new ArrayList<>();

            boolean[] mask = roi.MakeMask();
            // erase first line and first row to ensure closed contours
            Arrays.fill(mask, 0, roi.dims.width - 1, false);
            for (int index = 0;
                 index < mask.length;
                 index += roi.dims.width)
            {
                mask[index] = false;
            }

            for (int y = 0, index = 0;
                 y < roi.dims.height;
                 ++y)
            {
                for (int x = 0;
                     x < roi.dims.width;
                     ++x, ++index)
                {
                    // The image is divided into square cells containing two
                    // triangles each:
                    //
                    // a---b---
                    // |../|../
                    // |./.|./.
                    // |/..|/..
                    // c---d---
                    //
                    // By convention I choose to turn around the object in a
                    // clockwise fashion
                    // Warning: to ensure connectivity, the objects must NOT touch
                    // the image border, strange behavior may occur otherwise

                    boolean a = mask[index];
                    boolean b = (x + 1 < roi.dims.width) && mask[index + 1];
                    boolean c = (y + 1 < roi.dims.height) && mask[index + roi.dims.width];
                    boolean d = (x + 1 < roi.dims.width) && (y + 1 < roi.dims.height)
                                && mask[index + 1 + roi.dims.width];

                    // For each triangle, check for difference between image values
                    // to determine the contour location
                    // => there are 6 possible combinations in each triangle, that
                    // is 12 per cube
                    if (a != b)
                    {
                        if (b == c) // diagonal edge
                        {
                            if (a == false) // b,c are inside
                            {
                                CreateEdge(segments, x, y + 0.5, x + 0.5, y);
                            }
                            else
                            // b,c are outside
                            {
                                CreateEdge(segments, x + 0.5, y, x, y + 0.5);
                            }
                        }
                        else // a = c -> vertical edge
                        if (a == false) // a,c are outside
                        {
                            CreateEdge(segments, x + 0.5, y + 0.5, x + 0.5, y);
                        }
                        else
                        // a,c are inside
                        {
                            CreateEdge(segments, x + 0.5, y, x + 0.5, y + 0.5);
                        }
                    }
                    else // a = b -> horizontal edge only if c is different
                    {
                        if (a != c)
                        {
                            if (a == false) // a,b are outside
                            {
                                CreateEdge(segments, x, y + 0.5, x + 0.5, y + 0.5);
                            }
                            else
                            // a,b are inside
                            {
                                CreateEdge(segments, x + 0.5, y + 0.5, x, y + 0.5);
                            }
                        }
                    }

                    if (c != d)
                    {
                        if (b == c) // diagonal edge
                        {
                            if (c == false) // b,c are outside
                            {
                                CreateEdge(segments, x + 0.5, y + 1, x + 1, y + 0.5);
                            }
                            else
                            // b,c are inside
                            {
                                CreateEdge(segments, x + 1, y + 0.5, x + 0.5, y + 1);
                            }
                        }
                        else // b = d -> vertical edge
                        if (c == false) // b,d are inside
                        {
                            CreateEdge(segments, x + 0.5, y + 1, x + 0.5, y + 0.5);
                        }
                        else
                        // b,d are outside
                        {
                            CreateEdge(segments, x + 0.5, y + 0.5, x + 0.5, y + 1);
                        }
                    }
                    else // c = d -> horizontal edge only if b is different
                    {
                        if (b != c)
                        {
                            if (b == false) // c,d are inside
                            {
                                CreateEdge(segments, x + 0.5, y + 0.5, x + 1, y + 0.5);
                            }
                            else
                            // c,d are outside
                            {
                                CreateEdge(segments, x + 1, y + 0.5, x + 0.5, y + 0.5);
                            }
                        }
                    }
                }
            }

            if (segments.isEmpty())
            {
                return;
            }

            for (Vector2 p : segments.get(0))
            {
                p.AddTo(roi.origin.ToVector2());
                AddPoint(p);
            }

            // at this point the triangulated contour has an actual resolution of halfgrid
            // if 2*resolution < desired_resolution, resample() will loop and destroy the contour
            // decimate the contour by a factor 2 recursively until 2*resolution >= desired_resolution
            double current_resolution_doubled = 1;
            while (current_resolution_doubled < sampling * 0.7)
            {
                for (int index = 0;
                     index < points.size();
                     ++index)
                {
                    points.remove(index);
                }
                current_resolution_doubled *= 2;
            }
        }

        private

                void UpdateNormals()
        {
            int n = points.size();

            clockWise = GetAlgebraicInterior() < 0;

            // first point
            {
                Vector2 p1 = points.get(n - 1);
                Vector2 p2 = points.get(1);
                contourNormals[0] = p1.Subtract(p2);
                contourNormals[0].RotatateByMinus90By();
                contourNormals[0].Normalise();

                if (clockWise)
                {
                    contourNormals[0].MultiplyBy(-1);
                }
            }

            // middle points
            for (int index = 1;
                 index < n - 1;
                 ++index)
            {
                Vector2 p1 = points.get(index - 1);
                Vector2 p2 = points.get(index + 1);
                contourNormals[index] = p1.Subtract(p2);
                contourNormals[index].RotatateByMinus90By();
                contourNormals[index].Normalise();

                if (clockWise)
                {
                    contourNormals[index].MultiplyBy(-1);
                }
            }

            // last point
            {
                Vector2 p1 = points.get(n - 2);
                Vector2 p2 = points.get(0);
                contourNormals[n - 1] = p1.Subtract(p2);
                contourNormals[n - 1].RotatateByMinus90By();
                contourNormals[n - 1].Normalise();

                if (clockWise)
                {
                    contourNormals[n - 1].MultiplyBy(-1);
                }
            }
        }

        ROIExtractor.ROI ToROI(Dims imageDims)
        {
            // NOTE: need to check this works, and can be optimised?
            int minX = Integer.MAX_VALUE;
            int minY = Integer.MAX_VALUE;
            int maxX = 0;
            int maxY = 0;

            for (Vector2 point : points)
            {
                Vector2i floor = point.ToVector2i_floor();

                if (floor.x > maxX)
                {
                    maxX = floor.x;
                }
                if (floor.y > maxY)
                {
                    maxY = floor.y;
                }
                if (floor.x < minX)
                {
                    minX = floor.x;
                }
                if (floor.y < minY)
                {
                    minY = floor.y;
                }
            }

            minX = Math.max(minX, 0);
            minY = Math.max(minY, 0);
            maxX = Math.min(maxX, imageDims.width);
            maxY = Math.min(maxY, imageDims.height);

            Vector2i origin = new Vector2i(minX, minY);
            Vector2i dimVec = new Vector2i(maxX, maxY).Subtract(origin).Add(Vector2i.One());
            Dims dims = new Dims(dimVec.x, dimVec.y);

            if (dims.width <= 0 || dims.height <= 0)
            {
                return null;
            }

            byte[] pixels = new byte[dims.width * dims.height];

            for (int y = 0, index = 0;
                 y < dims.height;
                 ++y)
            {
                for (int x = 0;
                     x < dims.width;
                     ++x, ++index)
                {
                    for (double addVal = 0;
                         addVal <= 1.0;
                         addVal += 0.1)
                    {
                        Vector2 testPoint = new Vector2i(x, y).Add(origin).
                                ToVector2().Add(addVal);

                        if (GetDistanceToEdge(testPoint) > 0)
                        {
                            pixels[index] = (byte) 1;
                            break;
                        }
                        else
                        {
                            pixels[index] = (byte) 0;
                        }
                    }
                }
            }

            ROIExtractor.ROI[] rois = ROIExtractor.SingleImageROIs(pixels, dims);

            if (rois.length > 1)
            {
                System.err.println("Converting contour at "
                                   + GetBoundingCircle().GetCentre().toString() + " to an ROI resulted "
                                   + "in a split");

                return null;
            }
            else if (rois.length == 0)
            {
                return null;
            }

            ROIExtractor.ROI roi = rois[0];

            roi.origin.AddTo(origin);

            return roi;
        }

        boolean DoesIntersect(Contour other)
        {
            return GetBoundingCircle().DoesIntersect(other.GetBoundingCircle());
        }

        // Copied from Icy plugins.adufour.activecontours
        public
                class ContourException extends Exception
        {

            public final
                    Contour source;

            public final
                    Contour[] children;

            /**
             * Creates a new exception for the specified contour
             *
             * @param contour the contour undergoing a topology break
             * @param children an array containing zero or more contours that
             * should replace the contour raising the exception
             */
            public
                    ContourException(Contour contour, Contour[] children)
            {
                super("Topology break detected in contour " + contour.hashCode());
                this.source = contour;
                this.children = children;
            }

        }
    }

    private final
            double EPSILON = 0.0000001;

    private
            ROIExtractor.ROI[] roiOutput = new ROIExtractor.ROI[0];

    private final
            TrackGroup trackGroup = new TrackGroup();

    /**
     * All contours present on the current time point
     */
    private final
            HashSet<ContourDetection> allContoursAtTimeT = new HashSet<>();

    /**
     * Set of contours that have not yet converged on the current time point
     */
    private final
            HashSet<ContourDetection> evolvingContoursAtTimeT = new HashSet<>();

    private
            Processor multiThreadService;

    private final
            HashMap<TrackSegment, Double> volumes = new HashMap<>();

    private
            int currentTime = 0;

    private final
            Dims dims;

    ///////////////////////  PARAMS START  ///////////////////////
    private final
            BinaryMask.FilterKernel erodeFilterKernel;
    private final
            int nErodeFilterPasses;

    private final
            double contourSampling;
    private final
            double contourDivisionSensitivity;
    private final
            double contourConvCritSq;
    private final
            int contourSlidingWindowSize;
    private final
            double contourTimeStep;
    private final
            int contourConvergence_nIter;

    private final
            double internalWeight;
    private final
            double edgeWeight;
    private final
            double regionWeight;
    private final
            double regionSensitivity;
    private final
            double axisWeight;
    private final
            double balloonWeight;
    private final
            double volumeConstraint;
    private final
            boolean couplingFlag;
    ///////////////////////  PARAMS  END   ///////////////////////

    public static
            class ContourDetection extends Tracking.Detection
    {

        private final
                Contour contour;

        public
                ContourDetection(Contour c, int t)
        {
            super(c.GetBoundingCircle().GetCentre(), t);
            contour = c;
        }

        public
                ContourDetection(ContourDetection c) throws CloneNotSupportedException
        {
            super(null, c.t);
            c.UpdatePosition();
            position = c.position.Copy();
            this.contour = c.contour.clone();
        }

        public
                ContourDetection(ContourDetection c, double convWinExpandFactor)
        {
            super(null, c.t);
            c.UpdatePosition();
            position = c.position.Copy();
            this.contour = new Contour(c.contour, convWinExpandFactor);
        }

        private
                void UpdatePosition()
        {
            position = contour.GetBoundingCircle().GetCentre();
        }

        Contour GetContour()
        {
            return contour;
        }

        @Override
        protected
                ContourDetection clone() throws CloneNotSupportedException
        {
            super.clone();
            return new ContourDetection(this);
        }
    }

    private final
            ContourUpdate updater;

    private
            Future<?> updaterRefreshJob;

    ContourTracker(Processor p, Contour[] contours, Dims dims,
                   final
                   double contourSampling, final
                   double contourDivisionSensitivity,
                   final
                   double contourConvCritSq, final
                   int contourSlidingWindowSize,
                   final
                   double contourTimeStep, final
                   int contourConvergence_nIter,
                   final
                   double internalWeight, final
                   double edgeWeight, final
                   double regionWeight,
                   final
                   double regionSensitivity, final
                   double axisWeight, final
                   double balloonWeight,
                   final
                   double volumeConstraint, final
                   boolean couplingFlag, final
                   BinaryMask.FilterKernel erodeFilterKernel, final
                   int nErodeFilterPasses)
    {
        this(p, contours, dims, null, contourSampling, contourDivisionSensitivity,
             contourConvCritSq, contourSlidingWindowSize, contourTimeStep,
             contourConvergence_nIter, internalWeight, edgeWeight,
             regionWeight, regionSensitivity, axisWeight, balloonWeight,
             volumeConstraint, couplingFlag, erodeFilterKernel, nErodeFilterPasses);
    }

    ContourTracker(Processor p, Contour[] contours, Dims dims, ContourUpdate updater,
                   final
                   double contourSampling, final
                   double contourDivisionSensitivity,
                   final
                   double contourConvCritSq, final
                   int contourSlidingWindowSize,
                   final
                   double contourTimeStep, final
                   int contourConvergence_nIter,
                   final
                   double internalWeight, final
                   double edgeWeight, final
                   double regionWeight,
                   final
                   double regionSensitivity, final
                   double axisWeight, final
                   double balloonWeight,
                   final
                   double volumeConstraint, final
                   boolean couplingFlag, final
                   BinaryMask.FilterKernel erodeFilterKernel, final
                   int nErodeFilterPasses)
    {
        ///////////////////////  PARAMS START  ///////////////////////
        this.nErodeFilterPasses = nErodeFilterPasses;
        this.erodeFilterKernel = erodeFilterKernel;
        this.contourSampling = contourSampling;
        this.contourDivisionSensitivity = contourDivisionSensitivity;
        this.contourConvCritSq = contourConvCritSq;
        this.contourSlidingWindowSize = contourSlidingWindowSize;
        this.contourTimeStep = contourTimeStep;
        this.contourConvergence_nIter = contourConvergence_nIter;
        this.internalWeight = internalWeight;
        this.edgeWeight = edgeWeight;
        this.regionWeight = regionWeight;
        this.regionSensitivity = regionSensitivity;
        this.axisWeight = axisWeight;
        this.balloonWeight = balloonWeight;
        this.volumeConstraint = volumeConstraint;
        this.couplingFlag = couplingFlag;
        ///////////////////////  PARAMS  END   ///////////////////////

        multiThreadService = p;
        this.dims = dims;
        this.updater = updater;
        updaterRefreshJob = null;

        Future<?>[] tasks = new Future<?>[contours.length];
        int index = -1;

        for (Contour con_in : contours)
        {

            Runnable initializer = ()
                    ->
            {
                ContourDetection contour_det = new ContourDetection(
                        new Contour(con_in), 0);

                TrackSegment segment = new TrackSegment();
                segment.AddDetection(contour_det);
                synchronized (trackGroup)
                {
                    trackGroup.AddTrackSegment(segment);
                }

                AddContourToUpdater(contour_det.GetContour());
            };

            tasks[++index] = multiThreadService.submit(initializer);
        }

        try
        {
            for (Future<?> future : tasks)
            {
                future.get();
            }
        }
        catch (InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }

        RefreshUpdater();
    }

    private
            void AddContourToUpdater(Contour contour)
    {
        if (updater != null)
        {
            synchronized (updater)
            {
                updater.AddContour(contour);
            }
        }
    }

    private
            void UpdateContourInUpdater(Contour contour)
    {
        if (updater != null)
        {
            synchronized (updater)
            {
                updater.UpdateContour(contour);
            }
        }
    }

    private
            void RemoveContourFromUpdater(Contour contour)
    {
        if (updater != null)
        {
            synchronized (updater)
            {
                updater.RemoveContour(contour);
            }
        }
    }

    private
            void RefreshUpdater()
    {
        if (updater != null)
        {
            if (updaterRefreshJob != null)
            {
                try
                {
                    updaterRefreshJob.get();
                }
                catch (InterruptedException | ExecutionException ex)
                {
                    Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            updaterRefreshJob = multiThreadService.submit(()
                    ->
            {
                synchronized (updater)
                {
                    updater.Refresh();
                }
            });
        }
    }

    public
            TrackGroup GetTrackGroup()
    {
        return trackGroup;
    }

    void EvolveContoursToNextFrame(final
            byte[] regionPixels, final
                                   byte[] edgePixels, String header)
    {
        final
                int t = currentTime++;

        InitContours(t);

        // evolve contours on the current image
        EvolveContours(regionPixels, edgePixels, t, header);

        FindNewObjects(regionPixels, t);

        EvolveContours(regionPixels, edgePixels, t, header);

        ClearDeadContours(regionPixels, t, header);

        // store detections and results
        StoreResults(t);

        // clean other non-necessary stuff
        allContoursAtTimeT.clear();
        ClearPreviousTimepoint(t);
    }

    private
            void ClearDeadContours(final
                    byte[] pixels, final
                                   int t, String header)
    {
        final
                AtomicInteger nKilled = new AtomicInteger();
        LinkedList<Future<?>> tasks = new LinkedList<>();

        new ArrayList<>(trackGroup.GetTrackSegmentList()).
                stream().map((segment) -> (Runnable) ()
                ->
        {
            if (segment == null)
            {
                synchronized (trackGroup)
                {
                    trackGroup.RemoveTrackSegment(segment);
                }
                synchronized (nKilled)
                {
                    nKilled.incrementAndGet();
                }
            }
        }).forEach((task)
                ->
        {
            tasks.add(multiThreadService.submit(task));
        });

        try
        {
            for (Future<?> future : tasks)
            {
                future.get();
            }
        }
        catch (InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }

        tasks.clear();
        new ArrayList<>(trackGroup.GetTrackSegmentList()).stream().map((segment)
                -> (Runnable) ()
                ->
        {
            Detection detection = segment.GetDetectionAtTime(t);

            if (detection != null)
            {
                Contour contour = ((ContourDetection) detection)
                        .GetContour();

                Vector2 centre = contour.GetBoundingCircle()
                        .GetCentre();
                double radius = contour.GetBoundingCircle()
                        .GetRadius();

                if (centre.x < -radius || centre.x > (dims.width
                                                      + radius) || centre.y < -radius || centre.y
                                                                                         > (dims.height + radius))
                {
                    segment.RemoveDetection(detection);
                    synchronized (nKilled)
                    {
                        nKilled.incrementAndGet();
                    }

                    RemoveContourFromUpdater(contour);

                }
                else
                {
                    int total = 0;
                    ROIExtractor.ROI roi = contour.ToROI(dims);

                    if (roi == null)
                    {
                        segment.RemoveDetection(detection);
                        synchronized (nKilled)
                        {
                            nKilled.incrementAndGet();
                        }

                        RemoveContourFromUpdater(contour);

                    }
                    else
                    {
                        boolean[] mask = roi.MakeMask();

                        for (int y = 0, index = 0;
                             y < roi.dims.height;
                             ++y)
                        {
                            for (int x = 0;
                                 x < roi.dims.width;
                                 ++x, ++index)
                            {
                                if (mask[index])
                                {
                                    int globalY = y + roi.origin.y;
                                    int globalX = x + roi.origin.x;

                                    if (globalY >= 0
                                        && globalY < dims.height
                                        && globalX >= 0
                                        && globalX < dims.width)
                                    {
                                        int globalIndex = (globalY
                                                           * dims.width)
                                                          + globalX;
                                        total += (pixels[globalIndex] & 0xff);
                                    }
                                }
                            }
                        }

                        if (total == 0)
                        {
                            segment.RemoveDetection(detection);
                            synchronized (nKilled)
                            {
                                nKilled.incrementAndGet();
                            }

                            RemoveContourFromUpdater(contour);
                        }
                    }
                }

                if (segment.GetDetectionList().isEmpty())
                {
                    synchronized (trackGroup)
                    {
                        trackGroup.RemoveTrackSegment(segment);
                    }
                }
            }
        }).forEach((task)
                ->
        {
            tasks.add(multiThreadService.submit(task));
        });

        try
        {
            for (Future<?> future : tasks)
            {
                future.get();
            }
        }
        catch (InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println(header + "Removed " + nKilled + " dead contours "
                           + "at timepoint " + t);

        RefreshUpdater();
    }

    private
            void ClearPreviousTimepoint(int currentT)
    {
        if (currentT > 0)
        {
            final
                    int t = currentT - 1;

            LinkedList<Future<?>> tasks = new LinkedList<>();

            new ArrayList<>(trackGroup.GetTrackSegmentList()).stream().map((segment)
                    -> (Runnable) ()
                    ->
            {
                Detection detection = segment.GetDetectionAtTime(t);

                if (detection != null)
                {
                    segment.RemoveDetection(detection);

                    if (segment.GetDetectionList().isEmpty())
                    {
                        synchronized (trackGroup)
                        {
                            trackGroup.RemoveTrackSegment(segment);
                        }
                    }
                }
            }).forEach((task)
                    ->
            {
                tasks.add(multiThreadService.submit(task));
            });

            try
            {
                for (Future<?> future : tasks)
                {
                    future.get();
                }
            }
            catch (InterruptedException | ExecutionException ex)
            {
                Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    private
            void InitContours(final
                    int t)
    {
        if (t > 0)
        {
            Future<?>[] tasks = new Future<?>[trackGroup.
                    GetTrackSegmentList().size()];
            int index = -1;

            for (TrackSegment segment : trackGroup.GetTrackSegmentList())
            {
                Runnable task = ()
                        ->
                {
                    ContourDetection previous = (ContourDetection) segment.GetDetectionAtTime(t - 1);

                    if (previous != null)
                    {
//                                ContourDetection next = new ContourDetection(previous, 2);
                        ContourDetection next = new ContourDetection(previous, 1.0);
                        next.t = t;

                        segment.AddDetection(next);
                    }
                };

                tasks[++index] = multiThreadService.submit(task);
            }

            try
            {
                for (Future<?> future : tasks)
                {
                    future.get();
                }
            }
            catch (InterruptedException | ExecutionException ex)
            {
                Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    private
            void EvolveContours(final
                    byte[] regionPixels, final
                                byte[] edgePixels, final
                                int t, String header)
    {
        // retrieve the contours on the current frame and store them in currentContours

        allContoursAtTimeT.clear();

        trackGroup.GetTrackSegmentList().stream().map((segment)
                -> segment.GetDetectionAtTime(t)).filter((det)
                -> (det != null)).forEach((det)
                ->
        {
            allContoursAtTimeT.add((ContourDetection) det);
        });

        if (allContoursAtTimeT.isEmpty())
        {
            return;
        }

        int iter = 0;
        int nbConvergedContours = 0;

        int nPrevious_resample = -1;
        int nPrevious_deform = -1;

        while (nbConvergedContours < allContoursAtTimeT.size())
        {
            nbConvergedContours = 0;

            // take a snapshot of the current list of evolving (i.e. non-converged) contours
            evolvingContoursAtTimeT.clear();

            for (ContourDetection contour : allContoursAtTimeT)
            {
                if (contour.GetContour().HasConverged())
                {
                    ++nbConvergedContours;
                    continue;
                }

                // if the contour hasn't converged yet, store it for the main loop
                evolvingContoursAtTimeT.add(contour);
            }

            if (evolvingContoursAtTimeT.isEmpty())
            {
                break;
            }

            // re-sample the contours to ensure homogeneous resolution
            nPrevious_resample = ResampleContours(t, nPrevious_resample, header);

            // compute deformations issued from the energy minimization
            nPrevious_deform = DeformContours(regionPixels, edgePixels, t,
                                              nPrevious_deform, header);

            // compute energy
            // computeEnergy(mainService, allContours);
            if (iter > contourConvergence_nIter)
            {
                System.out.println("Stopped converging on frame "
                                   + "" + t + " after " + iter + " iterations");
                return;
            }

            ++iter;
        }

        System.out.println("Converged on frame "
                           + t + " in " + iter + " iterations");
    }

    /**
     * Deform contours together
     *
     * @param field the allowed displacement field for all contours. By default
     * this field is the image volume, but this could be tweaked to restrict
     * evolution to a certain area
     */
    private
            int DeformContours(final
                    byte[] regionPixels, final
                               byte[] edgePixels,
                               final
                               int t, int nPrevious, String header)
    {
        int n = evolvingContoursAtTimeT.size();

        if (n != nPrevious)
        {
            System.out.println(header + "Deforming " + n + " contours at "
                               + "timepoint " + t);
        }

        if (evolvingContoursAtTimeT.size() == 1 && allContoursAtTimeT.size() == 1)
        {
            // no multi-threading needed

            ContourDetection contour = evolvingContoursAtTimeT.iterator().next();
            TrackSegment segment = trackGroup.GetTrackSegmentWithDetection(contour);

            if (Math.abs(edgeWeight) > EPSILON)
            {
                contour.GetContour().ComputeEdgeForces(edgePixels, dims, edgeWeight);
            }

            if (internalWeight > EPSILON)
            {
                contour.GetContour().ComputeInternalForces(internalWeight);
            }

            if (regionWeight > EPSILON)
            {
                contour.GetContour().ComputeRegionForces(regionPixels, dims,
                                                         regionWeight, regionSensitivity);
            }

            if (Math.abs(balloonWeight) > EPSILON)
            {
                contour.GetContour().ComputeBalloonForces(balloonWeight);
            }

            if (axisWeight > EPSILON)
            {
                contour.GetContour().ComputeAxisForces(axisWeight);
            }

            if (volumeConstraint > EPSILON && volumes.containsKey(segment))
            {
                contour.GetContour().ComputeVolumeConstraint(volumes.get(segment));
            }

            contour.GetContour().Move(dims, contourTimeStep);

            UpdateContourInUpdater(contour.GetContour());
        }
        else
        {
            ArrayList<Callable<ContourDetection>> tasks
                                                  = new ArrayList<>(evolvingContoursAtTimeT.size());

            for (final
                    ContourDetection contour : evolvingContoursAtTimeT)
            {
                final
                        TrackSegment segment = trackGroup.
                                GetTrackSegmentWithDetection(contour);

                tasks.add((Callable<ContourDetection>) ()
                        ->
                {
                    if (Math.abs(edgeWeight) > EPSILON)
                    {
                        contour.GetContour().ComputeEdgeForces(edgePixels,
                                                               dims, edgeWeight);
                    }

                    if (internalWeight > EPSILON)
                    {
                        contour.GetContour().ComputeInternalForces(internalWeight);
                    }

                    if (regionWeight > EPSILON)
                    {
                        contour.GetContour().ComputeRegionForces(regionPixels,
                                                                 dims, regionWeight, regionSensitivity);
                    }

                    if (Math.abs(balloonWeight) > EPSILON)
                    {
                        contour.GetContour().ComputeBalloonForces(balloonWeight);
                    }

                    if (axisWeight > EPSILON)
                    {
                        contour.GetContour().ComputeAxisForces(axisWeight);
                    }

                    if (couplingFlag)
                    {
                        // Don't move the contours just now: coupling feedback must be computed
                        // against ALL contours (including those which have already converged)
                        allContoursAtTimeT.stream().filter((otherContour)
                                -> !(otherContour == null || otherContour == contour))
                                .forEach((otherContour)
                                        ->
                                {
                                    contour.GetContour().ComputeFeedbackForces(otherContour.GetContour());
                                });

                        if (volumeConstraint > EPSILON && volumes.containsKey(segment))
                        {
                            contour.GetContour().ComputeVolumeConstraint(volumes.get(segment));
                        }
                    }
                    else
                    {
                        // move contours asynchronously
                        contour.GetContour().Move(dims, contourTimeStep);

                        UpdateContourInUpdater(contour.GetContour());
                    }

                    return contour;
                });
            }

            try
            {
                for (Future<ContourDetection> future
                     : multiThreadService.invokeAll(tasks))
                {
                    try
                    {
                        future.get();
                    }
                    catch (ExecutionException e)
                    {
                        throw new RuntimeException(e.getCause());
                    }
                }
            }
            catch (InterruptedException | RejectedExecutionException e)
            {
                Thread.currentThread().interrupt();
                return n;
            }

            if (couplingFlag)
            {
                tasks.clear();

                for (final
                        ContourDetection contour : evolvingContoursAtTimeT)
                {
                    tasks.add((Callable<ContourDetection>) ()
                            ->
                    {
                        contour.GetContour().Move(dims, contourTimeStep);

                        UpdateContourInUpdater(contour.GetContour());

                        return contour;
                    });
                }

                try
                {
                    for (Future<ContourDetection> future
                         : multiThreadService.invokeAll(tasks))
                    {
                        try
                        {
                            future.get();
                        }
                        catch (ExecutionException e)
                        {
                            throw new RuntimeException(e.getCause());
                        }
                    }
                }
                catch (InterruptedException | RejectedExecutionException e)
                {
                    Thread.currentThread().interrupt();
                }
            }
        }

        RefreshUpdater();

        return n;
    }

    private
            int ResampleContours(final
                    int t, int nPrevious, String header)
    {
        int n = evolvingContoursAtTimeT.size();

        if (n != nPrevious)
        {
            System.out.println(header + "Resampling " + n + " contours "
                               + "at timepoint " + t);
        }

        boolean loop = true;

        int maxIterations = 10000;
        int itCount = 0;

        while (loop)
        {
            if (itCount++ > maxIterations || Thread.currentThread().
                    isInterrupted())
            {
                break;
            }

            loop = false;

            if (evolvingContoursAtTimeT.size() == 1)
            {
                // no multi-threading needed
                Iterator<ContourDetection> iterator = evolvingContoursAtTimeT
                        .iterator();
                if (iterator.hasNext())
                {
                    ContourDetection contour = iterator.next();
                    ReSampler reSampler = new ReSampler(trackGroup, contour,
                                                        evolvingContoursAtTimeT, allContoursAtTimeT,
                                                        updater);
                    if (reSampler.call())
                    {
                        loop = true;
                    }
                }
            }
            else
            {
                ArrayList<ReSampler> tasks = new ArrayList<>(
                        evolvingContoursAtTimeT.size());

                evolvingContoursAtTimeT.stream().forEach((contour)
                        ->
                {
                    tasks.add(new ReSampler(trackGroup, contour,
                                            evolvingContoursAtTimeT, allContoursAtTimeT,
                                            updater));
                });

                try
                {
                    for (Future<Boolean> resampled : multiThreadService.
                            invokeAll(tasks))
                    {
                        if (resampled.get())
                        {
                            loop = true;
                        }
                    }
                }
                catch (InterruptedException e)
                {
                    Thread.currentThread().interrupt();
                }
                catch (ExecutionException e)
                {
                    throw new RuntimeException(e);
                }
                catch (RuntimeException e)
                {
                    throw e;
                }
            }

            RefreshUpdater();
        }

        return n;
    }

    private
            void FindNewObjects(final
                    byte[] pixels_in, final
                                int t)
    {
        byte[] pixels = new byte[pixels_in.length];
        System.arraycopy(pixels_in, 0, pixels, 0, pixels.length);

        LinkedList<Future<?>> tasks = new LinkedList<>();

        trackGroup.GetTrackSegmentList().stream().map((segment)
                -> (ContourDetection) (segment.GetDetectionAtTime(t))).
                filter((contour) -> !(contour == null)).map((contour)
                -> (Runnable) ()
                ->
        {
            ROIExtractor.ROI roi = contour.GetContour().ToROI(dims);

            if (roi != null)
            {
                boolean[] mask = roi.MakeMask();

                for (int y = 0, index = 0;
                     y < roi.dims.height;
                     ++y)
                {
                    for (int x = 0;
                         x < roi.dims.width;
                         ++x, ++index)
                    {
                        if (mask[index])
                        {
                            int globalY = y + roi.origin.y;
                            int globalX = x + roi.origin.x;

                            if (globalY >= 0
                                && globalY < dims.height
                                && globalX >= 0
                                && globalX < dims.width)
                            {
                                int globalIndex = (globalY
                                                   * dims.width)
                                                  + globalX;
                                pixels[globalIndex] = (byte) 0;
                            }
                        }
                    }
                }
            }

        }).
                forEach((task)
                        ->
                {
                    tasks.add(multiThreadService.submit(task));
                });

        try
        {
            for (Future<?> future : tasks)
            {
                future.get();
            }
        }
        catch (InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }

        for (int i = 0;
             i < nErodeFilterPasses;
             ++i)
        {
            BinaryMask.ErodeFilter(pixels, dims, erodeFilterKernel);
        }

        ROIExtractor.ROI[] rois = ROIExtractor.SingleImageROIs(pixels, dims);

        tasks.clear();
        ArrayList<Contour> newContours = new ArrayList<>(rois.length);

        for (ROIExtractor.ROI roi : rois)
        {
            Runnable task = ()
                    ->
            {
                Contour newContour = new Contour(contourSampling, contourDivisionSensitivity,
                                                 contourSlidingWindowSize, contourConvCritSq, roi);

                synchronized (newContours)
                {
                    newContours.add(newContour);
                }
            };

            tasks.add(multiThreadService.submit(task));
        }

        tasks.stream().forEach((future)
                ->
        {
            try
            {
                future.get();
            }
            catch (InterruptedException | ExecutionException ex)
            {
                System.err.println(ex.getMessage());
            }
        });

        tasks.clear();
        HashSet<Integer> toRem = new HashSet<>();

        trackGroup.GetTrackSegmentList().stream().map((segment)
                -> (ContourDetection) (segment.GetDetectionAtTime(t))).
                filter((contour) -> !(contour == null)).map((contour)
                -> (Runnable) ()
                ->
        {
            int index = -1;
            for (Contour newCon : newContours)
            {
                ++index;

                if (contour.GetContour().DoesIntersect(newCon))
                {
                    synchronized (toRem)
                    {
                        toRem.add(index);
                    }
                }
            }

        }).forEach((task)
                ->
        {
            tasks.add(multiThreadService.submit(task));
        });

        try
        {
            for (Future<?> future : tasks)
            {
                future.get();
            }
        }
        catch (InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }

        Integer[] toRem_array = toRem.toArray(new Integer[0]);
        Arrays.sort(toRem_array, Collections.reverseOrder());

        for (int index : toRem_array)
        {
            newContours.remove(index);
        }

        System.out.println("Found " + newContours.size() + " new contours in "
                           + "frame " + t);

        for (int index = 0;
             index < newContours.size();
             ++index)
        {
            for (TrackSegment segment : trackGroup.GetTrackSegmentList())
            {
                Detection detection = segment.GetLastDetection();

                if (segment.nextList.isEmpty() && detection.t == (t - 1) && 
                    ((ContourDetection) detection).GetContour().
                            DoesIntersect(newContours.get(index)))
                {   
                    segment.AddDetection(new ContourDetection(newContours.
                            get(index), t));

                    AddContourToUpdater(newContours.
                            get(index));

                    newContours.remove(index--);
                    break;
                }
            }
        }

        newContours.stream().map((newContour)
                ->
        {
            AddContourToUpdater(newContour);

            TrackSegment newSegment = new TrackSegment();
            newSegment.AddDetection(new ContourDetection(newContour, t));
            return newSegment;
        }).forEach((newSegment)
                ->
        {
            trackGroup.AddTrackSegment(newSegment);
        });

        RefreshUpdater();
    }

    private
            void StoreResults(final
                    int t)
    {
        ArrayList<TrackSegment> segments = trackGroup.GetTrackSegmentList();

        ArrayList<ROIExtractor.ROI> rois = new ArrayList<>(roiOutput.length);

        Future<?>[] tasks = new Future<?>[segments.size()];

        for (int index = 0;
             index < segments.size();
             ++index)
        {
            TrackSegment segment = segments.get(index);

            ContourDetection contour = (ContourDetection) segment.
                    GetDetectionAtTime(t);

            Runnable task = ()
                    ->
            {
                if (contour != null)
                {

                    if (!volumes.containsKey(segment))
                    {
                        double vol = Math.abs(contour.GetContour().
                                GetAlgebraicInterior());

                        synchronized (volumes)
                        {
                            if (!volumes.containsKey(segment))
                            {
                                volumes.put(segment, vol);
                            }
                        }
                    }

                    // output as ROIs
                    ROIExtractor.ROI roi = contour.GetContour().ToROI(dims);
                    if (roi != null)
                    {
                        roi.t = t;
                        roi.trackID = segment.GetId();

                        if (segment.previousList.isEmpty() == false)
                        {
                            roi.previousTrackID = segment.previousList.get(0).GetId();
                        }

                        roi.area = Math.abs(contour.GetContour().
                                GetAlgebraicInterior_OnlyPositive(dims));
                        roi.perimeter = contour.GetContour().
                                GetPerimeter_OnlyPositive(dims);
                        roi.centre = contour.GetContour().
                                GetBoundingCircle().GetCentre_OnlyPositive(dims);

                        synchronized (rois)
                        {
                            rois.add(roi);
                        }
                    }
                }
            };

            tasks[index] = multiThreadService.submit(task);
        }

        for (Future<?> future : tasks)
        {
            try
            {
                future.get();
            }
            catch (InterruptedException | ExecutionException ex)
            {
                System.err.println(ex.getMessage());
            }
        }

        if (rois.size() > 0)
        {
            roiOutput = rois.toArray(new ROIExtractor.ROI[0]);
        }
    }

    ROIExtractor.ROI[] GetROIs()
    {
        return roiOutput;
    }

//    Contour[] GetCurrentContours()
//    {
//        final
//                int t = currentTime - 1;
//
//        LinkedList<Contour> contours = new LinkedList<>();
//
//        trackGroup.GetTrackSegmentList().stream().map((segment)
//                -> segment.GetDetectionAtTime(t)).filter((detection)
//                -> (detection != null)).forEach((detection)
//                ->
//                {
//                    contours.add(((ContourDetection) detection).GetContour());
//        });
//
//        return contours.toArray(new Contour[0]);
//    }
// Copied from Icy plugins.adufour.activecontours
    private static
            class ReSampler implements Callable<Boolean>
    {

        private final
                TrackGroup trackGroup;

        private final
                ContourDetection contour;

        private final
                HashSet<ContourDetection> allContours;

        private final
                HashSet<ContourDetection> evolvingContours;

        private final
                ContourUpdate updater;

        ReSampler(TrackGroup trackGroup, ContourDetection contour,
                  HashSet<ContourDetection> evolvingContours,
                  HashSet<ContourDetection> allContours,
                  ContourUpdate updater)
        {
            this.trackGroup = trackGroup;
            this.contour = contour;
            this.allContours = allContours;
            this.evolvingContours = evolvingContours;
            this.updater = updater;
        }

        @Override
        public
                Boolean call()
        {
            boolean change = false;

            try
            {
                contour.GetContour().ReSample(0.6, 1.4);

                if (updater != null)
                {
                    synchronized (updater)
                    {
                        updater.UpdateContour(contour.GetContour());
                    }
                }
            }
            catch (ContourException e)
            {
                change = true;

                synchronized (allContours)
                {
                    allContours.remove(contour);
                }
                synchronized (evolvingContours)
                {
                    evolvingContours.remove(contour);
                }

                if (updater != null)
                {
                    synchronized (updater)
                    {
                        updater.RemoveContour(contour.GetContour());
                    }
                }

                TrackSegment currentSegment = null;

                for (TrackSegment segment : new ArrayList<>(trackGroup.
                        GetTrackSegmentList()))
                {
                    if (segment == null)
                    {
                        continue;
                    }

                    if (segment.ContainsDetection(contour))
                    {
                        currentSegment = segment;
                        break;
                    }
                }

                if (currentSegment != null)
                {
                    currentSegment.RemoveDetection(contour);

                    if (currentSegment.GetDetectionList().isEmpty())
                    {
                        // the current contour is the only detection in this segment
                        // => remove the whole segment
                        synchronized (trackGroup)
                        {
                            trackGroup.RemoveTrackSegment(currentSegment);
                        }

                        currentSegment = null;
                    }
                }

                if (!(e instanceof ContourException))
                {
                    Logger.getLogger(ContourTracker.class.getName()).log(Level.SEVERE, null, e);
                    return change;
                }

                // 3) Deal with the children
                Contour[] children = ((ContourException) e).children;

                if (children == null)
                {
                    return change;
                }

                for (Contour child_c : children)
                {
                    ContourDetection child = new ContourDetection(child_c,
                                                                  contour.t);
                    synchronized (allContours)
                    {
                        allContours.add(child);
                    }

                    synchronized (evolvingContours)
                    {
                        evolvingContours.add(child);
                    }

                    if (updater != null)
                    {
                        synchronized (updater)
                        {
                            updater.AddContour(child.GetContour());
                        }
                    }

                    // create the new track segment with the child contour
                    TrackSegment childSegment = new TrackSegment();
                    childSegment.AddDetection(child);

                    synchronized (trackGroup)
                    {
                        trackGroup.AddTrackSegment(childSegment);
                    }

                    if (currentSegment != null)
                    {
                        currentSegment.AddNext(childSegment);
                    }
                }
            }

            return change;
        }
    }

}
